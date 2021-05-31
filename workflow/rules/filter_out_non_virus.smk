rule virsorter2_create_db:
    output:
        directory("results/virsorter2/db")
    conda:
        "../envs/virsorter2.yaml"
    threads: threads
    log:
        "results/logs/virsorter2/db.log"
    shell:
        "virsorter setup -d {output} -j {threads} 2> {log}"

rule virsorter2_filter:
    input: 
        get_fasta,
        "results/virsorter2/db"
    output:
        "results/virsorter2_filter/{sample}/final-viral-combined.fa"
    conda:
        "../envs/virsorter2.yaml"
    threads: threads
    log:
        "results/logs/virsorter2_filter/{sample}.log"
    shell:
        "virsorter run --working-dir $(dirname {output}) --seqfile {input[0]} --jobs {threads} --min-score 0.7 --include-groups RNA --verbose 2> {log}"

rule get_palmscan:
    output: 
        directory("resources/palmscan")
    shell: 
        """
        git -C resources clone https://github.com/rcedgar/palmscan.git
        chmod +x ./resources/palmscan/bin/palmscan
        """

rule run_palmscan:
    input: 
        "results/virsorter2_filter/{sample}/final-viral-combined.fa",
        "resources/palmscan"
    output:
        ppout="results/palmscan/{sample}_ppout.faa",
        ppout_nt="results/palmscan/{sample}_ppout_nt.fa",
        report="results/palmscan/{sample}_report.txt",
        fevout="results/palmscan/{sample}_fevout.fev"
    log:
        "results/logs/palmscan/{sample}.log"
    shell:
        """
        ./resources/palmscan/bin/palmscan -search_pp {input[0]} -rdrp -ppout {output.ppout} -ppout_nt {output.ppout_nt} -report {output.report} -fevout {output.fevout} 2> {log}
        """

rule get_palmdb:
    output:
        directory("resources/palmdb")
    conda:
        "../envs/blast.yaml"
    shell:
        """
        git -C resources clone https://github.com/rcedgar/palmdb.git
        gunzip resources/palmdb/2021-03-14/named.fa.gz
        makeblastdb -in resources/palmdb/2021-03-14/named.fa -dbtype prot
        cd-hit -i resources/palmdb/2021-03-14/named.fa -o resources/palmdb/2021-03-14/named.fa.cdhit.90 -c 0.9
        """

rule align_palmprint:
    input: 
        "results/palmscan/{sample}_ppout.faa",
        "resources/palmdb"
    output: 
        "results/blast/{sample}_named_blast.tsv"
    threads: threads
    conda:
        "../envs/blast.yaml"
    log:
        "results/logs/blast/{sample}.log"
    shell:
        """
        blastp -query {input[0]} -db {input[1]}/2021-03-14/named.fa -out {output} -outfmt 6 -num_alignments 10 -num_threads {threads} 2> {log}
        """

rule annotate_palmprint:
    input: 
        blast="results/blast/{sample}_named_blast.tsv"
    output: 
        "results/annotate_palmprint/{sample}_annotated.tsv"
    run: 
        from Bio import SeqIO
        import re

        tax_dict = dict()
        records = SeqIO.parse("resources/palmdb/2021-03-14/named.fa", "fasta")
        for record in records:
            tax_dict[record.id] = record.description

        blast_dict = dict()
        with open(input.blast) as blast_handle:
            for line in blast_handle:
                fields = line.split("\t")
                identity = float(fields[2])
                if not fields[0] in blast_dict:
                    blast_dict[fields[0]] = {90: list(), 75: list(), 40: list()}
                if identity >= 90:
                    blast_dict[fields[0]][90].append(fields[1])
                if identity >= 75:
                    blast_dict[fields[0]][75].append(fields[1])
                if identity >= 40:
                    blast_dict[fields[0]][40].append(fields[1])

        annotated_dict = dict()

        for key in blast_dict.keys():
            for thres in blast_dict[key]:
                blast_count = dict()
                if not blast_dict[key][thres]:
                    continue
                for blast_result in blast_dict[key][thres]:
                    tax_full = tax_dict[blast_result]
                    if thres == 90:
                        rank = re.findall(r'species:.+?$',tax_full)
                    elif thres == 75:
                        rank = re.findall(r'genus:.+?,',tax_full)
                    elif thres == 40:
                        rank = re.findall(r'family:.+?,',tax_full)
                    if not rank:
                        continue
                    rank = rank[0]
                    if rank in blast_count:
                        blast_count[rank]+=1
                    else:
                        blast_count[rank]=1

                if not blast_count:
                    continue
                max_key = max(blast_count, key=blast_count.get)
                b_out = blast_dict[key][thres]
                if blast_count[max_key] > len(b_out)/2:
                    for hit in blast_dict[key][thres]:
                        if max_key in tax_dict[hit]:
                            tax_full = tax_dict[hit]
                    annotated_dict[key] = f"{max_key}\t{tax_full}"
                    break
                else:
                    annotated_dict[key] = f".\t"

        with open(str(output), 'w') as out:
            for key in annotated_dict.keys():
                out.write(f"{key}\t{annotated_dict[key]}\n")

rule create_family_cluster:
    input: 
        "results/annotate_palmprint/{sample}_annotated.tsv",
        "results/blast/{sample}_named_blast.tsv",
        "results/palmscan/{sample}_ppout.faa",
        "resources/palmdb"
    output:
        "results/families/{sample}/{sample}.ok"
    run:
        from Bio import SeqIO
        import re
        
        with open(output[0], 'w') as out:
            out.write('ok')

        blast_dict = dict()
        with open(input[1]) as blast_handle:
            for line in blast_handle:
                fields = line.split("\t")
                if not fields[1] in blast_dict:
                    blast_dict[fields[1]]=1

        seqs_id_dict = dict()

        with open(input[0]) as ann_handle:
            for line in ann_handle:
                fields = line.split("\t")
                if fields[1] == ".":
                    continue
                family = re.findall(r'family:(.+?),',fields[2])
                if not family:
                    continue
                family = family[0]
                if not family in seqs_id_dict:
                    seqs_id_dict[family] = list()
                seqs_id_dict[family].append(fields[0])

        tax_dict = dict()
        records = SeqIO.parse("resources/palmdb/2021-03-14/named.fa", "fasta")
        for record in records:
            if record.id in blast_dict:
                family = re.findall(r'family:(.+?),',record.description)
                if not family:
                    continue
                family = family[0]
                if family in seqs_id_dict:
                    seqs_id_dict[family].append(record.id)

        seq_dict_ann = SeqIO.to_dict(SeqIO.parse(input[2], "fasta"))
        seq_dict_named = SeqIO.to_dict(SeqIO.parse("resources/palmdb/2021-03-14/named.fa.cdhit.90", "fasta"))

        seq_dict_ann.update(seq_dict_named)

        for family in seqs_id_dict:
            family_records = list()
            for seq_id in seqs_id_dict[family]:
                if seq_id in seq_dict_ann:
                    record = seq_dict_ann[seq_id]
                    specie = re.findall(r'species:(.+?)$',record.description)
                    if specie:
                        specie = specie[0]
                        record.id += f"_{specie}"
                        record.id = record.id.replace(" ", "_")
                        print(record.id)
                    family_records.append(record)
            
            out_dir = '/'.join(output[0].split('/')[:-1])
            with open(f"{out_dir}/{family}.fasta", 'w') as fout:
                SeqIO.write(family_records, fout, "fasta")

rule phylogenic_tree:
    input: 
        "results/families/{sample}/{sample}.ok"
    output: 
        "results/trees/{sample}/trees.ok"
    conda:
        "../envs/trees.yaml"
    shell:
        """
        for fasta in `ls $(dirname {input[0]})/*.fasta`
        do
            out_msa=`echo $(dirname {output[0]})/$(basename $fasta).msa`
            muscle -in $fasta -verbose  -quiet -out $out_msa
            fasttree_out=$(dirname {output[0]})/$(basename $fasta).tree
            FastTree $out_msa > $fasttree_out
            touch {output[0]}
        done
        """ 

rule plot_trees:
    input: 
        "results/trees/{sample}/trees.ok"
    output: 
        "results/trees/{sample}/trees_plot.ok"
    shell:
        """
        for tree in `ls $(dirname {input[0]})/*.tree`
        do
        png=`echo $(dirname {output[0]})/$(basename $tree).png`
        python workflow/scripts/plot_tree.py --tree_file $tree --out_file $png
        touch {output[0]}
        done
        """
