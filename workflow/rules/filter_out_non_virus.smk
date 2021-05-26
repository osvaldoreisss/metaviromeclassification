rule filter_out_non_virus:
    input: get_fasta
    output: "results/filter_out_non_virus/{sample}.fa"
    shell:
        """
        echo {input}
        cat {input} > {output}
        """ 


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
        "results/filter_out_non_virus/{sample}.fa",
        "resources/palmscan"
    output:
        ppout="results/palmscan/{sample}_ppout.faa",
        ppout_nt="results/palmscan/{sample}_ppout_nt.fa",
        report="results/palmscan/{sample}_report.txt",
        fevout="results/palmscan/{sample}_fevout.fev"
    shell:
        """
        ./resources/palmscan/bin/palmscan -search_pp {input[0]} -rdrp -ppout {output.ppout} -ppout_nt {output.ppout_nt} -report {output.report} -fevout {output.fevout}
        """

rule get_palmdb:
    output:
        directory("resources/palmdb")
    conda:
        "../envs/blast.yaml"
    shell:
        """
        git -C resources clone git@github.com:rcedgar/palmdb.git
        gunzip resources/palmdb/2021-03-14/named.fa.gz
        makeblastdb -in resources/palmdb/2021-03-14/named.fa -dbtype prot
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
    shell:
        """
        blastp -query {input[0]} -db {input[1]}/2021-03-14/named.fa -out {output} -outfmt 6 -num_alignments 10 -num_threads {threads}
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
                    print(tax_full)
                    if thres == 90:
                        rank = re.findall(r'species:.+?$',tax_full)
                    elif thres == 75:
                        rank = re.findall(r'genus:.+?,',tax_full)
                    elif thres == 40:
                        rank = re.findall(r'family:.+?,',tax_full)
                    print(rank)
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
                if blast_count[max_key] >= len(blast_result)/2:
                    annotated_dict[key] = max_key
                    break
                else:
                    annotated_dict[key] = "."
        print(output)
        with open(str(output), 'w') as out:
            for key in annotated_dict.keys():
                out.write(f"{key}\t{annotated_dict[key]}\n")