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
        ppout="results/palmscan/{sample}/{sample}_ppout.faa",
        ppout_nt="results/palmscan/{sample}/{sample}_ppout_nt.fa",
        report="results/palmscan/{sample}/{sample}_report.txt",
        fevout="results/palmscan/{sample}/{sample}_fevout.fev"
    shell:
        """
        ./resources/palmscan/bin/palmscan -search_pp {input[0]} -rdrp -ppout {output.ppout} -ppout_nt {output.ppout_nt} -report {output.report} -fevout {output.fevout}
        """