from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"

samples = pd.read_csv(config["samples"], sep=",").set_index("sample", drop=False)

def get_fasta(wildcards):
    fasta = samples.loc[(wildcards.sample), ["assembly"]].dropna()
    return f"{fasta.assembly}"