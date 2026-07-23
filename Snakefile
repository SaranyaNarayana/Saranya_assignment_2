# Snakefile

SRA = "SRR1972739"
REF_ID = "AF086833.2"

RESULTS = "results"
RAW = f"{RESULTS}/raw"
ALIGNED = f"{RESULTS}/aligned"
VARIANTS = f"{RESULTS}/variants"
ANNOTATED = f"{RESULTS}/annotated"
QC = f"{RESULTS}/qc"
SNPEFF = f"{RESULTS}/snpEff"
SNPEFF_DATA = f"{SNPEFF}/data/reference_db"


rule all:
    input:
        f"{ANNOTATED}/annotated_variants.vcf",
        f"{QC}/{SRA}_fastqc.html"

# Step 1: Download reference genome
rule download_reference:
    output:
        f"{RAW}/reference.fasta"
    shell:
        "efetch -db nucleotide -id {REF_ID} -format fasta > {output}"
