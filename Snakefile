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


# Step 2: Download sequencing data 
rule download_sra:
    output:
        f"{RAW}/{SRA}.fastq"
    params:
        sra=SRA,
        raw=RAW
    shell:
        """
        prefetch {params.sra} -O {params.raw}
        fastq-dump -X 10000 {params.raw}/{params.sra}/{params.sra}.sra -O {params.raw}
        """


# Step 3: FastQC on raw reads 
rule fastqc_raw:
    input:
        f"{RAW}/{SRA}.fastq"
    output:
        f"{QC}/{SRA}_fastqc.html"
    params:
        qc=QC
    shell:
        "fastqc -o {params.qc} {input}"


# Step 4: Index reference genome 
rule index_reference_samtools:
    input:
        f"{RAW}/reference.fasta"
    output:
        f"{RAW}/reference.fasta.fai"
    shell:
        "samtools faidx {input}"

# Step 5: Build BWA index 
rule bwa_index:
    input:
        f"{RAW}/reference.fasta"
    output:
        f"{RAW}/reference.fasta.bwt"
    shell:
        "bwa index {input}"


# Step 6: Create FASTA sequence dictionary (GATK)
rule create_dict:
    input:
        f"{RAW}/reference.fasta"
    output:
        f"{RAW}/reference.dict"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"
