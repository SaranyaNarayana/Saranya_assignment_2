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


# Step 7: Align reads with BWA (with read groups)
rule align_reads:
    input:
        ref=f"{RAW}/reference.fasta",
        fastq=f"{RAW}/{SRA}.fastq",
        bwt=f"{RAW}/reference.fasta.bwt"
    output:
        f"{ALIGNED}/aligned.sam"
    shell:
        r"""
        bwa mem -R '@RG\tID:1\tLB:lib1\tPL:illumina\tPU:unit1\tSM:sample1' {input.ref} {input.fastq} > {output}
        """


# Step 8: Convert SAM to sorted BAM
rule sort_bam:
    input:
        f"{ALIGNED}/aligned.sam"
    output:
        f"{ALIGNED}/aligned.sorted.bam"
    shell:
        "samtools view -b {input} | samtools sort -o {output}"


# Step 9: Validate BAM file
rule validate_bam:
    input:
        f"{ALIGNED}/aligned.sorted.bam"
    output:
        touch(f"{ALIGNED}/.validated")
    shell:
        "gatk ValidateSamFile -I {input} -MODE SUMMARY"


# Step 10: Mark duplicates 
rule mark_duplicates:
    input:
        bam=f"{ALIGNED}/aligned.sorted.bam",
        check=f"{ALIGNED}/.validated"
    output:
        bam=f"{ALIGNED}/dedup.bam",
        metrics=f"{ALIGNED}/dup_metrics.txt"
    shell:
        "gatk MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics}"
