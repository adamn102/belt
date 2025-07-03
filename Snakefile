configfile: "configs/20250703.yaml"

rule all:
    input:
        #expand("data/5-calls/{sample}.vcf", sample=config["samples"]),
        expand("data/1-fastq_files/{sample}.fastq",sample=config["samples"]),
        expand("data/3-mapped_reads/{sample}.sam",sample=config["samples"]),
        expand("data/4-sorted_reads/{sample}.bam.bai", sample=config["samples"]),
        expand("data/5-barcode_counts/{sample}",sample=config["samples"])




#filter sequences
rule filter_sequences:
    input:
        R1 = "data/1-fastq_files/{sample}.fastq",
    params:
        tm="0",
        qs="15", # the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.
    log:
        "logs/fastq/{sample}.html"
    output:
        R1 = "data/2-fastq_filtered/{sample}.fastq",
    shell:
        "fastplong --in {input.R1}  --out {output.R1} -h "
        " -q {params.qs} --html {log}"


def get_reference(wildcards):
    reference = config["samples"][wildcards.sample].split("/")[2].split(".")[0]
    return "data/0-reference_index/" + reference + ".fa",

rule minimap2_map:
    input:
        ref = get_reference,
        R = "data/2-fastq_filtered/{sample}.fastq",
    output:
        directory("data/3-mapped_reads/{sample}.sam"),
    log:
        "logs/minimap2/{sample}.log"
    shell:
        """
        minimap2 -ax map-ont {input.ref} {input.R} > {output}
        """

rule samtools_samToBam:
    input:
        "data/3-mapped_reads/{sample}.sam"
    output:
        "data/3-mapped_reads/{sample}.bam"
    shell:
        " samtools view -b -S  {input} > {output}"

rule samtools_sort:
    input:
        "data/3-mapped_reads/{sample}.bam"
    output:
        "data/4-sorted_reads/{sample}.bam"
    shell:
        "samtools sort {input} -o {output}"

rule samtools_index:
    input:
        "data/4-sorted_reads/{sample}.bam"
    output:
        "data/4-sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

def get_barcode(wildcards):
    return config["barcodes"][wildcards.sample],

def get_min_length(wildcards):
    return config["min_lengths"][wildcards.sample],

rule barcode_report:
    input:
        reference = get_reference,
        bam = "data/4-sorted_reads/{sample}.bam",
        bai = "data/4-sorted_reads/{sample}.bam.bai"
    output:
        directory("data/5-barcode_counts/{sample}")
    params:
        barcode = get_barcode,
        min_length = get_min_length,
    shell:
        "python3 scripts/belt.py extract --reference {input.reference} --min_read_length {params.min_length} --alignment {input.bam} --motif {params.barcode} --output {output}"
