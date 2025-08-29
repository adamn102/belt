# Snakefile
import os
import glob

configfile: "configs/XTNNF9.yaml"

#########################################################
# Rule: all
# This rule collects all final consensus FASTA outputs.
#########################################################

import os
import glob

def list_all_compiled_consensus():
    return [f"data/8-compiled_consensus/{sample}-compiled-consensus.fa" for sample in config["samples"]]

def list_all_consensus_stats():
    return [f"data/10-transgene_mutations/per_read/{sample}_consensus_summary.csv" for sample in config["samples"]]

rule all:
    input:
        expand("data/2-fastq_filtered/{sample}.fastq",sample=config["samples"]),
        expand("data/4-sorted_reads/{sample}.bam.bai", sample=config["samples"]),
        expand("data/5-umi_maps/{sample}_Paired_UMIs_read_names.csv",sample=config["samples"]),
        expand("data/6-umi_groups/{sample}.grouped.bam",sample=config["samples"]),
        expand("data/7-consensus_calls/{sample}/consensus",sample=config["samples"]),
        #expand("data/9-mapped_consensus/{sample}.bam",sample=config["samples"]),
        #list_all_consensus_stats(),


#########################################################
# Rule: filter_sequences
#########################################################

rule filter_sequences:
    input:
        R = "data/1-fastq_files/{sample}.fastq"
    params:
        l = "1000",
        q = "25",
        u = "50"
    log:
        h = "logs/fastq/{sample}.html",
        j = "logs/fastq/{sample}.json"
    output:
        R = "data/2-fastq_filtered/{sample}.fastq"
    shell:
        """
        fastplong --in {input.R} --out {output.R} \
        --length_required {params.l} \
        --qualified_quality_phred {params.q} \
        --unqualified_percent_limit {params.u} \
        --disable_adapter_trimming \
        --html {log.h} --json {log.j}
        """

#########################################################
# Reference and Bed helper functions
#########################################################

def get_reference(wildcards):
    reference = config["samples"][wildcards.sample].split("/")[2].split("_R1")[0].split(".")[0].split("-")[0]
    return "data/0-reference_sequence/" + reference + ".fa",

def get_bed(wildcards):
    reference = config["samples"][wildcards.sample].split("-")[-1].split('_')[0]
    return "data/0-bed_file/" + reference + ".bed"

#########################################################
# Rule: map_reads
#########################################################

rule map_reads:
    input:
        ref = get_reference,
        R = "data/2-fastq_filtered/{sample}.fastq"
    log:
        L = "logs/minimap2/{sample}.log"
    params:
        threads = 8
    output:
        O = "data/3-mapped_reads/{sample}.bam"
    shell:
        """
        minimap2 -t {params.threads} -a {input.ref} {input.R} > {output.O} 2> {log.L}
        """

#########################################################
# Rule: samtools_sort
#########################################################

rule samtools_sort:
    input:
        "data/3-mapped_reads/{sample}.bam"
    output:
        "data/4-sorted_reads/{sample}.bam"
    shell:
        "samtools sort {input} -o {output}"

#########################################################
# Rule: samtools_index
#########################################################

rule samtools_index:
    input:
        "data/4-sorted_reads/{sample}.bam"
    output:
        "data/4-sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

#########################################################
# Rule: umi_map
#########################################################

def get_barcode(wildcards):
    return config["barcodes"][wildcards.sample]

def get_whitelist(wildcards):
    file = config["whitelist"][wildcards.sample]
    if file == "cluster":
        return "cluster"
    elif file == "false":
        return "false"
    else:
        return f"data/0-library_whitelist/{file}.csv"

rule umi_map:
    input:
        reference = get_reference,
        alignment = "data/4-sorted_reads/{sample}.bam",
        bai       = "data/4-sorted_reads/{sample}.bam.bai"
    params:
        minReadLength = 1000,
        whitelist = get_whitelist,
        motif = get_barcode
    output:
        outfile = "data/5-umi_maps/{sample}_Paired_UMIs_read_names.csv",
        analysis = directory("data/11-analyses/umi_maps/{sample}")
    shell:
        """
        belt umi-map \
        --reference {input.reference} \
        --alignment {input.alignment} \
        --min_read_length {params.minReadLength} \
        --motif {params.motif} \
        --whitelist {params.whitelist} \
        --output_dir data/5-umi_maps \
        --analysis_out {output.analysis}
        """

#########################################################
# Rule: umi_group
#########################################################

rule umi_group:
    input:
        reference = get_reference,
        alignment = "data/4-sorted_reads/{sample}.bam",
        correctedUMIs = "data/5-umi_maps/{sample}_Paired_UMIs_read_names.csv"
    params:
        grouping = "5pUMI",
        threshold = 10
    output:
        outfile = "data/6-umi_groups/{sample}.grouped.bam"
    shell:
        """
        belt umi \
        --reference {input.reference} \
        --alignment {input.alignment} \
        --corrected_umis {input.correctedUMIs} \
        --grouping {params.grouping} \
        --threshold {params.threshold} \
        --output_dir data/6-umi_groups
        """


#########################################################
# Rule: consensus_call (per-UMI)
#########################################################

def get_consensus_reference(wildcards):
    reference = config["samples"][wildcards.sample].split("/")[2].split("_R1")[0].split(".")[0].split("-")[0]
    return "data/0-consensus_reference/" + reference + ".fa"
    
rule consensus_from_grouped_bam:
    input:
        grouped_bam = "data/6-umi_groups/{sample}.grouped.bam",
        ref = get_consensus_reference
    output:
        manifest = "data/7-consensus_calls/{sample}/consensus/manifest.txt",
        consensus_dir = directory("data/7-consensus_calls/{sample}/consensus")
    threads: 16
    log:
        "logs/consensus/{sample}.log"
    shell:
        r"""
        python3 scripts/consensus_from_grouped_bam.py \
            --grouped-bam {input.grouped_bam} \
            --reference {input.ref} \
            --outdir data/7-consensus_calls \
            --sample {wildcards.sample} \
            --procs {threads} \
            --umi-threads 1 \
            --manifest {output.manifest} \
            --log-file {output.consensus_dir}/runner.log \
            >> {log} 2>&1
        """

#########################################################
# Rule: consensus_call compile by sample
#########################################################


def list_sample_consensus_fastas(wildcards):
    """Return list of consensus fasta files for this sample."""
    return sorted(
        glob.glob(f"data/7-consensus_calls/{wildcards.sample}/consensus/*_consensus.fasta")
    )

rule compile_consensus:
    input:
        files = list_sample_consensus_fastas,
        manifest = "data/7-consensus_calls/{sample}/consensus/manifest.txt"
    output:
        "data/8-compiled_consensus/{sample}-compiled-consensus.fa"
    shell:
        """
        cat {input.files} > {output}
        """
    
#########################################################
# Rule: map consensus sequences to reference
#########################################################
def get_reference_map_consensus(wildcards):
    reference = config["samples"][wildcards.sample].split("/")[2].split("_R1")[0].split(".")[0].split("-")[0]
    return "data/0-consensus_map_reference/" + reference + ".fa",
    

rule map_compiled_consensus:
    input:
        consensus = "data/8-compiled_consensus/{sample}-compiled-consensus.fa",
        ref = get_reference_map_consensus,
    output:
        "data/9-mapped_consensus/{sample}.bam"
    shell:
        """
        minimap2 -a --MD --eqx {input.ref} {input.consensus} | samtools sort > {output}
        samtools index {output}
        """

#########################################################
# Rule: Tabulate mutations by consensus read
#########################################################
def get_bed_annotations(wildcards):
    reference = config["samples"][wildcards.sample].split("/")[2].split("_R1")[0].split(".")[0].split("-")[0]
    return "data/0-bed_annotations/" + reference + ".bed",
    
    
rule compute_consensus_stats:
    input:
        bam = "data/9-mapped_consensus/{sample}.bam",
        ref = get_reference_map_consensus,
    params:
        bed = get_bed_annotations
    output:
        per_read = "data/10-transgene_mutations/per_read/{sample}_consensus_summary.csv",
        variants = "data/10-transgene_mutations/variants/{sample}_variant_table.csv",
        site_summary = "data/10-transgene_mutations/site_summary/{sample}_site_summary.csv"
    shell:
        """
        python3 scripts/full_compute_mutation_stats_from_bam.py \
            --input_bam {input.bam} \
            --reference_fasta {input.ref} \
            --per_read_csv {output.per_read} \
            --variants_csv {output.variants} \
            --site_summary_csv {output.site_summary}
        """
