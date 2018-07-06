__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-07-06"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for processing SAM/BAM files
For usage, include this in your workflow.
"""

# import other packages
import os
import fnmatch
from snakemake.exceptions import MissingInputException

# set some local variables
home = os.environ['HOME']

rule:
    version: 0.3

# rules
rule bam_quality_filter:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bowtie2_pe.output
    output:
        temp("{runID}/{outdir}/{reference_version}/bowtie2/quality_filtered/{unit}.bam")
    shell:
        "samtools view -b -h -q {params.qual} {input} > {output}"

rule bam_sort:
    params:
        qual = config["alignment_quality"]
    threads:
        4
    input:
        rules.bam_quality_filter.output
    output:
        "{runID}/{outdir}/{reference_version}/bowtie2/{unit}.sorted.bam"
    shell:
        "samtools sort -@ {threads} {input} -T {wildcards.unit}.Q{params.qual}.sorted -o {output}"

rule bam_mark_duplicates:
    params:
        qual = config["alignment_quality"],
        picard = home + config["picard"],
        temp = home + config["temp_dir"]
    input:
        rules.bam_sort.output
    output:
        temp("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.duplicates_marked.bam")
    shell:
        """
            java -Djava.io.tmpdir={params.temp} \
            -Xmx24G \
            -jar {params.picard} MarkDuplicates \
            INPUT={input}\
            OUTPUT={output}\
            ASSUME_SORTED=TRUE\
            METRICS_FILE={output}.metrics.txt
        """

rule bam_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_mark_duplicates.output
    output:
        temp("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.duplicates_marked.bam.bai")
    shell:
        "samtools index {input} {output}"

rule bam_rmdup:
    input:
        rules.bam_mark_duplicates.output
    output:
        protected("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.bam")
    shell:
        "samtools rmdup {input} {output}"

rule bam_rmdup_index:
    params:
        qual = config["alignment_quality"]
    input:
        rules.bam_rmdup.output
    output:
        protected("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.bam.bai")
    shell:
        "samtools index {input} {output}"
