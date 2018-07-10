__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-09-07"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for running MACS2 for ChIP peak calling.
For usage, include this in your workflow.
"""

def get_replicates_input(wildcards):
    fn = []
    for i in config["samples"][][wildcards.Input][wildcards.sample + "_" + wildcards.Input + "_" + wildcards.digest]:
        fn.append("processed_data/hg38/duplicates_removed/" + i + ".DeDup.sorted.fastq_q20.bam")
    return(fn)

def get_replicates_chip(wildcards):
    fn = []
    for i in config["samples"][wildcards.digest][wildcards.ChIP][wildcards.sample + "_" + wildcards.ChIP + "_" + wildcards.digest]:
        fn.append("processed_data/hg38/duplicates_removed/" + i + ".DeDup.sorted.fastq_q20.bam")
    return(fn)

rule macs2_predictd:
    params:
        macs2_dir = config["macs2_dir"]
    input:
        chip = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{ChIP}.bam"
    output:
        "{runID}/{outdir}/{reference_version}/macs2/predictd/{ChIP}"
    log:
        "{runID}/{outdir}/{reference_version}/macs2/predictd/{ChIP}/log.txt"
    shell:
        """
            {params.macs2_dir}/macs2 predictd --ifile {input.chip}\
                                              --gsize hs\
                                              --mfold 3 50\
                                              --outdir {output} 1 >> {log} 2 >> {log}
        """


rule macs2_callpeak:
    params:
        extsize = config["parameters"]["macs2"]["extsize"],
        macs2_dir = config["macs2_dir"]
    input:
        input = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{Input}.bam",
        chip = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{ChIP}.bam"
    output:
        "{runID}/{outdir}/{reference_version}/macs2/callpeak/{digest}/{ChIP}_vs_{Input}/{sample}"
    shell:
        """
            {params.macs2_dir}/macs2 callpeak -B \
                                              -t {input.chip}\
                                              -c {input.input}\
                                              -n {wildcards.sample}\
                                              --nomodel\
                                              --extsize {params.extsize}\
                                              --outdir {output}
        """
