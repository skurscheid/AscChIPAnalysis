__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2016-09-07"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for running MACS2 for ChIP peak calling.
For usage, include this in your workflow.
"""

rule macs2_predictd:
    params:
        macs2_dir = home + config["macs2_dir"]
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
                                              --outdir {output}
        """


rule macs2_callpeak_cutoff_analysis:
    params:
        extsize = 180, # estimated from predictd output
        macs2_dir = home + config["macs2_dir"]
    input:
        input = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{Input}.bam",
        chip = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{ChIP}.bam"
    output:
        "{runID}/{outdir}/{reference_version}/macs2/callpeak/{ChIP}_vs_{Input}/"
    shell:
        """
            {params.macs2_dir}/macs2 callpeak -B \
                                              -t {input.chip}\
                                              -c {input.input}\
                                              -n {wildcards.ChIP}\
                                              --nomodel\
                                              --extsize {params.extsize}\
                                              --cutoff-analysis\
                                              --outdir {output}
        """
rule macs2_call_broadpeaks:
    params:
        extsize = 180, # estimated from predictd output
        macs2_dir = home + config["macs2_dir"]
    input:
        input = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{Input}.bam",
        chip = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{ChIP}.bam"
    output:
        "{runID}/{outdir}/{reference_version}/macs2/callpeak_broad/{ChIP}_vs_{Input}/"
    shell:
        """
            {params.macs2_dir}/macs2 callpeak -B \
                                              -t {input.chip}\
                                              -c {input.input}\
                                              -n {wildcards.ChIP}\
                                              --nomodel\
                                              --extsize {params.extsize}\
                                              --broad\
                                              --mfold 3 50\
                                              --outdir {output}
        """
