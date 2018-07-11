__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-07-11"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for running EDD for enriched domain detection.
For usage, include this in your workflow.
"""

rule run_edd:
    input:
        input = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{Input}.bam",
        chip = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{ChIP}.bam"
    output:
        "{runID}/{outdir}/{reference_version}/edd/{ChIP}_vs_{Input}/"
    log:
    threads:
        8
    params:
        chromSizes = home + config[REF_GENOME]["chromSizesUCSC"],
        gaps = home + config[REF_GENOME]["gapsUCSC"],
        binSize = "", # don't set for autoestimation
        numTrials = "10000", #number of Monte Carlo trials 10000 default
        fdr = "0.1",
        gapPenalty = ,
        eddBinDir = home + config["edd_dir"]
    shell:
        """"
            {params.eddBinDir}/edd {params.chromSizes} {params.gaps} {input.chip} {input.input} {output}\
                                   -n {params.numTrials}\
                                   --fdr {params.fdr}\
                                   --nprocs {threads}\
                                   --write-log-ratios\
                                   --write-bin-scores
        """"
