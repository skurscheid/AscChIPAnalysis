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
        control = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{Input}.bam",
        chip = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{ChIP}.bam"
    output:
        "{runID}/{outdir}/{reference_version}/edd/{ChIP}_vs_{Input}/"
    log:
        "run_edd_{ChIP}_vs_{Input}.log"
    threads:
        8
    params:
        chromSizes = home + config["references"][REF_GENOME]["chromSizesUCSC"],
        gaps = home + config["references"][REF_GENOME]["gapsUCSC"],
        numTrials = "10000",
        fdr = "0.1",
        eddBinDir = home + config["edd_dir"]
    shell:
        """"
            {params.eddBinDir}/edd {params.chromSizes} {params.gaps} {input.chip} {input.control} {output}\
                                   -n {params.numTrials}\
                                   --fdr {params.fdr}\
                                   --nprocs {threads}\
                                   --write-log-ratios\
                                   --write-bin-scores 1>>{log} 2>>{log}\
        """"
