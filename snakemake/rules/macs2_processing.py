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
            if [ ! -d {output} ]; then mkdir -p {output}; fi;\
            {params.macs2_dir}/macs2 predictd --ifile {input.chip}\
                                              --gsize hs\
                                              --mfold 3 50\
                                              --outdir {output}\
                                              $1>>{log}\
                                              $2>>{log}
        """


# rule macs2_callpeak:
#     params:
#         extsize = config["parameters"]["macs2"]["extsize"],
#         macs2_dir = config["macs2_dir"]
#     input:
#         input = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{Input}.bam",
#         chip = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{ChIP}.bam"
#     output:
#         "{runID}/{outdir}/{reference_version}/macs2/callpeak/{digest}/{ChIP}_vs_{Input}/{sample}"
#     shell:
#         """
#             {params.macs2_dir}/macs2 callpeak -B \
#                                               -t {input.chip}\
#                                               -c {input.input}\
#                                               -n {wildcards.sample}\
#                                               --nomodel\
#                                               --extsize {params.extsize}\
#                                               --outdir {output}
#         """
