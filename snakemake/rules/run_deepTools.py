__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-01-17"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools analysis on ChIP-Seq data
For usage, include this in your workflow.
"""

def cli_parameters_computeMatrix(wildcards):
    a = config["program_parameters"][wildcards["application"]][wildcards["tool"]][wildcards["command"]]
    if wildcards["command"] == "reference-point":
        a["--referencePoint"] = wildcards.referencePoint
    return(a)

def cli_parameters_normalization(wildcards):
    if wildcards["norm"] == "RPKM":
        a = "--normalizeUsingRPKM"
    elif wildcards["norm"] == "1xcoverage":
        a = " ".join(("--normalizeTo1x", config["references"][REF_GENOME]["effectiveSize"]))
    return(a)

def cli_parameters_bamCoverage(wildcards):
    a = config["program_parameters"][wildcards["application"]][wildcards["tool"]][wildcards["mode"]]
    b = str()
    for (key, val) in a.items():
        if val == " ":
            f = key + " "
            b = b + f
        else:
            f = key + "=" + val + " "
            b = b + f
    if wildcards["mode"] == "MNase":
        b = b + "--MNase"
    return(b.rstrip())

def get_computeMatrix_input(wildcards):
    fn = []
    path = "/".join((wildcards["assayID"],
                     wildcards["runID"],
                     config["processed_dir"],
                     config["references"][REF_GENOME]["version"][0],
                     wildcards["application"],
                     "bamCoverage",
                     wildcards["mode"],
                     wildcards["duplicates"]))
    for i in config["samples"][wildcards["assayID"]][wildcards["runID"]]:
        fn.append("/".join((path, "_".join((i, wildcards["mode"], "RPKM.bw")))))
    return(fn)


rule bamCoverage:
    version:
        0.1
    params:
        deepTools_dir = home + config["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cli_parameters_bamCoverage
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        bam = lambda wildcards: wildcards.assayID + "/" + wildcards.runID + "/" + wildcards.outdir + "/" + wildcards.reference_version + "/bowtie2/" + wildcards.duplicates + "/" +  wildcards.sample + ".Q" + config["alignment_quality"] + ".sorted.bam"
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{sample}_{mode}_{norm}.bw"
    shell:
        """
        {params.deepTools_dir}/bamCoverage --bam {input.bam} \
                                           --outFileName {output} \
                                           --outFileFormat bigwig \
                                           {params.program_parameters} \
                                           --numberOfProcessors {threads} \
                                           --normalizeUsingRPKM \
                                           --ignoreForNormalization {params.ignore}
        """

rule computeMatrix:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        file = get_computeMatrix_input,
        region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards.reference_version][wildcards.region]
    output:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{region}_{mode}.matrix.gz"
    shell:
        """
            {params.deepTools_dir}/computeMatrix {wildcards.command} \
                                                 --regionsFileName {input.region} \
                                                 --scoreFileName {input.file} \
                                                 --missingDataAsZero \
                                                 --skipZeros \
                                                 --numberOfProcessors {threads} \
                                                 {params.program_parameters} \
                                                 --outFileName {output.matrix_gz}
        """


rule plotProfile:
    version:
        0.1
    params:
        deepTools_dir = home + config["deepTools_dir"],
    input:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{region}_{mode}.matrix.gz"
    output:
        figure = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{region}.pdf",
        data = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{region}.data",
        regions = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/{plotType}.{mode}.{region}.bed"
    shell:
        """
            {params.deepTools_dir}/plotProfile --matrixFile {input.matrix_gz} \
                                               --outFileName {output.figure} \
                                               --outFileNameData {output.data} \
                                               --outFileSortedRegions {output.regions} \
                                               --plotType {wildcards.plotType}
        """

rule bam_compare_pooled_replicates:
    version:
        0.1
    params:
        deepTools_dir = home + config["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cli_parameters_bamCoverage
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        control = "{assayID}/{runID}/{outdir}/{reference_version}/samtools/merge/{duplicates}/{control}.bam",
        treatment = "{assayID}/{runID}/{outdir}/{reference_version}/samtools/merge/{duplicates}/{treatment}.bam"
    output:
        "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{scaleFactors}/{treatment}_vs_{control}_{mode}_{ratio}_RPKM.bw"
    shell:
        """
            {params.deepTools_dir}/bamCompare --bamfile1 {input.treatment} \
                                              --bamfile2 {input.control} \
                                              --outFileName {output} \
                                              --scaleFactorsMethod {wildcards.scaleFactors} \
                                              --ratio {wildcards.ratio} \
                                              --numberOfProcessors {threads} \
                                              --normalizeUsingRPKM \
                                              --ignoreForNormalization {params.ignore}
        """

rule bam_coverage_pooled_replicates:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        ignore = config["program_parameters"]["deepTools"]["ignoreForNormalization"],
        program_parameters = cli_parameters_bamCoverage,
        normalization = lambda wildcards: cli_parameters_normalization(wildcards)
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        bam = "{assayID}/{runID}/{outdir}/{reference_version}/samtools/merge/{duplicates}/{sample_group}.bam",
        bai = "{assayID}/{runID}/{outdir}/{reference_version}/samtools/merge/{duplicates}/{sample_group}.bam.bai"
    output:
        bigwig = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{mode}/{duplicates}/{sample_group}_{mode}_{norm}.bw"
    shell:
        """
            {params.deepTools_dir}/bamCoverage --bam {input.bam} \
                                               --outFileName {output.bigwig} \
                                               --outFileFormat bigwig \
                                               {params.program_parameters} \
                                               --numberOfProcessors {threads} \
                                               {params.normalization} \
                                               --ignoreForNormalization {params.ignore}
        """
