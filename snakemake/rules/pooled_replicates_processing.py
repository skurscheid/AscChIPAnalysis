from snakemake.exceptions import MissingInputException
import os

home = os.environ['HOME']

REF_GENOME = config["references"]["genomes"][1]

# run parameters as variables
ASSAYID = "ChIP-Seq"
RUNID = "merged"
OUTDIR = config["processed_dir"]
REFVERSION = reference_version = config["references"][REF_GENOME]["version"][0]
QUALITY = config["alignment_quality"]

def cli_parameters_computeMatrix(wildcards):
    a = config["program_parameters"][wildcards["application"]]["computeMatrix"][wildcards["command"]]
    if wildcards["command"] == "reference-point":
        a["--referencePoint"] = wildcards.referencePoint
    return(a)


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

rule run_computeMatrix_pooled_replicates:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{sampleGroup}_{region}_{mode}.{norm}.matrix.gz",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               command = ["reference-point", "scale-regions"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               referencePoint = "TSS",
               sampleGroup = ["H2AZ_10A_high", "H2AZ_TGFb_10A", "Inp_10A_WT_high", "Inp_10A_TGFb_high", "Inp_shZ_10A_high"],
               region = ["allGenes", "TanEMTup", "TanEMTdown"],
               mode = ["MNase", "normal"],
               norm = ["RPKM", "1xcoverage"])

rule run_computeMatrix_pooled_replicates_single_matrix:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{sampleGroup}_{region}_{mode}.{norm}.matrix.gz",
               assayID = ASSAYID,
               runID = RUNID,
               outdir = OUTDIR,
               reference_version = REFVERSION,
               application = "deepTools",
               command = ["reference-point", "scale-regions"],
               duplicates = ["duplicates_marked", "duplicates_removed"],
               referencePoint = "TSS",
               sampleGroup = "allSamples",
               region = ["allGenes", "TanEMTup", "TanEMTdown"],
               mode = ["MNase", "normal"],
               norm = ["RPKM", "1xcoverage"])

rule run_plotProfile_pooled_replicates:
    input:
        expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/allSamples_{plotType}.{mode}.{region}.{norm}.{suffix}",
                assayID = ASSAYID,
                runID = RUNID,
                outdir = OUTDIR,
                reference_version = REFVERSION,
                application = "deepTools",
                tool = "plotProfile",
                command = ["reference-point", "scale-regions"],
                duplicates = ["duplicates_marked", "duplicates_removed"],
                referencePoint = "TSS",
                plotType = "se",
                region = ["allGenes", "TanEMTup", "TanEMTdown"],
                mode = ["MNase", "normal"],
                suffix = ["pdf", "bed", "data"],
                norm = ["RPKM", "1xcoverage"])

rule computeMatrix_pooled_replicates:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        file = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/bamCoverage/{mode}/{duplicates}/{sampleGroup}_{mode}_{norm}.bw",
        region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards.region]
    output:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/{sampleGroup}_{region}_{mode}.{norm}.matrix.gz"
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

rule computeMatrix_pooled_replicates_single_matrix:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        program_parameters = lambda wildcards: ' '.join("{!s}={!s}".format(key, val.strip("\\'")) for (key, val) in cli_parameters_computeMatrix(wildcards).items())
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        file = lambda wildcards: expand("{assayID}/{runID}/{outdir}/{reference_version}/{application}/bamCoverage/{mode}/{duplicates}/{sampleGroup}_{mode}_{norm}.bw",
                                       assayID = ASSAYID,
                                       runID = RUNID,
                                       outdir = OUTDIR,
                                       reference_version = REFVERSION,
                                       application = "deepTools",
                                       duplicates = wildcards["duplicates"],
                                       sampleGroup = ["H2AZ_10A_high", "H2AZ_TGFb_10A", "Inp_10A_WT_high", "Inp_10A_TGFb_high", "Inp_shZ_10A_high"],
                                       region = wildcards["region"],
                                       mode = wildcards["mode"],
                                       norm = wildcards["norm"]),
        region = lambda wildcards: home + config["program_parameters"]["deepTools"]["regionFiles"][wildcards["reference_version"]][wildcards.region]
    output:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/allSamples_{region}_{mode}.{norm}.matrix.gz"
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

rule plotProfile_pooled_replicates:
    version:
        0.1
    params:
        deepTools_dir = home + config["deepTools_dir"],
    input:
        matrix_gz = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/computeMatrix/{command}/{duplicates}/{referencePoint}/allSamples_{region}_{mode}.{norm}.matrix.gz"
    output:
        figure = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/allSamples_{plotType}.{mode}.{region}.{norm}.pdf",
        data = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/allSamples_{plotType}.{mode}.{region}.{norm}.data",
        regions = "{assayID}/{runID}/{outdir}/{reference_version}/{application}/{tool}/{command}/{duplicates}/{referencePoint}/allSamples_{plotType}.{mode}.{region}.{norm}.bed"
    shell:
        """
            {params.deepTools_dir}/plotProfile --matrixFile {input.matrix_gz} \
                                               --outFileName {output.figure} \
                                               --outFileNameData {output.data} \
                                               --outFileSortedRegions {output.regions} \
                                               --plotType {wildcards.plotType}
        """
