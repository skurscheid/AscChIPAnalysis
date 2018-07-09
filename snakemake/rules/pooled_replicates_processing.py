from snakemake.exceptions import MissingInputException
import os

def cli_parameters_computeMatrix(wildcards):
    a = config["program_parameters"][wildcards["application"]]["computeMatrix"][wildcards["command"]]
    if wildcards["command"] == "reference-point":
        a["--referencePoint"] = wildcards.referencePoint
    return(a)

rule bam_merge:
    version:
        0.4
    params:
        cwd = os.getcwd()
    threads:
        8
    input:
        bam = lambda wildcards: expand("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.{suffix}",
                                       runID = wildcards["runID"],
                                       outdir = wildcards["outdir"],
                                       reference_version = wildcards["reference_version"],
                                       unit = config["samples"]["ChIP-Seq"][wildcards["type"]][wildcards["sample"]],
                                       suffix = ["bam"]),
        index = lambda wildcards: expand("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.{suffix}",
                                       runID = wildcards["runID"],
                                       outdir = wildcards["outdir"],
                                       reference_version = wildcards["reference_version"],
                                       unit = config["samples"]["ChIP-Seq"][wildcards["type"]][wildcards["sample"]],
                                       suffix = ["bam.bai"]),
    output:
        protected("{runID}/{outdir}/{reference_version}/bowtie2/merged/{sample}.bam")
    run:
        if (len(input.bam) > 1):
            shell("samtools merge --threads {threads} {output} {input.bam}")
        else:
            shell("ln -s {params.cwd}/{input.bam} {output}")

rule index_merged:
    version:
        0.1
    input:
        rules.bam_merge.output
    output:
        protected("{runID}/{outdir}/{reference_version}/bowtie2/merged/{sample}.bam.bai")
    shell:
        "samtools index {input} {output}"

rule bam_compare_pooled_replicates:
    version:
        0.1
    params:
        deepTools_dir = home + config["deepTools_dir"],
        program_parameters = cli_parameters_bamCoverage,
        extension = "150",
        ratio = "log2",
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        control = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{control}.bam",
        treatment = "{runID}/{outdir}/{reference_version}/bowtie2/merged/{treatment}.bam"
    output:
        "{runID}/{outdir}/{reference_version}/deepTools/bamCompare/{treatment}_vs_{control}.bw"
    shell:
        """
            {params.deepTools_dir}/bamCompare --bamfile1 {input.treatment} \
                                              --bamfile2 {input.control} \
                                              --outFileName {output} \
                                              --extendReads {params.extension} \
                                              --ratio {params.ratio} \
                                              --numberOfProcessors {threads} \
                                              --normalizeUsingRPKM \
                                              --ignoreForNormalization {params.ignore}
        """
