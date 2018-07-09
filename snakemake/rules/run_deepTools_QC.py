__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2017-01-11"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

from snakemake.exceptions import MissingInputException
import os

"""
Rules for running deepTools QC/QC on ChIP-Seq data
For usage, include this in your workflow.
"""

# samples for processing
UNITS = ["ASC-ad-D0-H3K9me2-1_S7",
         "ASC-ad-D0-H3K9me2-2_S8",
         "ASC-ad-D0-H3K9me3-1_S7",
         "ASC-ad-D0-H3K9me3-2_S8",
         "ASC-ad-D0-INP-LMNB-1_S6",
         "ASC-ad-D1-H3K9me2-1_S1",
         "ASC-ad-D1-H3K9me2-2_S2",
         "ASC-ad-D1-H3K9me3-1_S3",
         "ASC-ad-D1-H3K9me3-2_S4",
         "ASC-ad-D1-INP-LMNB-1_S3",
         "ASC-ad-D3-H3K9me2-1_S5",
         "ASC-ad-D3-H3K9me2-2_S6",
         "ASC-ad-D3-H3K9me3-1_S7",
         "ASC-ad-D3-H3K9me3-2_S8",
         "ASC-ad-D3-INP-LMNB-1_S6"])

# rules
rule plotFingerprint_deduplicated:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = "BAM fingerprint",
        labels = UNITS,
        extension = "200"
    threads:
        lambda wildcards: int(str(config["program_parameters"]["deepTools"]["threads"]).strip("['']"))
    input:
        bam = expand("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.bam",
                       runID = "run1",
                       outdir = config["processed_dir"],
                       reference_version = REF_VERSION,
                       unit = UNITS)
        index = expand("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.bam.bai",
                       runID = "run1",
                       outdir = config["processed_dir"],
                       reference_version = REF_VERSION,
                       unit = UNITS)
    output:
        "{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/fingerprints_duplicates_removed.png"
    shell:
        """
            {params.deepTools_dir}/plotFingerprint --bamfiles {input.bam} \
                                                   --numberOfProcessors {threads} \
                                                   --extendReads {params.extension}\
                                                   --plotTitle "{params.plotTitle}" \
                                                   --labels {params.labels} \
                                                   --skipZeros \
                                                   --plotFile {output}
        """

rule multiBamSummary_deduplicated:
    version:
        0.2
    params:
        deepTools_dir = home + config["deepTools_dir"],
        binSize = config["program_parameters"]["deepTools"]["binSize"],
        labels = UNITS
    threads:
        24
    input:
        bam = expand("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.bam",
                       runID = "run1",
                       outdir = config["processed_dir"],
                       reference_version = REF_VERSION,
                       unit = UNITS)
        index = expand("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.bam.bai",
                       runID = "run1",
                       outdir = config["processed_dir"],
                       reference_version = REF_VERSION,
                       unit = UNITS)
    output:
        npz = "{runID}/{outdir}/{reference_version}/deepTools/multiBamSummary/results.npz"
    shell:
        """
            {params.deepTools_dir}/multiBamSummary bins --bamfiles {input.bam} \
                                                        --labels {params.labels} \
                                                        --numberOfProcessors {threads} \
                                                        --centerReads \
                                                        --binSize {params.binSize} \
                                                        --outFileName {output.npz}
        """


rule plotCorrelation_heatmap:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = "Correlation heatmap"
    input:
        npz = rules.multiBamSummary_deduplicated.output
    output:
        "{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.png",
        "{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.tab"
    shell:
        """
            {params.deepTools_dir}/plotCorrelation --corData {input.npz} \
                                                   --corMethod spearman \
                                                   --skipZeros \
                                                   --plotTitle "{params.plotTitle}" \
                                                   --whatToPlot heatmap \
                                                   --colorMap RdYlBu \
                                                   --plotNumbers \
                                                   -o {output[0]} \
                                                   --outFileCorMatrix {output[1]}
        """

rule plotPCA:
    params:
        deepTools_dir = home + config["deepTools_dir"],
        plotTitle = "PCA - "
    input:
        npz = rules.multiBamSummary_deduplicated.output
    output:
        "{runID}/{outdir}/{reference_version}/deepTools/plotPCA/PCA_readCounts.png"
    shell:
        """
            {params.deepTools_dir}/plotPCA --corData {input.npz} \
                                           --plotFile {output} \
                                           --plotTitle "{params.plotTitle}"
        """
