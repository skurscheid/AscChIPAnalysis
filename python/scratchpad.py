from snakemake.io import expand, glob_wildcards
import json
from pprint import pprint
import os

REF_GENOME = config["references"]["genomes"][1]
REF_VERSION = config["references"][REF_GENOME]["version"][2]

chromSizes = home + config["references"][REF_GENOME]["chromSizesUCSC"],
gaps = home + config["references"][REF_GENOME]["gapsUCSC"],
binSize = "", # don't set for autoestimation
numTrials = "10000", #number of Monte Carlo trials 10000 default
fdr = "0.1",
gapPenalty = "", # using default - autodetect
eddBinDir = home + config["edd_dir"]

wildcards = {"assayID" : "ChIP-Seq",
             "runID" : "run1",
             "outdir" : "processed_data",
             "reference_version": REF_VERSION,
             "application" : "bowtie2",
             "sample" : "ASC-ad-D0-H3K9me2",
             "type" : "replicates"}

bam = lambda wildcards: expand("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.{suffix}",
                         outdir = wildcards["outdir"],
                         reference_version = wildcards["reference_version"],
                         unit = config["samples"]["ChIP-Seq"][wildcards["type"]][wildcards["sample"]],
                         suffix = ["bam"])

expand("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.{suffix}",
        runID = wildcards["runID"],
        outdir = wildcards["outdir"],
        reference_version = wildcards["reference_version"],
        unit = config["samples"]["ChIP-Seq"][wildcards["type"]][wildcards["sample"]],
        suffix = ["bam"])

expand("{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/fingerprints_duplicates_removed.png",
       runID = "run1",
       outdir = config["processed_dir"],
       reference_version = REF_VERSION),
expand("{runID}/{outdir}/{reference_version}/deepTools/plotPCA/PCA_readCounts.png",
        runID = "run1",
       outdir = config["processed_dir"],
       reference_version = REF_VERSION),
expand("{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.{suffix}",
        runID = "run1",
        outdir = config["processed_dir"],
        reference_version = REF_VERSION,
        suffix = ["png", "tab"])

expand("{runID}/{outdir}/{reference_version}/deepTools/bamCompare/{treatment}_vs_{control}.bw",
       runID = "run1",
       outdir = config["processed_dir"],
       reference_version = REF_VERSION,
       treatment = SAMPLES[0:2],
       control = SAMPLES[2])

               control = SAMPLES[2]),
expand("{runID}/{outdir}/{reference_version}/deepTools/bamCompare/{treatment}_vs_{control}.bw",
       runID = "run1",
       outdir = config["processed_dir"],
       reference_version = REF_VERSION,
       treatment = SAMPLES[3:5],
       control = SAMPLES[5])

expand("{runID}/{outdir}/{reference_version}/deepTools/bamCompare/{treatment}_vs_{control}.bw",
       runID = "run1",
       outdir = config["processed_dir"],
       reference_version = REF_VERSION,
       treatment = SAMPLES[6:8],
       control = SAMPLES[8])
