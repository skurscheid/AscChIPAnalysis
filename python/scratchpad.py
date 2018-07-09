from snakemake.io import expand, glob_wildcards
import json
from pprint import pprint

REF_GENOME = config["references"]["genomes"][1]
REF_VERSION = config["references"][REF_GENOME]["version"][2]

wildcards = {"assayID" : "ChIP-Seq",
             "runID" : "run1",
             "outdir" : "processed_data",
             "reference_version": REF_VERSION,
             "application" : "bowtie2",
             "sample" : "ASC-ad-D0-H3K9me2"}

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
