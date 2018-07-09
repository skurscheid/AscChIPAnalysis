__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-07-06"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.1

localrules:
    all, run_kallisto, run_STAR, run_htseq, run_cutadapt

# define global variables such as reference version of genome so that it can be accessed
# throughout the whole worfklow
REF_GENOME = config["references"]["genomes"][1]
REF_VERSION = config["references"][REF_GENOME]["version"][2]

home = os.environ['HOME']
WORKFLOWDIR = config["WORKFLOWDIR"]
wrapper_dir = home + "/Development/snakemake-wrappers/bio"
include_prefix = WORKFLOWDIR + "/snakemake/rules/"

include:
    include_prefix + "run_fastp.py"
include:
    include_prefix + "bam_processing.py"
include:
    include_prefix + "run_deepTools_QC.py"


rule all:
    input:
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
        # expand("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.final.bam.bai",
        #        runID = "run1",
        #        outdir = config["processed_dir"],
        #        reference_version = REF_VERSION,
        #        unit = ["ASC-ad-D0-H3K9me2-1_S7",
        #                 "ASC-ad-D0-H3K9me2-2_S8",
        #                 "ASC-ad-D0-H3K9me3-1_S7",
        #                 "ASC-ad-D0-H3K9me3-2_S8",
        #                 "ASC-ad-D0-INP-LMNB-1_S6",
        #                 "ASC-ad-D1-H3K9me2-1_S1",
        #                 "ASC-ad-D1-H3K9me2-2_S2",
        #                 "ASC-ad-D1-H3K9me3-1_S3",
        #                 "ASC-ad-D1-H3K9me3-2_S4",
        #                 "ASC-ad-D1-INP-LMNB-1_S3",
        #                 "ASC-ad-D3-H3K9me2-1_S5",
        #                 "ASC-ad-D3-H3K9me2-2_S6",
        #                 "ASC-ad-D3-H3K9me3-1_S7",
        #                 "ASC-ad-D3-H3K9me3-2_S8",
        #                 "ASC-ad-D3-INP-LMNB-1_S6"])
