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

# includes
include:
    include_prefix + "pooled_replicates_processing.py"
include:
    include_prefix + "macs2_processing.py"
include:
    include_prefix + "run_edd.py"

SAMPLES = ["ASC-ad-D0-H3K9me2",
          "ASC-ad-D0-H3K9me3",
          "ASC-ad-D0-INP-LMNB",
          "ASC-ad-D1-H3K9me2",
          "ASC-ad-D1-H3K9me3",
          "ASC-ad-D1-INP-LMNB",
          "ASC-ad-D3-H3K9me2",
          "ASC-ad-D3-H3K9me3",
          "ASC-ad-D3-INP-LMNB"]

rule all:
    input:
        expand("{runID}/{outdir}/{reference_version}/deepTools/bamCompare/{treatment}_vs_{control}.bw",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               treatment = SAMPLES[0:2],
               control = SAMPLES[2]),
        expand("{runID}/{outdir}/{reference_version}/deepTools/bamCompare/{treatment}_vs_{control}.bw",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               treatment = SAMPLES[3:5],
               control = SAMPLES[5]),
        expand("{runID}/{outdir}/{reference_version}/deepTools/bamCompare/{treatment}_vs_{control}.bw",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               treatment = SAMPLES[6:8],
               control = SAMPLES[8]),
        expand("{runID}/{outdir}/{reference_version}/macs2/predictd/{ChIP}",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               ChIP = SAMPLES[0:2] + SAMPLES[3:5] + SAMPLES[6:8]),
        # callpeak - cutoff analysis
        expand("{runID}/{outdir}/{reference_version}/macs2/callpeak/{ChIP}_vs_{Input}/",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               ChIP = SAMPLES[0:2],
               Input = SAMPLES[2]),
        expand("{runID}/{outdir}/{reference_version}/macs2/callpeak/{ChIP}_vs_{Input}/",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               ChIP = SAMPLES[3:5],
               Input = SAMPLES[5]),
        expand("{runID}/{outdir}/{reference_version}/macs2/callpeak/{ChIP}_vs_{Input}/",
           runID = "run1",
           outdir = config["processed_dir"],
           reference_version = REF_VERSION,
           ChIP = SAMPLES[6:8],
           Input = SAMPLES[8]),
        # callpeak, broad peaks
        expand("{runID}/{outdir}/{reference_version}/macs2/callpeak_broad/{ChIP}_vs_{Input}/",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               ChIP = SAMPLES[0:2],
               Input = SAMPLES[2]),
        expand("{runID}/{outdir}/{reference_version}/macs2/callpeak_broad/{ChIP}_vs_{Input}/",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               ChIP = SAMPLES[3:5],
               Input = SAMPLES[5]),
        expand("{runID}/{outdir}/{reference_version}/macs2/callpeak_broad/{ChIP}_vs_{Input}/",
           runID = "run1",
           outdir = config["processed_dir"],
           reference_version = REF_VERSION,
           ChIP = SAMPLES[6:8],
           Input = SAMPLES[8]),
        # run EDD
        expand("{runID}/{outdir}/{reference_version}/edd/{ChIP}_vs_{Input}/",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               ChIP = SAMPLES[0:2],
               Input = SAMPLES[2]),
        expand("{runID}/{outdir}/{reference_version}/edd/{ChIP}_vs_{Input}/",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = REF_VERSION,
               ChIP = SAMPLES[3:5],
               Input = SAMPLES[5]),
        expand("{runID}/{outdir}/{reference_version}/edd/{ChIP}_vs_{Input}/",
           runID = "run1",
           outdir = config["processed_dir"],
           reference_version = REF_VERSION,
           ChIP = SAMPLES[6:8],
           Input = SAMPLES[8])
