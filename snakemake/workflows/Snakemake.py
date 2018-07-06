_author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-07-06"

from snakemake.exceptions import MissingInputException
import os

rule:
    version: 0.1

localrules:
    all, run_kallisto, run_STAR, run_htseq, run_cutadapt

home = os.environ['HOME']
WORFKLOWDIR = config["WORFKLOWDIR"]

wrapper_dir = home + "/Development/snakemake-wrappers/bio"

include_prefix = WORFKLOWDIR + "/snakemake/rules/"

include:
    include_prefix + "run_fastp.py"

# define global variables such as reference version of genome so that it can be accessed
# throughout the whole worfklow
REF_GENOME = config["references"]["genomes"][1]

{outdir}/{reference_version}/bowtie2/{unit}.bam
rule all:
    input:
        expand("{runID}/{outdir}/{reference_version}/bowtie2/{unit}.bam",
               runID = "run1",
               outdir = config["processed_dir"],
               reference_version = config["references"][REF_GENOME]["version"][0],
               unit = ["ASC-ad-D0-H3K9me2-1",
                       "ASC-ad-D0-H3K9me2-2",
                       "ASC-ad-D0-H3K9me3-1",
                       "ASC-ad-D0-H3K9me3-2",
                       "ASC-ad-D0-INP-LMNB-1",
                       "ASC-ad-D1-H3K9me2-1",
                       "ASC-ad-D1-H3K9me2-2",
                       "ASC-ad-D1-H3K9me3-1",
                       "ASC-ad-D1-H3K9me3-2",
                       "ASC-ad-D1-INP-LMNB-1",
                       "ASC-ad-D3-H3K9me2-1",
                       "ASC-ad-D3-H3K9me2-2",
                       "ASC-ad-D3-H3K9me3-1",
                       "ASC-ad-D3-H3K9me3-2",
                       "ASC-ad-D3-INP-LMNB-1"]),
