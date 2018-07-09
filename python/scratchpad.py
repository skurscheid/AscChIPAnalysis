from snakemake.io import expand, glob_wildcards

expand("{runID}/{outdir}/{reference_version}/deepTools/plotFingerprint/fingerprints_duplicates_removed.png",
       "{runID}/{outdir}/{reference_version}/deepTools/plotPCA/PCA_readCounts.png"
       "{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.png",
       "{runID}/{outdir}/{reference_version}/deepTools/plotCorrelation/heatmap_SpearmanCorr_readCounts.tab",
       runID = "run1",
       outdir = config["processed_dir"],
       reference_version = REF_VERSION)
