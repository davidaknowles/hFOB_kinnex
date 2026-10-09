"""Joint isoform and gene count matrices (AnnData) from merged BAM + collapse group + pigeon filtered classification."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import pandas as pd, polars as pl
from hfob_kinnex.counts import read_to_cell, read_group, build_anndata, aggregate_genes

J = sys.argv[1]
rc = pl.read_parquet(f"{J}/read_cell.parquet") if os.path.exists(f"{J}/read_cell.parquet") else read_to_cell(f"{J}/merged.bam", f"{J}/read_cell.parquet")
print("reads", rc.height, flush=True)
grp = read_group(f"{J}/joint.group.txt")
cls = pd.read_csv(f"{J}/joint_classification.filtered_lite_classification.txt", sep="\t", low_memory=False)
cls = cls.drop(columns=[c for c in ["ORF_seq", "cell_barcodes"] if c in cls.columns])
iso = build_anndata(rc, grp, cls)
iso.write_h5ad(f"{J}/isoforms.h5ad")
aggregate_genes(iso).write_h5ad(f"{J}/genes.h5ad")
print(iso)
