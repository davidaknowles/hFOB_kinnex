"""Isoform structure + usage plots for example genes with differential isoform usage."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd, anndata as ad
from hfob_kinnex.annot import read_gtf_exons
from hfob_kinnex.plots import isoform_usage_plot

J, OUT = sys.argv[1], sys.argv[2]
genes = sys.argv[3].split(";")
iso = ad.read_h5ad(f"{J}/isoforms.h5ad")
qc = pd.read_csv(f"{OUT}/tables/cell_qc.tsv", sep="\t", index_col=0)
iso = iso[qc.index[qc.qc_pass]]
gene = iso.var.gene.astype(str)
ex = read_gtf_exons(f"{J}/events/expressed.gtf")
os.makedirs(f"{OUT}/figures/isoforms", exist_ok=True)
for g in genes:
    m = (gene == g).values
    sub = iso[:, m]
    pb = pd.DataFrame({s: np.asarray(sub.X[(sub.obs["sample"] == s).values].sum(0)).ravel() for s in ["day0", "day4"]},
                      index=sub.var_names)
    pb = pb[pb.index.isin(ex.transcript_id)]
    u = pb / pb.sum()
    top = u.max(1).sort_values(ascending=False).index[:8]
    u = u.loc[top]
    p = isoform_usage_plot(ex, u, f"{g}  (reads: day0 {int(pb.day0.sum())}, day4 {int(pb.day4.sum())})",
                           strand=sub.var.strand.iloc[0])
    p.save(f"{OUT}/figures/isoforms/{g}.png", width=10, height=1 + 0.45 * len(u), dpi=150)
