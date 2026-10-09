"""Export QC-pass cells x multi-isoform-gene isoforms for satuRn."""
import sys, os
import numpy as np, pandas as pd, anndata as ad, scipy.io as sio
J, OUT = sys.argv[1], sys.argv[2]
iso = ad.read_h5ad(f"{J}/isoforms.h5ad")
qc = pd.read_csv(f"{OUT}/tables/cell_qc.tsv", sep="\t", index_col=0)
iso = iso[qc.index[qc.qc_pass]]
tot = np.asarray(iso.X.sum(0)).ravel()
ok = tot >= 10
niso = iso.var.gene[ok].value_counts()
ok &= iso.var.gene.map(niso).fillna(0).values >= 2
sub = iso[:, ok]
os.makedirs(f"{J}/dtu", exist_ok=True)
sio.mmwrite(f"{J}/dtu/counts.mtx", sub.X.T.tocoo())  # isoforms x cells
sub.var[["gene"]].rename_axis("isoform").to_csv(f"{J}/dtu/txinfo.tsv", sep="\t")
sub.obs[["sample"]].rename_axis("cell").to_csv(f"{J}/dtu/cells.tsv", sep="\t")
print(sub.shape)
