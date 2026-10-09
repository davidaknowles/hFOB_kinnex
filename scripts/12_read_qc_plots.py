"""Read-length distributions (HiFi array reads, S-reads, deduplicated mapped reads) and exon-count summaries."""
import sys, glob
import numpy as np, pandas as pd
from plotnine import *
Q, OUT = sys.argv[1], sys.argv[2]
L = pd.concat([pd.read_parquet(f) for f in glob.glob(f"{Q}/*_len.parquet")])
M = pd.concat([pd.read_parquet(f) for f in glob.glob(f"{Q}/*_mapped_sample.parquet")])
L = pd.concat([L, M[["length", "sample"]].assign(type="mapped read (real cells)")])
s = L.groupby(["type", "sample"]).length.describe(percentiles=[0.1, 0.5, 0.9]).round(0)
s.to_csv(f"{OUT}/tables/read_length_summary.tsv", sep="\t"); print(s)
p = (ggplot(L, aes("length", color="sample")) + geom_density() + facet_wrap("~type", scales="free")
     + scale_x_log10() + theme_bw() + labs(x="read length (bp)"))
p.save(f"{OUT}/figures/read_lengths.png", width=12, height=3.5, dpi=150)
e = M.assign(n_exons=M.n_exons.clip(upper=20)).groupby(["sample", "n_exons"]).size().rename("n").reset_index()
e["frac"] = e.n / e.groupby("sample").n.transform("sum")
p = (ggplot(e, aes("n_exons", "frac", fill="sample")) + geom_col(position="dodge") + theme_bw()
     + labs(x="exons per read (20 = 20+)", y="fraction of reads"))
p.save(f"{OUT}/figures/exons_per_read.png", width=6, height=3.5, dpi=150)
print(M.groupby("sample").agg(mono_exon=("n_exons", lambda x: (x == 1).mean()), median_exons=("n_exons", "median")))
