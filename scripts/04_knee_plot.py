"""Barcode-rank (knee) plots from SMRT Link bcstats reports."""
import sys, pandas as pd
from plotnine import *
sl, out = sys.argv[1], sys.argv[2]
dfs = []
for day in ["day0", "day4"]:
    d = pd.read_csv(f"{sl}/{day}/2-SLoutputs/bcstats_report.tsv.gz", sep="\t")
    d.columns = ["barcode", "reads", "rank", "molecules", "real_cell", "type"]
    d = d[d.type.str.startswith("Group")].sort_values("molecules", ascending=False)
    d["rank"] = range(1, len(d) + 1)
    d["sample"] = day
    dfs.append(d[d["rank"] <= 200_000])
    print(day, "cells:", (d.real_cell == "cell").sum(), "median UMI in cells:", d.loc[d.real_cell == "cell", "molecules"].median())
d = pd.concat(dfs)
p = (ggplot(d, aes("rank", "molecules", color="real_cell")) + geom_point(size=0.3) + facet_wrap("~sample")
     + scale_x_log10() + scale_y_log10() + labs(x="barcode rank", y="UMIs (deduplicated molecules)", color="")
     + theme_bw())
p.save(out, width=8, height=3.5, dpi=150)
