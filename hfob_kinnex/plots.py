"""Isoform structure plots with per-group usage."""
import numpy as np
import pandas as pd
from plotnine import *


def compress_coords(exons, intron_len=150):
    """Map genomic coordinates to a compressed axis where gaps between merged exon blocks have fixed length."""
    blocks = exons[["start", "end"]].sort_values("start").values
    merged = []
    for s, e in blocks:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    offs, pos = [], 0
    for s, e in merged:
        offs.append((s, e, pos))
        pos += (e - s) + intron_len
    def f(x):
        for s, e, o in offs:
            if x <= e:
                return o + max(x - s, 0)
        return offs[-1][2] + offs[-1][1] - offs[-1][0]
    return np.vectorize(f)


def isoform_usage_plot(exons, usage, title, strand="+"):
    """exons: transcript_id, start, end. usage: DataFrame indexed by transcript with one column per group (fractions)."""
    ex = exons[exons.transcript_id.isin(usage.index)].copy()
    f = compress_coords(ex)
    ex["x0"], ex["x1"] = f(ex.start.values), f(ex.end.values)
    lab = {t: f"{t}  " + "  ".join(f"{c} {usage.loc[t, c]:.0%}" for c in usage.columns) for t in usage.index}
    order = list(usage.index)[::-1]
    ex["y"] = ex.transcript_id.map({t: i for i, t in enumerate(order)})
    tx = ex.groupby("transcript_id").agg(x0=("x0", "min"), x1=("x1", "max"), y=("y", "first")).reset_index()
    diff = usage.iloc[:, -1] - usage.iloc[:, 0]
    ex["dUsage"] = ex.transcript_id.map(diff)
    p = (ggplot() + geom_segment(tx, aes(x="x0", xend="x1", y="y", yend="y"), color="grey")
         + geom_rect(ex, aes(xmin="x0", xmax="x1", ymin="y-0.3", ymax="y+0.3", fill="dUsage"))
         + scale_fill_gradient2(low="#2e86c1", mid="#dddddd", high="#c0392b", midpoint=0, limits=(-1, 1),
                                name="Δ usage\n(day4−day0)")
         + scale_y_continuous(breaks=list(range(len(order))), labels=[lab[t] for t in order])
         + labs(x=f"compressed coordinate ({'5′→3′' if strand == '+' else '3′←5′'})", y="", title=title)
         + theme_bw() + theme(panel_grid=element_blank(), axis_text_x=element_blank()))
    return p
