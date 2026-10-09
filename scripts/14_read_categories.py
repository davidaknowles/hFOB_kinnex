"""Per-read breakdown (QC-pass cells) of pigeon filter outcome and structural category / subcategory by day."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import pandas as pd, polars as pl
from plotnine import *
from hfob_kinnex.counts import read_group

J, OUT = sys.argv[1], sys.argv[2]
cls = pl.read_csv(f"{J}/joint_classification.txt", separator="\t", columns=["isoform", "structural_category", "subcategory"],
                  infer_schema_length=0)
reasons = pl.read_csv(f"{J}/joint_classification.filtered_lite_reasons.txt", comment_prefix="#").rename(
    {"filtered_isoform": "isoform", "filter": "filter"})
cls = cls.join(reasons, on="isoform", how="left").with_columns(pl.col("filter").fill_null("pass"))
grp = read_group(f"{J}/joint.group.txt").with_columns(pl.col("read").str.split("_").list.first().alias("sample"))
rc = pl.read_parquet(f"{J}/read_cell.parquet")
qc = pd.read_csv(f"{OUT}/tables/cell_qc.tsv", sep="\t", index_col=0)
rc = rc.filter(pl.col("CB").is_in(list(qc.index[qc.qc_pass])))
d = grp.join(rc, on="read", how="inner").join(cls, left_on="pbid", right_on="isoform", how="left")
t = d.group_by(["sample", "filter", "structural_category", "subcategory"]).len().to_pandas()
t.to_csv(f"{OUT}/tables/read_categories.tsv", sep="\t", index=False)
f = t.groupby(["sample", "filter"])["len"].sum().unstack(); print((f.T / f.sum(1)).round(3))
p = t[t["filter"] == "pass"].groupby(["sample", "structural_category", "subcategory"])["len"].sum().reset_index()
p["frac"] = p["len"] / p.groupby("sample")["len"].transform("sum")
print(p.sort_values("frac", ascending=False).head(25).to_string())
# 5'/3' completeness among reads assigned to known transcripts
k = p[p.structural_category.isin(["full-splice_match", "incomplete-splice_match"])]
print(k.pivot_table(index="subcategory", columns="sample", values="frac").round(4))
a = t.groupby(["sample", "filter"])["len"].sum().reset_index()
a["frac"] = a["len"] / a.groupby("sample")["len"].transform("sum")
g = (ggplot(a, aes("sample", "frac", fill="filter")) + geom_col() + theme_bw() + labs(y="fraction of reads", x=""))
g.save(f"{OUT}/figures/read_filter_reasons.png", width=4.5, height=3.5, dpi=150)
c = p.groupby(["sample", "structural_category"])["frac"].sum().reset_index()
g = (ggplot(c, aes("sample", "frac", fill="structural_category")) + geom_col() + theme_bw() + labs(y="fraction of passing reads", x=""))
g.save(f"{OUT}/figures/read_structural_category.png", width=5.5, height=3.5, dpi=150)
