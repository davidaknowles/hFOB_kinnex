"""Compare satuRn DTU with the closed-form quasi-binomial isoform-usage test."""
import sys
import numpy as np, pandas as pd
from plotnine import *
J, OUT = sys.argv[1], sys.argv[2]
s = pd.read_csv(f"{J}/dtu/saturn_day4_vs_day0.tsv.gz", sep="\t")
q = pd.read_csv(f"{OUT}/tables/diff_isoform_usage_qb.tsv.gz", sep="\t")
m = q.merge(s, on=["isoform", "gene"], how="inner")
m.to_csv(f"{OUT}/tables/dtu_qb_vs_saturn.tsv.gz", sep="\t", index=False)
t = m[m.testable & m.empirical_FDR.notna()]
print("isoforms compared", len(t))
print("spearman logOR vs satuRn estimate", round(t[["logOR", "estimates"]].corr("spearman").iloc[0, 1], 3))
s_sig = (t.empirical_FDR < 0.05) & (t.dPSI.abs() >= 0.1)
r_sig = (t.regular_FDR < 0.05) & (t.dPSI.abs() >= 0.1)
print(pd.crosstab(t.sig, s_sig, rownames=["qb sig"], colnames=["satuRn empirical sig"]))
print(pd.crosstab(t.sig, r_sig, rownames=["qb sig"], colnames=["satuRn regular sig"]))
print("genes: qb", t[t.sig].gene.nunique(), "satuRn empirical", t[s_sig].gene.nunique(), "both", t[t.sig & s_sig].gene.nunique())
p = (ggplot(t.sample(min(len(t), 30000), random_state=0), aes("logOR", "estimates")) + geom_point(size=0.2, alpha=0.3)
     + geom_abline(color="red") + theme_bw() + labs(x="quasi-binomial log OR", y="satuRn estimate"))
p.save(f"{OUT}/figures/dtu_qb_vs_saturn.png", width=4.5, height=4.5, dpi=150)
