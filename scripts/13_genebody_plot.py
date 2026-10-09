"""Plot RSeQC gene-body coverage (5'->3' percentile) per sample, normalised to max."""
import sys, pandas as pd
from plotnine import *
d = pd.read_csv(sys.argv[1], sep="\t", index_col=0).T
d.index = d.index.astype(int); d.columns = [c.replace(".sub", "") for c in d.columns]
n = (d / d.max()).reset_index(names="percentile").melt(id_vars="percentile", var_name="sample", value_name="coverage")
print("coverage at 5' 10th pct / 3' 90th pct:", (d.loc[10] / d.loc[90]).round(2).to_dict())
p = (ggplot(n, aes("percentile", "coverage", color="sample")) + geom_line() + theme_bw()
     + labs(x="gene body percentile (5' → 3')", y="relative coverage", title="Canonical protein-coding transcripts"))
p.save(sys.argv[2], width=6, height=3.5, dpi=150)
