"""Extract read -> cell barcode table from the merged BAM."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from hfob_kinnex.counts import read_to_cell
d = read_to_cell(f"{sys.argv[1]}/merged.bam", f"{sys.argv[1]}/read_cell.parquet")
print(d.height, d["CB"].n_unique())
