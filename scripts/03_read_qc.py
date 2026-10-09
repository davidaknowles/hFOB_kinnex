"""Collect read-length samples per sample: HiFi array reads, S-reads, mapped dedup reads."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import glob
from hfob_kinnex.readqc import read_lengths, sample_mapped

sl, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
for day in ["day0", "day4"]:
    ccs = glob.glob(f"{sl}/{day}/1-CCS/*.hifi_reads.*.bam")[0]
    read_lengths(ccs, 200_000).to_frame().assign(sample=day, type="HiFi array read").to_parquet(f"{out}/{day}_hifi_len.parquet")
    read_lengths(f"{sl}/{day}/2-SLoutputs/segmented.bam", 1_000_000).to_frame().assign(sample=day, type="S-read").to_parquet(f"{out}/{day}_sread_len.parquet")
    sample_mapped(f"{sl}/{day}/2-SLoutputs/scisoseq.mapped.bam", 0.005).assign(sample=day).to_parquet(f"{out}/{day}_mapped_sample.parquet")
