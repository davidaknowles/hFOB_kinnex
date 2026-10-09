"""Read-level QC helpers for PacBio Kinnex single-cell BAMs."""
import random
import pysam
import pandas as pd


def read_lengths(bam_path, n=200_000, require_unaligned=True):
    """Lengths of the first n records of an (unaligned) PacBio BAM."""
    out = []
    with pysam.AlignmentFile(bam_path, check_sq=False) as f:
        for i, r in enumerate(f.fetch(until_eof=True)):
            if i >= n:
                break
            out.append(r.query_length)
    return pd.Series(out, name="length")


def sample_mapped(bam_path, frac=0.005, seed=0, real_cells_only=True):
    """Random subsample of primary alignments: length, aligned span, exon count, cell barcode."""
    rng = random.Random(seed)
    rows = []
    with pysam.AlignmentFile(bam_path) as f:
        for r in f.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if rng.random() > frac:
                continue
            if real_cells_only and r.has_tag("rc") and r.get_tag("rc") != 1:
                continue
            n_intron = sum(1 for op, l in r.cigartuples if op == 3)
            rows.append((r.query_name, r.reference_name, r.reference_start, r.reference_end,
                         "-" if r.is_reverse else "+", r.query_length, n_intron + 1,
                         r.get_tag("CB") if r.has_tag("CB") else None))
    return pd.DataFrame(rows, columns=["read", "chrom", "start", "end", "strand", "length", "n_exons", "CB"])
