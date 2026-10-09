"""Annotation helpers: GTF parsing and BED12 export."""
import gzip
import re
import pandas as pd


def read_gtf_exons(path, feature="exon"):
    """Exon records of a GTF/GFF (GENCODE or isoseq/pigeon style) as a DataFrame."""
    op = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != feature:
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', p[8]))
            tags = re.findall(r'tag "([^"]*)"', p[8])
            rows.append((p[0], int(p[3]) - 1, int(p[4]), p[6], a.get("transcript_id"), a.get("gene_id"),
                         a.get("gene_name"), a.get("gene_type"), ",".join(tags)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "transcript_id", "gene_id",
                                       "gene_name", "gene_type", "tags"])


def exons_to_bed12(ex):
    """Collapse exon rows (0-based half-open) to one BED12 line per transcript."""
    out = []
    for tid, g in ex.sort_values(["transcript_id", "start"]).groupby("transcript_id", sort=False):
        s, e = g.start.min(), g.end.max()
        sizes = ",".join(map(str, g.end - g.start)) + ","
        starts = ",".join(map(str, g.start - s)) + ","
        out.append((g.chrom.iat[0], s, e, tid, 0, g.strand.iat[0], s, e, 0, len(g), sizes, starts))
    return pd.DataFrame(out)
