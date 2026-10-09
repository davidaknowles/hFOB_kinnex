"""Alternative event definitions on an isoform set: SUPPA2 local events and TSS / polyA-site clusters."""
import numpy as np
import pandas as pd
import scipy.sparse as sp


def transcript_structure(ex):
    """Per-transcript summary from exon rows: TSS, TES, first and last exon coordinates (strand aware)."""
    ex = ex.sort_values(["transcript_id", "start"])
    g = ex.groupby("transcript_id", sort=False)
    t = g.agg(chrom=("chrom", "first"), strand=("strand", "first"), start=("start", "min"), end=("end", "max"),
              n_exons=("start", "size"))
    first_lo = g.nth(0).set_index("transcript_id")[["start", "end"]]
    last_lo = g.nth(-1).set_index("transcript_id")[["start", "end"]]
    plus = t.strand == "+"
    t["tss"] = np.where(plus, t.start, t.end)
    t["tes"] = np.where(plus, t.end, t.start)
    # first exon = 5'-most in transcript orientation
    t["fe_start"] = np.where(plus, first_lo.loc[t.index, "start"], last_lo.loc[t.index, "start"])
    t["fe_end"] = np.where(plus, first_lo.loc[t.index, "end"], last_lo.loc[t.index, "end"])
    t["le_start"] = np.where(plus, last_lo.loc[t.index, "start"], first_lo.loc[t.index, "start"])
    t["le_end"] = np.where(plus, last_lo.loc[t.index, "end"], first_lo.loc[t.index, "end"])
    return t


def cluster_1d(pos, window):
    """Single-linkage clusters of sorted positions with gaps <= window. Returns labels in input order."""
    order = np.argsort(pos)
    p = np.asarray(pos)[order]
    lab = np.concatenate([[0], np.cumsum(np.diff(p) > window)])
    out = np.empty(len(pos), int)
    out[order] = lab
    return out


def end_clusters(tx, gene, counts, which="tss", window=50):
    """Cluster transcript 5' (tss) or 3' (tes) ends within each gene.

    tx : transcript_structure output, gene : Series transcript -> gene, counts : Series transcript -> pooled reads.
    Returns per-transcript cluster id and per-cluster table with usage fraction and whether the cluster is a distinct
    terminal exon (no overlap with other clusters' terminal exons) or lies in a shared exon (tandem).
    """
    t = tx.join(gene.rename("gene")).join(counts.rename("count")).dropna(subset=["gene"])
    t["count"] = t["count"].fillna(0)
    ex_s, ex_e = ("fe_start", "fe_end") if which == "tss" else ("le_start", "le_end")
    labs = []
    for gname, d in t.groupby("gene"):
        lab = cluster_1d(d[which].values, window)
        labs.append(pd.Series([f"{gname}|{which}{l}" for l in lab], index=d.index))
    t["cluster"] = pd.concat(labs)
    cl = t.groupby("cluster").agg(gene=("gene", "first"), chrom=("chrom", "first"), strand=("strand", "first"),
                                  pos=(which, "median"), count=("count", "sum"),
                                  ex_start=(ex_s, "min"), ex_end=(ex_e, "max"), n_tx=("gene", "size"))
    cl["frac"] = cl["count"] / cl.groupby("gene")["count"].transform("sum").replace(0, np.nan)
    return t["cluster"], cl


def classify_terminal(cl, min_frac=0.1, min_count=10):
    """Among used clusters (frac >= min_frac, count >= min_count) of genes with >=2 such clusters, label each gene
    as having distinct terminal exons (alternative first / last exon) and/or tandem ends within one exon."""
    u = cl[(cl.frac >= min_frac) & (cl["count"] >= min_count)].copy()
    u = u[u.groupby("gene")["gene"].transform("size") >= 2]
    rows = []
    for gname, d in u.groupby("gene"):
        d = d.sort_values("ex_start")
        overl = (d.ex_start.values[1:] < d.ex_end.values[:-1])  # neighbouring terminal exons overlap
        rows.append((gname, len(d), bool((~overl).any()), bool(overl.any())))
    return pd.DataFrame(rows, columns=["gene", "n_clusters", "distinct_exons", "tandem"])


def read_ioe(path):
    """SUPPA2 .ioe file -> DataFrame with event_id, gene, type, inclusion and total transcript lists."""
    d = pd.read_csv(path, sep="\t")
    d["type"] = d.event_id.str.split(";").str[1].str.split(":").str[0]
    d["incl"] = d.alternative_transcripts.str.split(",")
    d["total"] = d.total_transcripts.str.split(",")
    return d


def incidence(lists, var_names):
    """Sparse transcripts x events incidence matrix from per-event transcript lists."""
    idx = pd.Index(var_names)
    r, c = [], []
    for j, l in enumerate(lists):
        k = idx.get_indexer(l)
        k = k[k >= 0]
        r.extend(k)
        c.extend([j] * len(k))
    return sp.csr_matrix((np.ones(len(r), np.float32), (r, c)), shape=(len(idx), len(lists)))


def group_incidence(labels, var_names):
    """Transcripts x groups indicator for a Series transcript -> group label; returns matrix and group names."""
    lab = labels.reindex(var_names)
    ok = lab.notna().values
    groups = pd.Index(sorted(lab[ok].unique()))
    M = sp.csr_matrix((np.ones(ok.sum(), np.float32), (np.where(ok)[0], groups.get_indexer(lab[ok]))),
                      shape=(len(var_names), len(groups)))
    return M, groups


def knn_smooth_ratio(conn, y, n, self_weight=1.0):
    """kNN-smoothed proportion: (A y) / (A n) with A = connectivities + self_weight * I. Returns (psi, smoothed n)."""
    A = conn + self_weight * sp.identity(conn.shape[0], format="csr")
    ys, ns = A @ y, A @ n
    return np.divide(ys, ns, out=np.full_like(ys, np.nan, dtype=float), where=ns > 0), ns
