"""Build cell x isoform and cell x gene count matrices from a joint isoseq collapse."""
import subprocess
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp
import anndata as ad


def read_to_cell(bam, out_parquet, threads=8):
    """Table of read name -> CB tag from an aligned BAM (primary alignments only)."""
    cmd = (f"samtools view -@ {threads} -F 0x904 {bam} | "
           "awk -F'\\t' '{for(i=12;i<=NF;i++) if(substr($i,1,5)==\"CB:Z:\"){print $1\"\\t\"substr($i,6); break}}'")
    p = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE)
    df = pl.read_csv(p.stdout, separator="\t", has_header=False, new_columns=["read", "CB"])
    df.write_parquet(out_parquet)
    return df


def read_group(group_txt):
    """isoseq collapse group file -> long table (pbid, read)."""
    g = pl.read_csv(group_txt, separator="\t", has_header=False, new_columns=["pbid", "reads"])
    return g.with_columns(pl.col("reads").str.split(",")).explode("reads").rename({"reads": "read"})


def build_anndata(read_cell, group_long, classification, gene_col="associated_gene"):
    """Isoform-level AnnData (cells x isoforms) restricted to isoforms in `classification`."""
    cls = classification.set_index("isoform")
    j = group_long.join(read_cell, on="read", how="inner").filter(pl.col("pbid").is_in(list(cls.index)))
    cnt = j.group_by(["CB", "pbid"]).len().to_pandas()
    cells = pd.Index(sorted(cnt.CB.unique()))
    isos = pd.Index(cls.index[cls.index.isin(cnt.pbid.unique())])
    X = sp.csr_matrix((cnt["len"].to_numpy(np.float32), (cells.get_indexer(cnt.CB), isos.get_indexer(cnt.pbid))),
                      shape=(len(cells), len(isos)))
    var = cls.loc[isos].copy()
    var["gene"] = var[gene_col].astype(str)
    obs = pd.DataFrame(index=cells)
    obs["sample"] = [c.split("_")[0] for c in cells]
    return ad.AnnData(X=X, obs=obs, var=var)


def aggregate_genes(iso):
    """Sum isoform counts to genes."""
    genes = pd.Index(sorted(iso.var.gene.unique()))
    M = sp.csr_matrix((np.ones(iso.n_vars, np.float32), (np.arange(iso.n_vars), genes.get_indexer(iso.var.gene))),
                      shape=(iso.n_vars, len(genes)))
    g = ad.AnnData(X=(iso.X @ M).tocsr(), obs=iso.obs.copy(), var=pd.DataFrame(index=genes))
    return g
