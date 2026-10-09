"""BED12 of GENCODE Ensembl_canonical protein-coding transcripts (for gene-body coverage)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from hfob_kinnex.annot import read_gtf_exons, exons_to_bed12
ex = read_gtf_exons(sys.argv[1])
ex = ex[(ex.gene_type == "protein_coding") & ex.tags.str.contains("Ensembl_canonical") & ~ex.chrom.isin(["chrM", "chrY"])]
exons_to_bed12(ex).to_csv(sys.argv[2], sep="\t", header=False, index=False)
