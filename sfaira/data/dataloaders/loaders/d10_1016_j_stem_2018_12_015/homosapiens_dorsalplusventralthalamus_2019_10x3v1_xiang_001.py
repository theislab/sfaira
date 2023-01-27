import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import os
import scipy
import gzip
from tempfile import TemporaryDirectory
import shutil


def load(data_dir, **kwargs):
    fn_matrix = os.path.join(data_dir, "GSE122342_matrix.mtx.gz")
    fn_genes = os.path.join(data_dir, "GSE122342_genes.tsv.gz")
    fn_barcodes = os.path.join(data_dir, "GSE122342_barcodes.tsv.gz")
    fn_schThO_class = os.path.join(data_dir, "schThO_class.txt")  # annotated cell clusters 
    fn_GSEA_enrich = os.path.join(data_dir, "GSEA_category_enrich.txt")  # gene set enrichment analysis (GSEA) 

    uncompressed_file_type = "mtx"
    with TemporaryDirectory() as tmpdir:
        tmppth = tmpdir + f"/decompressed.{uncompressed_file_type}"
        with gzip.open(fn_matrix, "rb") as input_f, open(tmppth, "wb") as output_f:
            shutil.copyfileobj(input_f, output_f)
        matrix = scipy.io.mmread(fn_matrix)

    genes = pd.read_csv(fn_genes, delimiter='\t', compression='gzip', index_col=0, header = None, names=["ensembl_gene_id", "gene_id"])
    barcodes = pd.read_csv(fn_barcodes, delimiter='\t', compression='gzip', index_col=0, header=None)
    barcodes.index.name = None
    metadata1 = pd.read_csv(fn_schThO_class, delimiter='\t', index_col=0)
    metadata1["organoid_age_days"] = "34"
    metadata1.loc[metadata1.Time == "late", "organoid_age_days"] = "89"
    obs = barcodes.join(metadata1)
    metadata2 = pd.read_csv(fn_GSEA_enrich, delimiter='\t', index_col=0)
    obs = obs.join(metadata2)

    adata = ad.AnnData(
        X=matrix.T,
        obs=obs,
        var=genes,)
    adata

    return adata
