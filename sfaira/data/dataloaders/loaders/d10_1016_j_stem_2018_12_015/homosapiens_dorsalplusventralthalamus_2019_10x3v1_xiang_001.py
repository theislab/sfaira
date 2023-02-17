import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import os
import scipy.io


def load(data_dir, **kwargs):
    fn_matrix = os.path.join(data_dir, "GSE122342_matrix.mtx.gz")
    fn_genes = os.path.join(data_dir, "GSE122342_genes.tsv.gz")
    fn_barcodes = os.path.join(data_dir, "GSE122342_barcodes.tsv.gz")
    fn_schThO_class = os.path.join(data_dir, "schThO_class.txt")
    fn_GSEA_enrich = os.path.join(data_dir, "GSEA_category_enrich.txt")

    barcodes = pd.read_csv(fn_barcodes, delimiter='\t', compression='gzip', index_col=0, header=None, names=[None])
    metadata1 = pd.read_csv(fn_schThO_class, delimiter='\t', index_col=0)
    metadata1["organoid_age_days"] = metadata1["Time"].replace({"early": "34", "late": "89"}).astype("category")
    metadata2 = pd.read_csv(fn_GSEA_enrich, delimiter='\t', index_col=0)
    obs = barcodes.join(metadata1).join(metadata2)

    x = scipy.io.mmread(fn_matrix).T.tocsr().astype(np.float32)
    var = pd.read_csv(fn_genes, delimiter='\t', compression='gzip', index_col=0, header=None, names=[None, "gene_id"])

    adata = ad.AnnData(
        X=x,
        obs=obs,
        var=var,
    )

    return adata
