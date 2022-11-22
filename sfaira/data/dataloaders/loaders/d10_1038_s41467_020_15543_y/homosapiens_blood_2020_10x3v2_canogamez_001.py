import anndata
import os
import pandas as pd
import scipy.io


def load(data_dir, **kwargs):
    fn = [os.path.join(data_dir, "NCOMMS-19-7936188_scRNAseq_raw_UMIs.mtx"),
          os.path.join(data_dir, 'NCOMMS-19-7936188_scRNAseq_barcodes.tsv'),
          os.path.join(data_dir, 'NCOMMS-19-7936188_scRNAseq_genes.tsv'),
          os.path.join(data_dir, "NCOMMS-19-7936188_metadata.txt")]

    with open(fn[0], 'rb') as mm:
        X = scipy.io.mmread(mm).T.tocsr()
    obs = pd.read_csv(fn[1], header=None, sep="\t", index_col=0)
    obs.index.name = None
    var = pd.read_csv(fn[2], header=None, sep="\t", names=['names'])
    var.index = var['names'].values

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    metadata = pd.read_csv(fn[3], sep="\t")
    adata.obs = adata.obs.join(metadata, how="inner")
    return adata
