import anndata
import os
import scipy.sparse
import pandas as pd
import gzip
import scipy.io


def load(data_dir, **kwargs):
    fn = [os.path.join(data_dir, "GSE158055_covid19_counts.mtx.gz"),
          os.path.join(data_dir, "GSE158055_covid19_barcodes.tsv.gz"),
          os.path.join(data_dir, "GSE158055_covid19_features.tsv.gz")]
    fn_meta = [os.path.join(data_dir, "GSE158055_cell_annotation.csv.gz"),
               os.path.join(data_dir, "GSE158055_sample_metadata.xlsx")]

    # Data is not in sc.read_10x_mtx format because GSE158055_covid19_counts does not end on "_matrix".
    with gzip.open(fn[0], 'rb') as mm:
        X = scipy.io.mmread(mm).T.tocsr()
    obs = pd.read_csv(fn[1], header=None, sep="\t", index_col=0)
    obs.index.name = None
    var = pd.read_csv(fn[2], header=None, sep="\t", names=['names'])
    var.index = var['names'].values
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    annotation = pd.read_csv(fn_meta[0])
    metadata = pd.read_excel(fn_meta[1], skiprows=20, skipfooter=37)
    metadata.columns = metadata.columns.str.split('characteristics: ').str[-1]
    tmp = annotation.merge(metadata, left_on='sampleID', right_on='Sample name').set_index('cellName')
    adata.obs = adata.obs.join(tmp).astype(str)

    return adata
