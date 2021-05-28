import anndata
import gzip
import os
import pandas as pd
import scipy.io
import tarfile


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE164378_RAW.tar")
    adatas = []
    with tarfile.open(fn) as tar:
        samples = ['GSM5008737_RNA_3P', 'GSM5008738_ADT_3P']
        for sample in samples:
            with gzip.open(tar.extractfile(sample + '-matrix.mtx.gz'), 'rb') as mm:
                x = scipy.io.mmread(mm).T.tocsr()
            obs = pd.read_csv(tar.extractfile(sample + '-barcodes.tsv.gz'), compression='gzip',
                              header=None, sep='\t', index_col=0)
            obs.index.name = None
            var = pd.read_csv(tar.extractfile(sample + '-features.tsv.gz'), compression='gzip',
                              header=None, sep='\t').iloc[:, :1]
            var.columns = ['names']
            var.index = var['names'].values
            adata = anndata.AnnData(X=x, obs=obs, var=var)
            adata.var_names_make_unique()
            adatas.append(adata)
        tar.close()

    adata = adatas[0]
    protein = adatas[1]
    meta = pd.read_csv(os.path.join(data_dir, 'GSE164378_sc.meta.data_3P.csv.gz'), index_col=0)
    adata.obs = adata.obs.join(meta)
    adata.obsm['protein_expression'] = pd.DataFrame(protein.X.A, columns=protein.var_names, index=protein.obs_names)

    return adata
