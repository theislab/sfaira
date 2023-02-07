import os
import numpy as np
import pandas as pd
import gzip
import scipy
import tarfile
import anndata as ad
from scipy.io import mmread


# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
# the sample_fn argument will be automatically set by sfaira to each of the sample_fns provided in the yaml top section
def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, "GSE164102_RAW.tar")

    age_dict = {'GSM4996689_DGE67anew.txt.gz': '20',
                'GSM4996690_DGE67bnew.txt.gz': '20',
                'GSM4996691_DGE69new.txt.gz': '20',
                'GSM4996692_DGE68new.txt.gz': '20',
                'GSM4996693_DGE70new.txt.gz': '20',
                'GSM4996694_DGE66new.txt.gz': '20',
                'GSM4996695_DGE71new.txt.gz': '20',
                'GSM4996696_DGE7374.txt.gz': '40',
                'GSM4996697_DGE7778.txt.gz': '40',
                'GSM4996698_DGE7576.txt.gz': '40',
                'GSM4996699_DGE7980.txt.gz': '40',
                'GSM4996700_DGE72.txt.gz': '40',
                'GSM4996701_DGE81.txt.gz': '40'}

    with tarfile.open(fn) as tar:
        with gzip.open(tar.extractfile(sample_fn), 'rb') as df:
            d = pd.read_csv(df, delimiter='\t', index_col=0).T

    adata = ad.AnnData(X=scipy.sparse.csr_matrix(d), var=pd.DataFrame(index=d.columns.values),
                       obs=pd.DataFrame(index=d.index))

    adata.obs['organoid_age_days'] = age_dict[sample_fn]

    return adata
