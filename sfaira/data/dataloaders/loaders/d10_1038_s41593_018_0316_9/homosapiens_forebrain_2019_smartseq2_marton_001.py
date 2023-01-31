import anndata as ad
import os
import numpy as np
import pandas as pd
import scipy.sparse


# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
def load(data_dir, **kwargs):

    fn = os.path.join(data_dir, "GSE115011_marton_all_cells.csv.gz")
    d = pd.read_csv(fn, index_col = 0).T

    adata = ad.AnnData(X=scipy.sparse.csr_matrix(d),
                       var=pd.DataFrame(index=d.columns.values),
                       obs=pd.DataFrame(index=d.index))

    return adata
