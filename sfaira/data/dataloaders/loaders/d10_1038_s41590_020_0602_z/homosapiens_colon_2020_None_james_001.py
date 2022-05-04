import anndata
import os
import numpy as np
import scipy.sparse


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "james20.processed.h5ad")
    adata = anndata.read(fn)
    adata.X = np.expm1(adata.X)
    adata.X = adata.X.multiply(scipy.sparse.csc_matrix(adata.obs["n_counts"].values[:, None])).multiply(1 / 10000)
    # Assay maps are described here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7212050/
    adata.obs["assay"] = [{
        '290b': "10x 3' transcription profiling",
        '298c': "10x 3' transcription profiling",
        '302c': "10x 3' transcription profiling",
        '390c': "10x 5' transcription profiling",
        '417c': "10x 5' transcription profiling",
    }[x] for x in adata.obs["donor"].values]
    return adata
