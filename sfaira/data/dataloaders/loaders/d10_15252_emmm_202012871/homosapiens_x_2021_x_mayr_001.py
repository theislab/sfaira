import anndata
import os
import scipy.sparse

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    dir_tar = buffered_decompress(os.path.join(data_dir, "mayr_et_al.tar.gz"))
    fn = os.path.join(dir_tar, "storage", "groups", "ml01", "datasets", "projects",
                      "20200210_Schiller_dropseq_humanILD_meshal.ansari", "data_objects", sample_fn)
    adata = anndata.read(fn)
    if not isinstance(adata.X, scipy.sparse.csr_matrix):
        adata.X = scipy.sparse.csr_matrix(adata.X)
    return adata
