import anndata
import os


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE79363_first_dataset_read_count.txt.gz")
    adata = anndata.read_text(fn).transpose()
    return adata
