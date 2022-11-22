import os
import gzip
import anndata
import shutil


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE156793_S3_gene_count.loom.gz")
    fn_tmp = os.path.join(os.path.expanduser("~"), "tmp.loom")
    with gzip.open(fn, 'rb') as f_in:
        with open(fn_tmp, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    adata = anndata.read_loom(fn_tmp)
    os.remove(fn_tmp)

    adata.obs['Development_day'] = adata.obs['Development_day'].astype(str) + '-days'

    return adata
