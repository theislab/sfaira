import anndata
import os
# import any packages you require for dataloading here. you can assume packages like numpy and pandas being available


# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
def load(data_dir, **kwargs):
    # replace my-data-file.h5ad with the filename you're loading
    fn = os.path.join(data_dir, "becker_crc.h5ad")
    # replace the simple data loading code below with the code required to load your data file(s)
    adata = anndata.read_h5ad(fn)
    adata.obs["age"] = adata.obs["age"].astype(str)
    return adata  # your load function needs to return an AnnData object
