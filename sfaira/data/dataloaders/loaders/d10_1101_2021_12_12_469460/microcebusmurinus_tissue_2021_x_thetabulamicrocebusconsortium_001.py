import anndata
import os


# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
def load(data_dir,**kwargs):
    sample_fn="31777475"
    adata=anndata.read(data_dir+os.sep+sample_fn)
    return adata  # your load function needs to return an AnnData object


