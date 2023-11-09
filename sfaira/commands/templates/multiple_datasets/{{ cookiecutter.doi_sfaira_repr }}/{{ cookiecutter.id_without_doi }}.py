import anndata
import os
# import any packages you require for dataloading here. you can assume packages like numpy and pandas being available


# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
# the sample_fn argument will be automatically set by sfaira to each of the sample_fns provided in the yaml top section
def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    # replace the simple data loading code below with the code required to load your data file(s)
    adata = anndata.read_h5ad(fn)

    return adata  # your load function needs to return an AnnData object
