import anndata
import os

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This dataloader requires manual preprocessing of the Rdata file that can be obtained from the link in the
    `download_website` attribute of this class. The preprocessing code below uses the rpy2 and anndata2ri python
    packages to convert the R object to anndata (pip install anndata2ri), run it in a jupyter notebook:

    ## Notebook Cell 1
    import anndata2ri
    anndata2ri.activate()
    %load_ext rpy2.ipython

    ## Notebook Cell 2
    %%R -o sce
    library(Seurat)
    load("tissue.rdata")
    new_obj = CreateSeuratObject(counts = tissue@raw.data)
    new_obj@meta.data = tissue@meta.data
    sce <- as.SingleCellExperiment(new_obj)

    ## Notebook cell 3
    sce.write("ramachandran.h5ad")

    :param data_path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://datashare.is.ed.ac.uk/bitstream/handle/10283/3433/tissue.rdata"
        self.download_url_meta = None

        self.author = "Ramachandran"
        self.doi = "10.1038/s41586-019-1631-3"
        self.normalization = "raw"
        self.organ = "liver"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.year = 2019

        self.var_symbol_col = "index"

        self.obs_key_cellontology_original = "annotation_lineage"
        self.obs_key_state_exact = "condition"
        self.obs_key_healthy = self.obs_key_state_exact
        self.healthy_state_healthy = "Uninjured"

        self.set_dataset_id(idx=1)

    def _load(self):
        fn = os.path.join(self.data_dir, "ramachandran.h5ad")
        adata = anndata.read(fn)

        return adata
