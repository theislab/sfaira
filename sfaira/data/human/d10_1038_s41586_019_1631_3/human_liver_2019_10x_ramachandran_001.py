import anndata
import os
from typing import Union

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
    load('tissue.rdata')
    new_obj = CreateSeuratObject(counts = tissue@raw.data)
    new_obj@meta.data = tissue@meta.data
    sce <- as.SingleCellExperiment(new_obj)

    ## Notebook cell 3
    sce.write('ramachandran.h5ad')

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.organism = "human"
        self.id = "human_liver_2019_10x_ramachandran_001_10.1038/s41586-019-1631-3"
        self.download = "https://datashare.is.ed.ac.uk/bitstream/handle/10283/3433/tissue.rdata"
        self.download_meta = None
        self.organ = "liver"
        self.sub_tissue = "liver"
        self.author = 'Henderson'
        self.year = 2019
        self.doi = '10.1038/s41586-019-1631-3'
        self.protocol = '10x'
        self.normalization = 'raw'
        self.var_symbol_col = 'index'
        self.obs_key_cellontology_original = 'annotation_lineage'
        self.obs_key_state_exact = 'condition'
        self.obs_key_healthy = self.obs_key_state_exact
        self.healthy_state_healthy = 'Uninjured'

        self.class_maps = {
            "0": {
                'MPs': 'MP',
                'Tcells': 'Tcells',
                'ILCs': 'ILC',
                'Endothelia': 'Endothelia',
                'Bcells': 'Bcells',
                'pDCs': 'pDCs',
                'Plasma Bcells': 'Plasma B cell',
                'Mast cells': 'Mast cell',
                'Mesenchyme': 'Mesenchyme',
                'Cholangiocytes': 'Cholangiocytes',
                'Hepatocytes': 'Hepatocytes',
                'Mesothelia': 'Mesothelia',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "liver", "ramachandran.h5ad")
        self.adata = anndata.read(fn)
