import anndata
import os
from typing import Union
from .external import DatasetBase


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
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_liver_2019_10x_ramachandran_001_10.1038/s41586-019-1631-3"
        self.download = "https://datashare.is.ed.ac.uk/bitstream/handle/10283/3433/tissue.rdata"
        self.download_meta = None
        self.organ = "liver"
        self.sub_tissue = "liver"
        self.annotated = True

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
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "liver", "ramachandran.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = 'Henderson'
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = '10.1038/s41586-019-1631-3'
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs["annotation_lineage"]
        self.adata.obs[self._ADATA_IDS_SFAIRA.healthy] = [i == 'Uninjured' for i in self.adata.obs["condition"]]
        self.adata.obs[self._ADATA_IDS_SFAIRA.state_exact] = ['healthy' if i == 'Uninjured' else i for i in self.adata.obs["condition"]]

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
