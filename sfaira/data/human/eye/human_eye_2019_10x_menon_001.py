import anndata
import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data file which can be obtained from the `download_website` attribute of
    this class.

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
        self.id = "human_eye_2019_10x_menon_001_10.1038/s41467-019-12780-8"
        self.download_website = "https://covid19.cog.sanger.ac.uk/menon19.processed.h5ad"
        self.download_website_meta = None
        self.organ = "eye"
        self.sub_tissue = "retina"
        self.annotated = True

        self.class_maps = {
            "0": {
                'ACs': 'Amacrine cell',
                'BPs': 'BPs',
                'Cones': 'Retinal cone cell',
                'Endo': 'Endothelial cell',
                'HCs': 'Horizontal cells',
                'Macroglia': 'Macroglia',
                'Microglia': 'Microglia',
                'RGCs': 'Retinal ganglion cell',
                'Rods': 'Rods',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "eye", "menon19.processed.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[ADATA_IDS_SFAIRA.author] = 'Hafler'
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = '10.1038/s41467-019-12780-8'
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['CellType']
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
