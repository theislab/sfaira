import anndata
import os
from typing import Union
from .external import DatasetBase
import pandas as pd


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data files which can be obtained from the `download_website` and
    `download_website_meta` attributes of this class.

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
        super().__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_pancreas_2016_smartseq2_segerstolpe_001_10.1016/j.cmet.2016.08.020"
        self.download = "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.processed.1.zip"
        self.download_meta = "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.sdrf.txt"
        self.organ = "pancreas"
        self.sub_tissue = "pancreas"
        self.author = "Sandberg"
        self.year = 2016
        self.doi = "10.1016/j.cmet.2016.08.020"
        self.protocol = 'Smartseq2'
        self.normalization = 'raw'
        self.var_symbol_col = 'index'
        self.obs_key_cellontology_original = 'Characteristics[cell type]'
        self.obs_key_state_exact = 'Characteristics[disease]'
        self.obs_key_healthy = self.obs_key_state_exact
        self.healthy_state_healthy = 'normal'

        self.class_maps = {
            "0": {
                'alpha cell': 'Alpha cell',
                'ductal cell': 'Ductal cell',
                'beta cell': 'Beta cell',
                'gamma cell': 'Gamma cell',
                'acinar cell': 'Acinar cell',
                'delta cell': 'Delta cell',
                'PSC cell': 'PSC cell',
                'unclassified endocrine cell': 'Unclassified endocrine cell',
                'co-expression cell': 'Co-expression cell',
                'endothelial cell': 'Endothelial cell',
                'epsilon cell': 'Epsilon cell',
                'mast cell': 'Mast cell',
                'MHC class II cell': 'MHC class II cell',
                'unclassified cell': 'Unknown',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = [
                os.path.join(self.path, "human", "pancreas", "E-MTAB-5061.processed.1.zip"),
                os.path.join(self.path, "human", "pancreas", "E-MTAB-5061.sdrf.txt")
            ]
        df = pd.read_csv(fn[0], sep='\t')
        df.index = df.index.get_level_values(0)
        df = df.drop('#samples', axis=1)
        df = df.T.iloc[:, :26178]
        self.adata = anndata.AnnData(df)
        self.adata.obs = pd.read_csv(fn[1], sep='\t').set_index('Source Name').loc[self.adata.obs.index]
        # filter observations which are not cells (empty wells, low quality cells etc.)
        self.adata = self.adata[self.adata.obs['Characteristics[cell type]'] != 'not applicable'].copy()
