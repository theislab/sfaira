import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS
import anndata
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
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_pancreas_2016_smartseq2_segerstolpe_001_10.1016/j.cmet.2016.08.020"
        self.download_website = "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.processed.1.zip"
        self.download_website_meta = "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5061/E-MTAB-5061.sdrf.txt"
        self.organ = "pancreas"
        self.sub_tissue = "pancreas"
        self.has_celltypes = True

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
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = [
                    os.path.join(self.path, "human/pancreas/E-MTAB-5061.processed.1.zip"),
                    os.path.join(self.path, "human/pancreas/E-MTAB-5061.sdrf.txt")
                ]
            df = pd.read_csv(fn[0], sep='\t')
            df.index = df.index.get_level_values(0)
            df = df.drop('#samples', axis=1)
            df = df.T.iloc[:, :26178]
            self.adata = anndata.AnnData(df)
            self.adata.obs = pd.read_csv(fn[1], sep='\t').set_index('Source Name').loc[self.adata.obs.index]
            # filter observations which are not cells (empty wells, low quality cells etc.)
            self.adata = self.adata[self.adata.obs['Characteristics[cell type]'] != 'not applicable'].copy()

        self.adata.uns[ADATA_IDS.lab] = "Sandberg"
        self.adata.uns[ADATA_IDS.year] = 2016
        self.adata.uns[ADATA_IDS.doi] = "10.1016/j.cmet.2016.08.020"
        self.adata.uns[ADATA_IDS.protocol] = 'Smartseq2'
        self.adata.uns[ADATA_IDS.organ] = self.organ
        self.adata.uns[ADATA_IDS.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS.animal] = "human"
        self.adata.uns[ADATA_IDS.id] = self.id
        self.adata.uns[ADATA_IDS.wget_download] = self.download_website
        self.adata.uns[ADATA_IDS.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS.normalization] = 'raw'

        self.adata.obs[ADATA_IDS.healthy] = [True if line == 'normal' else False for line in self.adata.obs['Characteristics[disease]']]
        self.adata.obs[ADATA_IDS.state_exact] = self.adata.obs['Characteristics[disease]'].astype('category')
        self.adata.obs[ADATA_IDS.state_exact] = self.adata.obs[ADATA_IDS.state_exact].cat.rename_categories({'normal':'healthy', 'type II diabetes mellitus':'type II diabetes mellitus'})

        self.adata.obs[ADATA_IDS.cell_ontology_class] = self.adata.obs['Characteristics[cell type]']
        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS.gene_id_ensembl)
