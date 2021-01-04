import os
from typing import Union
import pandas as pd
import anndata

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data file which can be obtained from the `download_website` and
    `download_website_meta` attributes of this class.

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
        self.organism = "loaders"
        self.id = "human_placenta_2018_10x_ventotormo_10.1038/s41586-018-0698-6"
        self.download = 'https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6701/E-MTAB-6701.processed.1.zip'
        self.download_meta = 'https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6701/E-MTAB-6701.processed.2.zip'
        self.organ = "placenta"
        self.sub_tissue = "placenta, decidua, blood"
        self.author = 'Teichmann'
        self.year = 2018
        self.doi = '10.1038/s41586-018-0698-6'
        self.protocol = "10x"
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = "healthy"
        self.var_symbol_col = 'names'
        self.var_ensembl_col = 'ensembl'
        self.obs_key_cellontology_original = 'annotation'
        self._obs_key_subtissue = 'location'

        self.class_maps = {
            "0": {
                'DC1': 'Dendritic Cells 1',
                'DC2': 'Dendritic Cells 2',
                'EVT': 'Extravillous Trophoblasts',
                'Endo (f)': 'Endothelial Cells f',
                'Endo (m)': 'Endothelial Cells m',
                'Endo L': 'Endothelial Cells L',
                'Epi1': 'Epithelial Glandular Cells 1',
                'Epi2': 'Epithelial Glandular Cells 2',
                'Granulocytes': 'Granulocytes',
                'HB': 'Hofbauer Cells',
                'ILC3': 'ILC3',
                'MO': 'Monocyte',
                'NK CD16+': 'NK Cells CD16+',
                'NK CD16-': 'NK Cells CD16-',
                'Plasma': 'B cell (Plasmocyte)',
                'SCT': 'Syncytiotrophoblasts',
                'Tcells': 'T cell',
                'VCT': 'Villous Cytotrophoblasts',
                'dM1': 'Decidual Macrophages 1',
                'dM2': 'Decidual Macrophages 2',
                'dM3': 'Decidual Macrophages 3',
                'dNK p': 'Decidual NK Cells p',
                'dNK1': 'Decidual NK Cells 1',
                'dNK2': 'Decidual NK Cells 2',
                'dNK3': 'Decidual NK Cells 3',
                'dP1': 'Perivascular Cells 1',
                'dP2': 'Perivascular Cells 2',
                'dS1': 'Decidual Stromal Cells 1',
                'dS2': 'Decidual Stromal Cells 2',
                'dS3': 'Decidual Stromal Cells 3',
                'fFB1': 'Fibroblasts 1',
                'fFB2': 'Fibroblasts 2',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = [
                os.path.join(self.path, "loaders", "placenta", "E-MTAB-6701.processed.1.zip"),
                os.path.join(self.path, "loaders", "placenta", "E-MTAB-6701.processed.2.zip"),
            ]
        self.adata = anndata.AnnData(pd.read_csv(fn[0], sep='\t', index_col='Gene').T)
        df = pd.read_csv(fn[1], sep='\t')
        for i in df.columns:
            self.adata.obs[i] = [df.loc[j][i] for j in self.adata.obs.index]

        self.adata.var['ensembl'] = [i.split("_")[1] for i in self.adata.var.index]
        self.adata.var['names'] = [i.split("_")[0] for i in self.adata.var.index]
        self.adata.var = self.adata.var.reset_index().reset_index().drop('index', axis=1)
        self.adata = self.adata[:, ~self.adata.var.index.isin(
            ['', '-1', '-10', '-11', '-2', '-3', '-4', '-5', '-6', '-7', '-8', '-9', 'A.2', 'A.3'])].copy()
