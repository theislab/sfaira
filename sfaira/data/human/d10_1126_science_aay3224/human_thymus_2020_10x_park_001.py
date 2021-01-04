import anndata
import os
from typing import Union
import numpy as np

from sfaira.data import DatasetBase


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
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.organism = "human"
        self.id = "human_thymus_2020_10x_park_001_10.1126/science.aay3224"
        self.download = "https://covid19.cog.sanger.ac.uk/park20.processed.h5ad"
        self.download_meta = None
        self.organ = "thymus"
        self.sub_tissue = "fetal thymus"
        self.author = "Teichmann"
        self.year = 2020
        self.doi = "10.1126/science.aay3224"
        self.protocol = '10x'
        self.normalization = 'norm'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'
        self.obs_key_cellontology_original = 'Anno_level_fig1'

        self.class_maps = {
            "0": {
                'B_memory': 'B_memory',
                'B_naive': 'B_naive',
                'B_plasma': 'B_plasma',
                'B_pro/pre': 'B_pro/pre',
                'CD4+T': 'CD4+T',
                'CD4+Tmem': 'CD4+Tmem',
                'CD8+T': 'CD8+T',
                'CD8+Tmem': 'CD8+Tmem',
                'CD8αα': 'CD8αα',
                'DC1': 'DC1',
                'DC2': 'DC2',
                'DN': 'DN',
                'DP': 'DP',
                'ETP': 'ETP',
                'Endo': 'Endo',
                'Epi_GCM2': 'Epi_GCM2',
                'Ery': 'Ery',
                'Fb_1': 'Fb_1',
                'Fb_2': 'Fb_2',
                'Fb_cycling': 'Fb_cycling',
                'ILC3': 'ILC3',
                'Lymph': 'Lymph',
                'Mac': 'Mac',
                'Mast': 'Mast',
                'Mgk': 'Mgk',
                'Mono': 'Mono',
                'NK': 'NK',
                'NKT': 'NKT',
                'NMP': 'NMP',
                'T(agonist)': 'T(agonist)',
                'TEC(myo)': 'TEC(myo)',
                'TEC(neuro)': 'TEC(neuro)',
                'Treg': 'Treg',
                'VSMC': 'VSMC',
                'aDC': 'aDC',
                'cTEC': 'cTEC',
                'mTEC(I)': 'mTEC(I)',
                'mTEC(II)': 'mTEC(II)',
                'mTEC(III)': 'mTEC(III)',
                'mTEC(IV)': 'mTEC(IV)',
                'mcTEC': 'mcTEC',
                'pDC': 'pDC',
                'αβT(entry)': 'alpha_beta_T(entry)',
                'γδT': 'gamma_delta_T',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = os.path.join(self.path, "human", "thymus", "park20.processed.h5ad")
        self.adata = anndata.read(fn)
        self.adata.X = np.expm1(self.adata.X)
