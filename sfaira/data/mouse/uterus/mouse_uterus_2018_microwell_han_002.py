import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetBase


class Dataset(DatasetBase):

    id: str

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "mouse"
        self.id = "mouse_uterus_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "uterus"
        self.sub_tissue = "uterus"
        self.annotated = True

        self.class_maps = {
            "0": {
                'B cell(Uterus)': 'B cell',
                'Dendritic cell(Uterus)': 'dendritic cell',
                'Endothelial cell_Cldn5 high(Uterus)': 'endothelial cell',
                'Endothelial cell_Tm4sf1 high(Uterus)': 'endothelial cell',
                'Glandular epithelium_Ltf high(Uterus)': 'glandular epithelial cell',
                'Glandular epithelium_Sprr2f high(Uterus)': 'glandular epithelial cell',
                'Granulocyte(Uterus)': 'granulocyte',
                'Keratinocyte(Uterus)': 'keratinocyte',
                'Macrophage(Uterus)': 'macrophage',
                'Monocyte(Uterus)': 'monocyte',
                'Muscle cell_Mgp high(Uterus)': 'muscle cell',
                'Muscle cell_Pcp4 high(Uterus)': 'muscle cell',
                'Smooth muscle cell_Rgs5 high(Uterus)': 'smooth muscle cell',
                'NK cell(Uterus)': 'NK cell',
                'Stromal cell_Ccl11 high(Uterus)': 'stromal cell',
                'Stromal cell_Cxcl14 high(Uterus)': 'stromal cell',
                'Stromal cell_Gm23935 high(Uterus)': 'stromal cell',
                'Stromal cell_Has1 high(Uterus)': 'stromal cell',
                'Stromal cell_Hsd11b2 high(Uterus)': 'stromal cell',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Uterus2_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        celltypes = pandas.read_csv(fn_meta, index_col=1)
        celltypes = celltypes.drop(['Unnamed: 0'], axis=1)

        data = pandas.read_csv(fn, sep=' ', header=0)
        self.adata = anndata.AnnData(data.T)
        self.adata = self.adata[np.array([x in celltypes.index for x in self.adata.obs_names])].copy()
        self.adata.obs = celltypes.loc[self.adata.obs_names, :]

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = "Guo"
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = "2018"
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = "10.1016/j.cell.2018.02.001"
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = "microwell-seq"
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue  # TODO
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = "mouse"
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'raw'
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs["Annotation"].values.tolist()
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_types_original] = self.adata.obs["Annotation"].values.tolist()
        self.set_unkown_class_id(ids=[np.nan, "nan"])
        self.adata.obs[self._ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[self._ADATA_IDS_SFAIRA.state_exact] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
