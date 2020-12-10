import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA


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
        self.id = "mouse_blood_2018_microwell-seq_han_005_10.1016/j.cell.2018.02.001"
        self.download_website = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "blood"
        self.sub_tissue = "blood"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'B cell_Igha high(Peripheral_Blood)': 'B cell',
                'B cell_Ly6d high(Peripheral_Blood)': 'B cell',
                'B cell_Rps27rt high(Peripheral_Blood)': 'B cell',
                'B cell_Vpreb3 high(Peripheral_Blood)': 'B cell',
                'Basophil_Prss34 high(Peripheral_Blood)': 'basophil',
                'Dendritic cell_Siglech high(Peripheral_Blood)': 'dendritic cell',
                'Erythroblast_Car2 high(Peripheral_Blood)': 'erythroblast',
                'Erythroblast_Hba-a2 high(Peripheral_Blood)': 'erythroblast',
                'Macrophage_Ace high(Peripheral_Blood)': 'macrophage',
                'Macrophage_Flt-ps1 high(Peripheral_Blood)': 'macrophage',
                'Macrophage_Pf4 high(Peripheral_Blood)': 'macrophage',
                'Macrophage_S100a4 high(Peripheral_Blood)': 'macrophage',
                'Monocyte_Elane high(Peripheral_Blood)': 'monocyte',
                'Monocyte_F13a1 high(Peripheral_Blood)': 'monocyte',
                'NK cell_Gzma high(Peripheral_Blood)': 'NK cell',
                'Neutrophil_Camp high(Peripheral_Blood)': 'neutrophil',
                'Neutrophil_Il1b high(Peripheral_Blood)': 'neutrophil',
                'Neutrophil_Ltf high(Peripheral_Blood)': 'neutrophil',
                'Neutrophil_Retnlg high(Peripheral_Blood)': 'neutrophil',
                'T cell_Gm14303 high(Peripheral_Blood)': 'T cell',
                'T cell_Trbc2 high(Peripheral_Blood)': 'T cell'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "PeripheralBlood5_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse", "temp_mouse_atlas", "MCA_CellAssignments.csv")

        celltypes = pandas.read_csv(fn_meta, index_col=1)
        celltypes = celltypes.drop(['Unnamed: 0'], axis=1)

        data = pandas.read_csv(fn, sep=' ', header=0)
        self.adata = anndata.AnnData(data.T)
        self.adata = self.adata[np.array([x in celltypes.index for x in self.adata.obs_names])].copy()
        self.adata.obs = celltypes.loc[self.adata.obs_names, :]

        self.adata.uns[ADATA_IDS_SFAIRA.author] = "Guo"
        self.adata.uns[ADATA_IDS_SFAIRA.year] = "2018"
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = "10.1016/j.cell.2018.02.001"
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = "microwell-seq"
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue  # TODO
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "mouse"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.has_celltypes
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'
        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs["Annotation"].values.tolist()
        self.set_unkown_class_id(ids=[np.nan, "nan"])
        self.adata.obs[ADATA_IDS_SFAIRA.cell_types_original] = self.adata.obs["Annotation"].values.tolist()
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
