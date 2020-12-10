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
        self.id = "mouse_lung_2018_microwell-seq_han_003_10.1016/j.cell.2018.02.001"
        self.download_website = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "lung"
        self.sub_tissue = "lung"
        self.annotated = True

        self.class_maps = {
            "0": {
                'AT1 Cell(Lung)': 'alveolar epithelial cell type I',
                'AT2 Cell(Lung)': 'alveolar epithelial cell type II',
                'Alveolar bipotent progenitor(Lung)': 'alveolar bipotent progenitor',
                'Alveolar macrophage_Ear2 high(Lung)': 'alveolar macrophage',
                'Alveolar macrophage_Pclaf high(Lung)': 'alveolar macrophage',
                'B Cell(Lung)': 'B cell',
                'Basophil(Lung)': 'basophil',
                'Ciliated cell(Lung)': 'ciliated cell',
                'Clara Cell(Lung)': 'clara cell',
                'Conventional dendritic cell_Gngt2 high(Lung)': "dendritic cell",
                'Conventional dendritic cell_H2-M2 high(Lung)': "dendritic cell",
                'Conventional dendritic cell_Mgl2 high(Lung)': "dendritic cell",
                'Conventional dendritic cell_Tubb5 high(Lung)': "dendritic cell",
                'Dendritic cell_Naaa high(Lung)': "dendritic cell",
                'Dividing T cells(Lung)': "T cell",
                'Dividing cells(Lung)': 'unknown',
                'Dividing dendritic cells(Lung)': "dendritic cell",
                'Endothelial cell_Kdr high(Lung)': "endothelial cell",
                'Endothelial cell_Tmem100 high(Lung)': "endothelial cell",
                'Endothelial cells_Vwf high(Lung)': "endothelial cell",
                'Eosinophil granulocyte(Lung)': 'eosinophil',
                'Ig−producing B cell(Lung)': 'B cell',
                'Interstitial macrophage(Lung)': 'lung macrophage',
                'Monocyte progenitor cell(Lung)': 'monocyte progenitor',
                'NK Cell(Lung)': 'NK cell',
                'Neutrophil granulocyte(Lung)': 'neutrophil',
                'Nuocyte(Lung)': 'nuocyte',
                'Plasmacytoid dendritic cell(Lung)': "plasmacytoid dendritic cell",
                'Stromal cell_Acta2 high(Lung)': 'stromal cell',
                'Stromal cell_Dcn high(Lung)': 'stromal cell',
                'Stromal cell_Inmt high(Lung)': 'stromal cell',
                'T Cell_Cd8b1 high(Lung)': "CD8-positive, alpha-beta T cell",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Lung3_dge.txt.gz")
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
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'
        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs["Annotation"].values.tolist()
        self.set_unkown_class_id(ids=[np.nan, "nan"])
        self.adata.obs[ADATA_IDS_SFAIRA.cell_types_original] = self.adata.obs["Annotation"].values.tolist()
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)

