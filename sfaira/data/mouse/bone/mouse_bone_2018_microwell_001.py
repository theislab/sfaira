import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA


class Dataset(DatasetBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "mouse"
        self.id = "mouse_bone_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download_website = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "bone"
        self.sub_tissue = "marrow"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'B cell_Igkc high(Bone-Marrow)': 'naive B cell',
                'Dendritic cell_H2-Eb1 high(Bone-Marrow)': 'dendritic cell',
                'Dendritic cell_Siglech high(Bone-Marrow)': 'dendritic cell',
                'Macrophage_Ms4a6c high(Bone-Marrow)': 'macrophage',
                'Macrophage_S100a4 high(Bone-Marrow)': 'macrophage',
                'Erythroblast(Bone-Marrow)': 'erythroid progenitor',
                'Mast cell(Bone-Marrow)': 'mast cell',
                'Monocyte_Mif high(Bone-Marrow)': 'monocyte',
                'Monocyte_Prtn3 high(Bone-Marrow)': 'monocyte',
                'Neutrophil progenitor(Bone-Marrow)': 'neutrophil progenitor',
                'Neutrophil_Cebpe high(Bone-Marrow)': 'neutrophil',
                'Neutrophil_Fcnb high(Bone-Marrow)': 'neutrophil',
                'Neutrophil_Mmp8 high(Bone-Marrow)': 'neutrophil',
                'Neutrophil_Ngp high(Bone-Marrow)': 'neutrophil',
                'Hematopoietic stem progenitor cell(Bone-Marrow)': 'hematopoietic precursor cell',
                'Pre-pro B cell(Bone-Marrow)': 'early pro-B cell',
                'T cell_Ms4a4b high(Bone-Marrow)': 'CD4-positive, alpha-beta T cell'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "BoneMarrow1_dge.txt.gz")
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
        self.adata.obs[ADATA_IDS_SFAIRA.cell_types_original] = self.adata.obs['Annotation']
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)

