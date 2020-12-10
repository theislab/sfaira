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
        self.id = "mouse_placenta_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download_website = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "placenta"
        self.sub_tissue = "placenta"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'B cell(Placenta)': 'B cell',
                'Basophil(Placenta)': 'basophil',
                'Decidual stromal cell(Placenta)': 'decidual stromal cell',
                'Dendritic cell(Placenta)': 'dendritic cell',
                'Endodermal cell_Afp high(Placenta)': 'endodermal cell',
                'Endothelial cell_Maged2 high(Placenta)': 'endothelial cell',
                'Erythroblast_Hbb-y high(Placenta)': 'erythroblast',
                'Granulocyte monocyte progenitors(Placenta)': 'monocyte progenitor',
                'Granulocyte_Neat1 high(Placenta)': 'granulocyte',
                'Granulocyte_S100a9 high(Placenta)': 'granulocyte',
                'HSPC_Lmo2 high(Placenta)': 'HSPC',
                'Invasive spongiotrophoblast(Placenta)': 'invasive spongiotrophoblast',
                'Labyrinthine trophoblast(Placenta)': 'labyrinthine trophoblast',
                'Macrophage_Apoe high(Placenta)': 'macrophage',
                'Macrophage_Spp1 high(Placenta)': 'macrophage',
                'Megakaryocyte progenitor cell(Placenta)': 'megakaryocte',
                'Monocyte(Placenta)': 'monocyte',
                'NK cell(Placenta)': 'NK cell',
                'NKT cell(Placenta)': 'NKT cell',
                'PE lineage cell_Gkn2 high(Placenta)': 'PE lineage cell',
                'PE lineage cell_S100g high(Placenta)': 'PE lineage cell',
                'Progenitor trophoblast_Gjb3 high(Placenta)': 'trophoblast progenitor',
                'Spiral artery trophoblast giant cells(Placenta)': 'spiral artery trophoblast giant cells',
                'Spongiotrophoblast_Hsd11b2 high(Placenta)': 'spongiotrophoblast',
                'Spongiotrophoblast_Phlda2 high(Placenta)': 'spongiotrophoblast',
                'Stromal cell(Placenta)': 'stromal cell',
                'Stromal cell_Acta2 high(Placenta)': 'stromal cell',
                'Trophoblast progenitor_Taf7l high(Placenta)': 'trophoblast progenitor',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "PlacentaE14.1_dge.txt.gz")
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

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS_SFAIRA.gene_id_ensembl)
