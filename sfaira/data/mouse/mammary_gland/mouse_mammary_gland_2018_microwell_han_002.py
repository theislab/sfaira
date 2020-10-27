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
        self.id = "mouse_mammary_gland_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.download_website = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "mammary_gland"
        self.sub_tissue = "mammary_gland"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'B cell_Cd79a&Fcer2a high(Mammary-Gland-Virgin)': 'B cell',
                'B cell_Cd79a&Iglc2 high(Mammary-Gland-Virgin)': 'B cell',
                'B cell_Jchain high(Mammary-Gland-Virgin)': 'B cell',
                'Dendritic cell_Cst3 high(Mammary-Gland-Virgin)': 'dendritic cell',
                'Dendritic cell_Fscn1 high(Mammary-Gland-Virgin)': 'dendritic cell',
                'Dendritic cell_Siglech high(Mammary-Gland-Virgin)': 'dendritic cell',
                'Dividing cell(Mammary-Gland-Virgin)': 'proliferative cell',
                'Luminal cell_Krt19 high (Mammary-Gland-Virgin)': 'luminal epithelial cell of mammary gland',
                'Luminal progenitor(Mammary-Gland-Virgin)': 'luminal progenitor cell',
                'Macrophage_C1qc high(Mammary-Gland-Virgin)': 'macrophage',
                'Macrophage_Lyz1 high(Mammary-Gland-Virgin)': 'macrophage',
                'NK cell(Mammary-Gland-Virgin)': 'NK cell',
                'Stem and progenitor cell(Mammary-Gland-Virgin)': 'stem and progenitor cell',
                'Stromal cell_Col3a1 high(Mammary-Gland-Virgin)': 'stromal cell',
                'Stromal cell_Pi16 high(Mammary-Gland-Virgin)': 'stromal cell',
                'T cell_Cd8b1 high(Mammary-Gland-Virgin)': 'T cell',
                'T cell_Ly6c2 high(Mammary-Gland-Virgin)': 'T cell',
                'T-cells_Ctla4 high(Mammary-Gland-Virgin)': 'T cell'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse/temp_mouse_atlas/500more_dge/MammaryGland.Virgin2_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse/temp_mouse_atlas/MCA_CellAssignments.csv")

        celltypes = pandas.read_csv(fn_meta, index_col=1)
        celltypes = celltypes.drop(['Unnamed: 0'], axis=1)

        data = pandas.read_csv(fn, sep=' ', header=0)
        self.adata = anndata.AnnData(data.T)
        self.adata = self.adata[np.array([x in celltypes.index for x in self.adata.obs_names])].copy()
        self.adata.obs = celltypes.loc[self.adata.obs_names, :]

        self.adata.uns[ADATA_IDS.author] = "Guo"
        self.adata.uns[ADATA_IDS.year] = "2018"
        self.adata.uns[ADATA_IDS.doi] = "10.1016/j.cell.2018.02.001"
        self.adata.uns[ADATA_IDS.protocol] = "microwell-seq"
        self.adata.uns[ADATA_IDS.organ] = self.organ
        self.adata.uns[ADATA_IDS.subtissue] = self.sub_tissue  # TODO
        self.adata.uns[ADATA_IDS.animal] = "mouse"
        self.adata.uns[ADATA_IDS.id] = self.id
        self.adata.uns[ADATA_IDS.wget_download] = self.download_website
        self.adata.uns[ADATA_IDS.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS.normalization] = 'raw'
        self.adata.obs[ADATA_IDS.cell_ontology_class] = self.adata.obs["Annotation"].values.tolist()
        self.set_unkown_class_id(ids=[np.nan, "nan"])
        self.adata.obs[ADATA_IDS.cell_types_original] = self.adata.obs["Annotation"].values.tolist()
        self.adata.obs[ADATA_IDS.healthy] = True
        self.adata.obs[ADATA_IDS.state_exact] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS.gene_id_ensembl)

