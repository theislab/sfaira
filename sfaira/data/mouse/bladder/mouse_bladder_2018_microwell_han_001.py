import anndata
import numpy as np
import os
import pandas
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS


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
        self.id = "mouse_bladder_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download_website = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "bladder"
        self.sub_tissue = "bladder"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                "Endothelial cell_Ly6c1 high(Bladder)": 'endothelial cell',
                "Vascular endothelial cell(Bladder)": 'endothelial cell',
                'Urothelium(Bladder)': 'bladder urothelial cell',
                'Dendritic cell_Cd74 high(Bladder)': 'dendritic cell',
                'Dendritic cell_Lyz2 high(Bladder)': 'dendritic cell',
                'Macrophage_Pf4 high(Bladder)': 'macrophage',
                'NK cell(Bladder)': 'NK cell',
                'Basal epithelial cell(Bladder)': 'basal epithelial cell',
                'Epithelial cell_Upk3a high(Bladder)': 'epithelial cell',
                'Epithelial cell_Gm23935 high(Bladder)': 'epithelial cell',
                'Mesenchymal stromal cell(Bladder)': 'mesenchymal stromal cell',
                'Stromal cell_Dpt high(Bladder)': 'stromal cell',
                'Stromal cell_Car3 high(Bladder)': 'stromal cell',
                'Smooth muscle cell(Bladder)': 'smooth muscle cell',
                'Vascular smooth muscle progenitor cell(Bladder)': 'smooth muscle cell',
                'Umbrella cell(Bladder)': 'umbrella cell'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse/temp_mouse_atlas/500more_dge/Bladder_dge.txt.gz")
            fn_meta = os.path.join(self.path, "mouse/temp_mouse_atlas/MCA_CellAssignments.csv")

        celltypes = pandas.read_csv(fn_meta, index_col=1)
        celltypes = celltypes.drop(['Unnamed: 0'], axis=1)

        data = pandas.read_csv(fn, sep=' ', header=0)
        self.adata = anndata.AnnData(data.T)
        self.adata = self.adata[np.array([x in celltypes.index for x in self.adata.obs_names])].copy()
        self.adata.obs = celltypes.loc[self.adata.obs_names, :]

        self.adata.uns[ADATA_IDS.lab] = "Guo"
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
