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
        self.id = "mouse_pancreas_2018_microwell-seq_han_001_10.1016/j.cell.2018.02.001"
        self.download = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "pancreas"
        self.sub_tissue = "pancreas"
        self.annotated = True

        self.class_maps = {
            "0": {
                'Acinar cell(Pancreas)': 'pancreatic acinar cell',
                'Dendrtic cell(Pancreas)': 'dendritic cell',
                'Ductal cell(Pancreas)': 'pancreatic ductal cell',
                'Endocrine cell(Pancreas)': "endocrine cell",
                'Dividing cell(Pancreas)': "endocrine cell",
                'Endothelial cell_Fabp4 high(Pancreas)': 'endothelial cell',
                'Endothelial cell_Lrg1 high(Pancreas)': 'endothelial cell',
                'Endothelial cell_Tm4sf1 high(Pancreas)': 'endothelial cell',
                'Erythroblast_Hbb-bt high(Pancreas)': 'erythroblast',
                'Erythroblast_Igkc high(Pancreas)': 'erythroblast',
                'Granulocyte(Pancreas)': 'granulocyte',
                'Macrophage_Ly6c2 high(Pancreas)': 'macrophage',
                'Macrophage(Pancreas)': 'macrophage',
                'Glial cell(Pancreas)': 'glial cell',
                'Smooth muscle cell_Acta2 high(Pancreas)': 'smooth muscle cell',
                'Smooth muscle cell_Rgs5 high(Pancreas)': 'smooth muscle cell',
                'Stromal cell_Fn1 high(Pancreas)': 'stromal cell',
                'Stromal cell_Mfap4 high(Pancreas)': 'stromal cell',
                'Stromal cell_Smoc2 high(Pancreas)': 'stromal cell',
                'T cell(Pancreas)': 't cell',
                'B cell(Pancreas)': 'b cell',
                'β-cell(Pancreas)': "pancreatic B cell"
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "Pancreas_dge.txt.gz")
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
        self.set_unkown_class_id(ids=[np.nan, "nan"])
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_types_original] = self.adata.obs["Annotation"].values.tolist()
        self.adata.obs[self._ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[self._ADATA_IDS_SFAIRA.state_exact] = "healthy"

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
