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
        self.id = "mouse_rib_2018_microwell-seq_han_003_10.1016/j.cell.2018.02.001"
        self.download_website = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "rib"
        self.sub_tissue = "rib"
        self.annotated = True

        self.class_maps = {
            "0": {
                'B cell(Neonatal-Rib)': 'B cell',
                'Cartilage cell_Clu high(Neonatal-Rib)': 'cartilage cell',
                'Cartilage cell_Col2a1 high(Neonatal-Rib)': 'cartilage cell',
                'Cartilage cell_Cxcl14 high(Neonatal-Rib)': 'cartilage cell',
                'Cartilage cell_Ppa1 high(Neonatal-Rib)': 'cartilage cell',
                'Cartilage cell_Prg4 high(Neonatal-Rib)': 'cartilage cell',
                'Dividing cell(Neonatal-Rib)': 'proliferative cell',
                'Endothelial cell(Neonatal-Rib)': 'endothelial cell',
                'Erythroblast_Hba-a1 high(Neonatal-Rib)': 'erythroblast',
                'Erythroblast_Ttr high(Neonatal-Rib)': 'erythroblast',
                'Granulocyte(Neonatal-Rib)': 'granulocyte',
                'Macrophage_C1qc high(Neonatal-Rib)': 'macrophage',
                'Macrophage_Ctss high(Neonatal-Rib)': 'macrophage',
                'Muscle cell(Neonatal-Rib)': 'muscle cell',
                'Muscle cell_Acta2 high(Neonatal-Rib)': 'muscle cell',
                'Muscle cell_Actc1 high(Neonatal-Rib)': 'muscle cell',
                'Neuron_Mpz high(Neonatal-Rib)': 'neuron',
                'Neuron_Stmn2 high(Neonatal-Rib)': 'neuron',
                'Neutrophil(Neonatal-Rib)': 'neutrophil',
                'Neutrophil_Elane high(Neonatal-Rib)': 'neutrophil',
                'Oligodendrocyte(Neonatal-Rib)': 'oligodendrocyte',
                'Osteoblast(Neonatal-Rib)': 'osteoblast',
                'Osteoclast(Neonatal-Rib)': 'osteoclast',
                'Stromal cell_Acta1 high(Neonatal-Rib)': 'stromal cell',
                'Stromal cell_Tnmd high(Neonatal-Rib)': 'stromal cell',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas/500more_dge", "NeonatalRib3_dge.txt.gz")
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
