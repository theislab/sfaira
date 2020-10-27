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
        self.id = "mouse_testis_2018_microwell-seq_han_002_10.1016/j.cell.2018.02.001"
        self.download_website = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.organ = "testis"
        self.sub_tissue = "testis"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'Elongating spermatid(Testis)': 'elongating spermatid',
                'Erythroblast_Hbb-bs high(Testis)': 'erythroblast',
                'Leydig cell(Testis)': 'leydig cell',
                'Macrophage_Lyz2 high(Testis)': 'macrophage',
                'Pre-Sertoli cell_Cst9 high(Testis)': 'pre-sertoli cell',
                'Pre-Sertoli cell_Ctsl high(Testis)': 'pre-sertoli cell',
                'Preleptotene spermatogonia(Testis)': 'preleptotene spermatogonia',
                'Sertoli cell(Testis)': 'sertoli cell',
                'Spermatids_1700016P04Rik high(Testis)': 'spermatid',
                'Spermatids_Cst13 high(Testis)': 'spermatid',
                'Spermatids_Hmgb4 high(Testis)': 'spermatid',
                'Spermatids_Tnp1 high(Testis)': 'spermatid',
                'Spermatocyte_1700001F09Rik high(Testis)': 'spermatocyte',
                'Spermatocyte_Cabs1 high(Testis)': 'spermatocyte',
                'Spermatocyte_Calm2 high(Testis)': 'spermatocyte',
                'Spermatocyte_Mesp1 high(Testis)': 'spermatocyte',
                'Spermatocyte_Slc2a3 high(Testis)': 'spermatocyte',
                'Spermatogonia_1700001P01Rik high(Testis)': 'spermatogonia',
                'Spermatogonia_Tbc1d23 high(Testis)': 'spermatogonia'
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse/temp_mouse_atlas/500more_dge/Testis2_dge.txt.gz")
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
