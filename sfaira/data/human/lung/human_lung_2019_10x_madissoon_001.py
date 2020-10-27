import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS
import anndata


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data file which can be obtained from the `download_website` attribute of
    this class.

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_lung_2019_10x_madissoon_001._10.1186/s13059-019-1906-x"
        self.download_website = "https://covid19.cog.sanger.ac.uk/madissoon19_lung.processed.h5ad"
        self.organ = "lung"
        self.sub_tissue = "parenchyma"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'T_CD4': 'T cell lineage',
                'Mast_cells': 'Mast cells',
                'Monocyte': 'Monocytes',
                'Blood_vessel': '2_Blood vessels',
                'Ciliated': 'Multiciliated lineage',
                'Macrophage_MARCOneg': 'Macrophages',
                'DC_plasmacytoid': 'Dendritic cells',
                'DC_1': 'Dendritic cells',
                'Muscle_cells': '2_Smooth Muscle',
                'Macrophage_MARCOpos': 'Macrophages',
                'T_cells_Dividing': 'T cell lineage',
                'DC_Monocyte_Dividing': 'Dendritic cells',
                'B_cells': 'B cell lineage',
                'T_CD8_CytT': 'T cell lineage',
                'NK_Dividing': 'Innate lymphoid cells',
                'T_regulatory': 'T cell lineage',
                'DC_2': 'Dendritic cells',
                'Alveolar_Type2': 'AT2',
                'Plasma_cells': 'B cell lineage',
                'NK': 'Innate lymphoid cells',
                'Alveolar_Type1': 'AT1',
                'Fibroblast': '2_Fibroblast lineage',
                'DC_activated': 'Dendritic cells',
                'Macrophage_Dividing': 'Macrophages',
                'Lymph_vessel': 'Lymphatic EC',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human/lung/madissoon19_lung.processed.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[ADATA_IDS.lab] = 'Meyer'
        self.adata.uns[ADATA_IDS.year] = 2020
        self.adata.uns[ADATA_IDS.doi] = "10.1186/s13059-019-1906-x"
        self.adata.uns[ADATA_IDS.protocol] = '10x'
        self.adata.uns[ADATA_IDS.organ] = self.organ
        self.adata.uns[ADATA_IDS.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS.animal] = "human"
        self.adata.uns[ADATA_IDS.id] = self.id
        self.adata.uns[ADATA_IDS.wget_download] = self.download_website
        self.adata.uns[ADATA_IDS.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS.normalization] = 'raw'

        self.adata.obs[ADATA_IDS.cell_ontology_class] = self.adata.obs['CellType']
        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
        self.adata.obs[ADATA_IDS.healthy] = True
        self.adata.obs['state_exact'] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col='gene.ids.HCATisStab7509734', new_index=ADATA_IDS.gene_id_ensembl)
