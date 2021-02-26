import anndata
import os
from typing import Union

from sfaira.data import DatasetBase
from sfaira.consts import ADATA_IDS_CELLXGENE


class Dataset(DatasetBase):
    """
    This is a dataloader for downloaded h5ad from cellxgene.

    :param data_path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            data_path: Union[str, None],
            fn: str,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(data_path=data_path, meta_path=meta_path, **kwargs)
        self.fn = fn

        self.obs_key_cellontology_class = ADATA_IDS_CELLXGENE.cell_ontology_class
        self.obs_key_cellontology_id = ADATA_IDS_CELLXGENE.cell_ontology_id
        self.obs_key_cellontology_original = ADATA_IDS_CELLXGENE.cell_types_original
        self.obs_key_dev_stage = ADATA_IDS_CELLXGENE.dev_stage
        self.obs_key_ethnicity = ADATA_IDS_CELLXGENE.ethnicity
        self.obs_key_healthy = ADATA_IDS_CELLXGENE.healthy
        self.obs_key_sex = ADATA_IDS_CELLXGENE.sex
        self.obs_key_organism = ADATA_IDS_CELLXGENE.organism
        self.obs_key_state_exact = ADATA_IDS_CELLXGENE.state_exact

        self.healthy_state_healthy = ADATA_IDS_CELLXGENE.disease_state_healthy

        self.var_ensembl_col = ADATA_IDS_CELLXGENE.gene_id_ensembl
        self.var_symbol_col = ADATA_IDS_CELLXGENE.gene_id_names

        self.class_maps = {
            "0": {},
        }

    def _load(self):
        """
        Note that in contrast to data set specific data loaders, here, the core attributes are only identified from
        the data in this function and are not already set in the constructor. These attributes can still be
        used through meta data containers after the data was loaded once.

        :return:
        """
        fn = os.path.join(self.data_dir_base, self.fn)
        adata = anndata.read(fn)
        adata.X = adata.raw.X
        # TODO delete raw?

        self.author = adata.uns[ADATA_IDS_CELLXGENE.author][ADATA_IDS_CELLXGENE.author_names]
        self.doi = adata.uns[ADATA_IDS_CELLXGENE.doi]
        self.download_url_data = self.download_url_data
        self.id = self.id
        self.normalization = 'raw'
        self.organ = str(self.fn).split("_")[3]  # TODO interface this properly
        # self.organ = adata.obs["tissue"].values[0]
        self.organism = adata.obs[ADATA_IDS_CELLXGENE.organism].values[0]
        self.protocol = adata.obs[ADATA_IDS_CELLXGENE.protocol].values[0]
        self.year = adata.uns[ADATA_IDS_CELLXGENE.year]
