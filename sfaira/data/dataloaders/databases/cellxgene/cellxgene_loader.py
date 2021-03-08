import anndata
import os
from typing import Union

from sfaira.data import DatasetBase
from sfaira.consts import AdataIdsCellxgene


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
        self._self._adata_ids_cellxgene = AdataIdsCellxgene()
        self.fn = fn

        self.obs_key_cellontology_class = self._adata_ids_cellxgene.cell_ontology_class
        self.obs_key_cellontology_id = self._adata_ids_cellxgene.cell_ontology_id
        self.obs_key_cellontology_original = self._adata_ids_cellxgene.cell_types_original
        self.obs_key_dev_stage = self._adata_ids_cellxgene.development_stage
        self.obs_key_ethnicity = self._adata_ids_cellxgene.ethnicity
        self.obs_key_healthy = self._adata_ids_cellxgene.healthy
        self.obs_key_sex = self._adata_ids_cellxgene.sex
        self.obs_key_organism = self._adata_ids_cellxgene.organism
        self.obs_key_state_exact = self._adata_ids_cellxgene.state_exact

        self.healthy_state_healthy = self._adata_ids_cellxgene.disease_state_healthy

        self.var_ensembl_col = self._adata_ids_cellxgene.gene_id_ensembl
        self.var_symbol_col = self._adata_ids_cellxgene.gene_id_names

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

        self.author = adata.uns[self._adata_ids_cellxgene.author][self._adata_ids_cellxgene.author_names]
        self.doi = adata.uns[self._adata_ids_cellxgene.doi]
        self.download_url_data = self.download_url_data
        self.id = self.id
        self.normalization = 'raw'
        self.organ = str(self.fn).split("_")[3]  # TODO interface this properly
        # self.organ = adata.obs["tissue"].values[0]
        self.organism = adata.obs[self._adata_ids_cellxgene.organism].values[0]
        self.assay_sc = adata.obs[self._adata_ids_cellxgene.assay].values[0]
        self.year = adata.uns[self._adata_ids_cellxgene.year]
