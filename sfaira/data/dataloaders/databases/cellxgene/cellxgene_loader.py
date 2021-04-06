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
        self._adata_ids_cellxgene = AdataIdsCellxgene()
        self.fn = fn

        self.cellontology_class_obs_key = self._adata_ids_cellxgene.cell_ontology_class
        self.cellontology_id_obs_key = self._adata_ids_cellxgene.cell_ontology_id
        self.cellontology_original_obs_key = self._adata_ids_cellxgene.cell_types_original
        self.development_stage_obs_key = self._adata_ids_cellxgene.development_stage
        self.disease_obs_key = self._adata_ids_cellxgene.disease
        self.ethnicity_obs_key = self._adata_ids_cellxgene.ethnicity
        self.sex_obs_key = self._adata_ids_cellxgene.sex
        self.organ_obs_key = self._adata_ids_cellxgene.organism
        self.state_exact_obs_key = self._adata_ids_cellxgene.state_exact

        self.var_ensembl_col = self._adata_ids_cellxgene.gene_id_ensembl
        self.var_symbol_col = self._adata_ids_cellxgene.gene_id_names

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
        self.assay_sc = adata.obs[self._adata_ids_cellxgene.assay_sc].values[0]
        self.year = adata.uns[self._adata_ids_cellxgene.year]
