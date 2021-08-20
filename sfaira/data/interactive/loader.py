import anndata
from typing import Union

from sfaira.data import DatasetBase


class DatasetInteractive(DatasetBase):

    def __init__(
            self,
            data: anndata.AnnData,
            organism: str,
            organ: str,
            gene_symbol_col: Union[str, None] = 'index',
            gene_ens_col: Union[str, None] = None,
            obs_key_celltypes: Union[str, None] = None,
            class_maps: dict = {},
            dataset_id: str = "interactive_dataset",
            data_path: Union[str, None] = ".",
            meta_path: Union[str, None] = ".",
            cache_path: Union[str, None] = ".",
    ):
        """
        Load data set into sfaira data format.

        :param data: Data set.
        :param organism: Organism of data set.
        :param organ: Organ of data set.
        :param gene_symbol_col: Column name in .var which contains gene symbols. Set to "index" to use the index.
        :param gene_ens_col:  Column name in .var which contains ENSG symbols. Set to "index" to use the index.
        :param obs_key_celltypes: .obs column name which contains cell type labels.
        :param class_maps: Cell type class maps.
        :param dataset_id: Identifer of data set.
        :param data_path:
        :param meta_path:
        :param cache_path:
        """
        super().__init__(data_path=data_path, meta_path=meta_path, cache_path=cache_path)
        self.id = dataset_id

        self.author = "interactive_dataset"
        self.doi_journal = "interactive_dataset"

        self.download_url_data = "."
        self.download_url_meta = "."

        # self.age  # not currently supported
        # self.assay_sc  # not currently supported
        # self.assay_differentiation  # not currently supported
        # self.assay_type_differentiation  # not currently supported
        # self.cell_line  # not currently supported
        # self.dev_stage  # not currently supported
        # self.ethnicity  # not currently supported
        # self.healthy  # not currently supported
        # self.normalisation  # not currently supported
        self.organ = organ
        self.organism = organism
        # self.sample_source # not currently supported
        # self.sex  # not currently supported
        # self.state_exact  # not currently supported
        # self.year  # not currently supported

        self.obs_key_cell_types_original = obs_key_celltypes

        # self.obs_key_age  # not currently supported
        # self.obs_key_assay_sc  # not currently supported
        # self.obs_key_assay_differentiation  # not currently supported
        # self.obs_key_assay_type_differentiation  # not currently supported
        # self.obs_key_cell_line  # not currently supported
        # self.obs_key_dev_stage  # not currently supported
        # self.obs_key_ethnicity  # not currently supported
        # self.obs_key_healthy  # not currently supported
        # self.obs_key_organ  # not currently supported
        # self.obs_key_organism  # not currently supported
        # self.obs_key_sample_source  # not currently supported
        # self.obs_key_sex  # not currently supported
        # self.obs_key_state_exact  # not currently supported

        self.gene_id_symbols_var_key = gene_symbol_col
        self.gene_id_ensembl_var_key = gene_ens_col

        self.class_maps = class_maps

        self.adata = data

    def _load(self):
        pass
