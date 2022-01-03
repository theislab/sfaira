import anndata
from typing import Union

from sfaira.data import DatasetBase


class DatasetInteractive(DatasetBase):

    def __init__(
            self,
            data: anndata.AnnData,
            gene_symbol_col: Union[str, None] = 'index',
            gene_ens_col: Union[str, None] = None,
            dataset_id: str = "interactive_dataset",
            data_path: Union[str, None] = ".",
            meta_path: Union[str, None] = ".",
            cache_path: Union[str, None] = ".",
    ):
        """
        Load data set into sfaira data format.

        :param data: Data set.
        :param gene_symbol_col: Column name in .var which contains gene symbols. Set to "index" to use the index.
        :param gene_ens_col:  Column name in .var which contains ENSG symbols. Set to "index" to use the index.
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

        self.gene_id_symbols_var_key = gene_symbol_col
        self.gene_id_ensembl_var_key = gene_ens_col

        self.adata = data

    def _load(self):
        pass
