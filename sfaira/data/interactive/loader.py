import anndata
from typing import Union
from .external import DatasetBase


class DatasetInteractive(DatasetBase):

    def __init__(
            self,
            data: anndata.AnnData,
            organism: str,
            organ: str,
            gene_symbol_col: Union[str, None] = 'index',
            gene_ens_col: Union[str, None] = None,
            class_maps: dict = {},
            dataset_id: str = "interactive",
            path: Union[str, None] = ".",
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
        :param class_maps: Cell type class maps.
        :param dataset_id: Identifer of data set.
        :param path:
        :param meta_path:
        :param cache_path:
        """
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path)
        self.adata = data

        self.id = dataset_id
        self.organ = organ
        self.organism = organism

        self.var_symbol_col = gene_symbol_col
        self.var_ensembl_col = gene_ens_col

        self.class_maps = class_maps

    def _load(self, fn=None):
        pass
