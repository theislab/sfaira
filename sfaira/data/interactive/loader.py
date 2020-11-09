import anndata
from typing import Union
from .external import DatasetBase


class Dataset(DatasetBase):

    def __init__(
            self,
            data: anndata.AnnData,
            species: str,
            organ: str,
            gene_symbol_col: Union[str, None] = 'index',
            gene_ensg_col: Union[str, None] = None,
            class_maps: dict = {},
            dataset_id: str = "interactive",
            **kwargs
    ):
        """

        :param data:
        :param species:
        :param organ:
        :param class_maps:
        :param id:
        :param kwargs:
        """
        DatasetBase.__init__(self=self, path=None, meta_path=None, **kwargs)
        self.adata = data
        self.species = species
        self.id = dataset_id
        self.organ = organ

        self.gene_symbol_col = gene_symbol_col
        self.gene_ensg_col = gene_ensg_col

        self.class_maps = class_maps

    def _load(self, fn=None):
        self._convert_and_set_var_names(
            symbol_col=self.gene_symbol_col,
            ensembl_col=self.gene_ensg_col,
            new_index='ensembl'
        )
