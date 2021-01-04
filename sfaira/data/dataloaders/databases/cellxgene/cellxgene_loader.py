import anndata
import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_CELLXGENE


class Dataset(DatasetBase):
    """
    This is a dataloader for downloaded h5ad from cellxgene.

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None],
            fn: str,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, **kwargs)
        self._ADATA_IDS_CELLXGENE = ADATA_IDS_CELLXGENE()
        self.fn = fn

        self.obs_key_cellontology_class = self._ADATA_IDS_CELLXGENE.cell_ontology_class
        self.obs_key_cellontology_id = self._ADATA_IDS_CELLXGENE.cell_ontology_id
        self.obs_key_cellontology_original = self._ADATA_IDS_CELLXGENE.cell_types_original
        self.obs_key_dev_stage = self._ADATA_IDS_CELLXGENE.dev_stage
        self.obs_key_ethnicity = self._ADATA_IDS_CELLXGENE.ethnicity
        self.obs_key_healthy = self._ADATA_IDS_CELLXGENE.healthy
        self.obs_key_sex = self._ADATA_IDS_CELLXGENE.sex
        self.obs_key_organism = self._ADATA_IDS_CELLXGENE.organism
        self.obs_key_state_exact = self._ADATA_IDS_CELLXGENE.state_exact
        self.obs_key_subtissue = self._ADATA_IDS_CELLXGENE.subtissue

        self.healthy_state_healthy = self._ADATA_IDS_CELLXGENE.disease_state_healthy

        self.var_ensembl_col = self._ADATA_IDS_CELLXGENE.gene_id_ensembl
        self.var_symbol_col = self._ADATA_IDS_CELLXGENE.gene_id_names

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        """
        Note that in contrast to data set specific data loaders, here, the core attributes are only identified from
        the data in this function and are not already set in the constructor. These attributes can still be
        used through meta data containers after the data was loaded once.

        :param fn:
        :return:
        """
        fn = os.path.join(self.path, self.fn)
        adata = anndata.read(fn)
        adata.X = adata.raw.X
        # TODO delete raw?

        self.author = adata.uns[self._ADATA_IDS_CELLXGENE.author][self._ADATA_IDS_CELLXGENE.author_names]
        self.doi = adata.uns[self._ADATA_IDS_CELLXGENE.doi]
        self.download = self.download
        self.id = self.id
        self.normalization = 'raw'
        self.organ = str(self.fn).split("_")[3]  # TODO interface this properly
        #self.organ = adata.obs["tissue"].values[0]
        self.organism = adata.obs[self._ADATA_IDS_CELLXGENE.organism].values[0]
        self.protocol = adata.obs[self._ADATA_IDS_CELLXGENE.protocol].values[0]
        self.year = adata.uns[self._ADATA_IDS_CELLXGENE.year]
