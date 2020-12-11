import anndata
import numpy as np
import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA, ADATA_IDS_CELLXGENE


class DatasetCellxgene(DatasetBase):
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
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self._ADATA_IDS_CELLXGENE = ADATA_IDS_CELLXGENE()
        self.fn = fn

        self.obs_key_ethnicity = self._ADATA_IDS_CELLXGENE.ethnicity
        self.obs_key_healthy = self._ADATA_IDS_CELLXGENE.healthy
        self.obs_key_species = self._ADATA_IDS_CELLXGENE.species
        self.obs_key_species = self._ADATA_IDS_CELLXGENE.species
        self.obs_key_sex = self._ADATA_IDS_CELLXGENE.sex
        self.obs_key_species = self._ADATA_IDS_CELLXGENE.species
        self.obs_key_subtissue = self._ADATA_IDS_CELLXGENE.subtissue
        self.obs_key_species = self._ADATA_IDS_CELLXGENE.species
        self.obs_key_state_exact = self._ADATA_IDS_CELLXGENE.disease

        self.healthy_state_healthy = self._ADATA_IDS_CELLXGENE.disease_state_healthy

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        fn = os.path.join(self.path, self.fn)
        adata = anndata.read(fn)
        adata.X = adata.raw.X

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = adata.uns[self._ADATA_IDS_CELLXGENE.author][self._ADATA_IDS_CELLXGENE.author_names]
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = adata.uns[self._ADATA_IDS_CELLXGENE.year]
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = adata.uns[self._ADATA_IDS_CELLXGENE.doi]
        if len(np.unique(adata.obs[ADATA_IDS_SFAIRA.species].values)) > 1:
            raise Warning("found multiple assay in data set %s" % self.fn)
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = adata.obs[self._ADATA_IDS_CELLXGENE.protocol].values[0]
        # Select tissue: blood is handled as a separate tissue in .obs
        #if len(np.unique(adata.obs["tissue"].values)) > 1:
        #    raise Warning("found multiple tissue in data set %s" % self.fn)
        #self.adata.uns["organ"] = adata.obs["tissue"].values[0]
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = str(self.fn).split("_")[3]
        if len(np.unique(adata.obs[ADATA_IDS_SFAIRA.species].values)) > 1:
            raise Warning("found multiple organisms in data set %s" % self.fn)
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = adata.obs[self._ADATA_IDS_CELLXGENE.species].values[0]
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[self._ADATA_IDS_SFAIRA.dev_stage] = adata.obs[self._ADATA_IDS_CELLXGENE.dev_stage].values

        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_id] = adata.obs[self._ADATA_IDS_CELLXGENE.cell_ontology_id].values.tolist()
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] = adata.obs[self._ADATA_IDS_CELLXGENE.cell_ontology_class].values.tolist()
        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_types_original] = adata.obs[self._ADATA_IDS_CELLXGENE.cell_types_original].values.tolist()

        self._convert_and_set_var_names(
            symbol_col=self._ADATA_IDS_CELLXGENE.gene_id_names,
            ensembl_col=self._ADATA_IDS_CELLXGENE.gene_id_ensembl
        )

