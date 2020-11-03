import anndata
import numpy as np
import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA, ADATA_IDS_CELLXGENE, META_DATA_FIELDS_CELLXGENE


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
        self.fn = fn
        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        fn = os.path.join(self.path, self.fn)
        adata = anndata.read(fn)
        adata.X = adata.raw.X

        self.adata.uns[ADATA_IDS_SFAIRA.author] = adata.uns[ADATA_IDS_CELLXGENE.author][ADATA_IDS_CELLXGENE.author_names]
        self.adata.uns[ADATA_IDS_SFAIRA.year] = adata.uns[ADATA_IDS_CELLXGENE.year]
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = adata.uns[ADATA_IDS_CELLXGENE.doi]
        if len(np.unique(adata.obs[ADATA_IDS_SFAIRA.species].values)) > 1:
            raise Warning("found multiple assay in data set %s" % self.fn)
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = adata.obs[ADATA_IDS_CELLXGENE.protocol].values[0]
        # Select tissue: blood is handled as a separate tissue in .obs
        #if len(np.unique(adata.obs["tissue"].values)) > 1:
        #    raise Warning("found multiple tissue in data set %s" % self.fn)
        #self.adata.uns["organ"] = adata.obs["tissue"].values[0]
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = str(self.fn).split("_")[3]
        if len(np.unique(adata.obs[ADATA_IDS_SFAIRA.species].values)) > 1:
            raise Warning("found multiple organisms in data set %s" % self.fn)
        self.adata.uns[ADATA_IDS_SFAIRA.species] = adata.obs[ADATA_IDS_CELLXGENE.species].values[0]
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[ADATA_IDS_SFAIRA.subtissue] = adata.obs[ADATA_IDS_CELLXGENE.subtissue].values
        self.adata.obs[ADATA_IDS_SFAIRA.dev_stage] = adata.obs[ADATA_IDS_CELLXGENE.dev_stage].values
        self.adata.obs[ADATA_IDS_SFAIRA.sex] = adata.obs[ADATA_IDS_CELLXGENE.sex].values
        self.adata.obs[ADATA_IDS_SFAIRA.ethnicity] = adata.obs[ADATA_IDS_CELLXGENE.ethnicity].values
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = adata.obs[ADATA_IDS_CELLXGENE.disease].values == ADATA_IDS_CELLXGENE.disease_state_healthy
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = adata.obs[ADATA_IDS_CELLXGENE.disease].values

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_id] = adata.obs[ADATA_IDS_CELLXGENE.cell_ontology_id].values.tolist()
        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = adata.obs[ADATA_IDS_CELLXGENE.cell_ontology_class].values.tolist()
        self.adata.obs[ADATA_IDS_SFAIRA.cell_types_original] = adata.obs[ADATA_IDS_CELLXGENE.cell_types_original].values.tolist()

        self._convert_and_set_var_names(
            symbol_col=ADATA_IDS_CELLXGENE.gene_id_names,
            ensembl_col=ADATA_IDS_CELLXGENE.gene_id_ensembl,
            new_index=ADATA_IDS_CELLXGENE.gene_id_ensembl
        )

