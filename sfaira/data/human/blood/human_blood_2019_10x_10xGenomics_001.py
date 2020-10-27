import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS
import anndata


class Dataset(DatasetBase):
    """
    This data loader requires manual preprocessing of the raw datafile. To download the data, use the link in the
    `.download_website` attribute of this class. To create the file required by this dataloader, run the following
    python code:

    import scanpy
    scanpy.read_10x_h5('pbmc_10k_v3_filtered_feature_bc_matrix.h5').write('pbmc_10k_v3_filtered_feature_bc_matrix.h5ad')

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_blood_2019_10x_10xGenomics_001_unknown"
        self.download_website = "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
        self.organ = "blood"
        self.sub_tissue = "pbmcs"
        self.has_celltypes = False

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human/blood/pbmc_10k_v3_filtered_feature_bc_matrix.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[ADATA_IDS.lab] = '10x Genomics'
        self.adata.uns[ADATA_IDS.year] = 2019
        self.adata.uns[ADATA_IDS.doi] = None
        self.adata.uns[ADATA_IDS.protocol] = '10x'
        self.adata.uns[ADATA_IDS.organ] = self.organ
        self.adata.uns[ADATA_IDS.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS.animal] = "human"
        self.adata.uns[ADATA_IDS.id] = self.id
        self.adata.uns[ADATA_IDS.wget_download] = "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
        self.adata.uns[ADATA_IDS.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS.normalization] = 'raw'

        self.adata.obs[ADATA_IDS.cell_ontology_class] = None
        self.adata.obs[ADATA_IDS.healthy] = True
        self.adata.obs[ADATA_IDS.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col='gene_ids', new_index=ADATA_IDS.gene_id_ensembl)
