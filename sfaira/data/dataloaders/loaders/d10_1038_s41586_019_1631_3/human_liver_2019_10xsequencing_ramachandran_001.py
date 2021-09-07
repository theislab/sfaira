import os
import numpy as np

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    TODO move state exact to disease, condition Uninjured is healthy
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://datashare.is.ed.ac.uk/bitstream/handle/10283/3433/tissue.rdata"
        self.download_url_meta = None

        self.assay_sc = "10x 3' v2"
        self.author = "Ramachandran"
        self.doi_journal = "10.1038/s41586-019-1631-3"
        self.normalization = "raw"
        self.organ = "liver"
        self.organism = "human"
        self.sample_source = "primary_tissue"
        self.year = 2019

        self.gene_id_symbols_var_key = "index"

        self.cell_type_obs_key = "annotation_lineage"
        self.state_exact_obs_key = "condition"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    import anndata2ri
    from rpy2.robjects import r

    fn = os.path.join(data_dir, "tissue.rdata")
    anndata2ri.activate()  # TODO: remove global activation of anndata2ri and use localconverter once it's fixed
    adata = r(
        f"library(Seurat)\n"
        f"load('{fn}')\n"
        f"new_obj = CreateSeuratObject(counts = tissue@raw.data)\n"
        f"new_obj@meta.data = tissue@meta.data\n"
        f"as.SingleCellExperiment(new_obj)\n"
    )
    adata.obs["nGene"] = adata.obs["nGene"].astype(np.int32)

    return adata
