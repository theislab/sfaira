import os

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = "https://datashare.is.ed.ac.uk/bitstream/handle/10283/3433/tissue.rdata"
        self.download_url_meta = None

        self.author = "Ramachandran"
        self.doi = "10.1038/s41586-019-1631-3"
        self.normalization = "raw"
        self.organ = "liver"
        self.organism = "human"
        self.protocol = "10X sequencing"
        self.year = 2019

        self.var_symbol_col = "index"

        self.cellontology_original_obs_key = "annotation_lineage"
        self.state_exact_obs_key = "condition"
        self.key_healthy_obs_key = "condition"
        self.healthy_state_healthy = "Uninjured"

        self.set_dataset_id(idx=1)

    def _load(self):
        import anndata2ri
        from rpy2.robjects import r

        fn = os.path.join(self.data_dir, "tissue.rdata")
        anndata2ri.activate()  # TODO: remove global activation of anndata2ri and use localconverter once it's fixed
        adata = r(
            f"library(Seurat)\n"
            f"load('{fn}')\n"
            f"new_obj = CreateSeuratObject(counts = tissue@raw.data)\n"
            f"new_obj@meta.data = tissue@meta.data\n"
            f"as.SingleCellExperiment(new_obj)\n"
        )

        return adata
