import anndata
from typing import Union
from sfaira.data import DatasetBase


class Dataset_d10_1101_661728(DatasetBase):
    """
    This is a dataloader template for tabula muris data.
    """

    def __init__(
            self,
            path: Union[str, None],
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            source: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.source = source
        if self.source == "aws":
            self.download = "https://czb-tabula-muris-senis.s3-us-west-2.amazonaws.com/Data-objects/"
        elif self.source == "figshare":
            self.download = "https://ndownloader.figshare.com/articles/8273102/versions/2"
        else:
            raise ValueError("source %s not recognized" % self.source)

        self.obs_key_cellontology_class = "cell_ontology_class"
        self.obs_key_cellontology_id = "cell_ontology_id"
        self.obs_key_cellontology_original = "free_annotation"
        self.obs_key_age = "age"
        self.obs_key_dev_stage = "development_stage"  # not given in all data sets
        self.obs_key_sex = "sex"
        self.obs_key_subtissue = "subtissue"

        self.author = "Quake"
        self.doi = "10.1101/661728"
        self.healthy = True
        self.normalization = 'norm'
        self.organism = "mouse"
        self.state_exact = "healthy"
        self.year = 2019

        self.var_ensembl_col = None
        self.var_symbol_col = "index"

    def _load_generalized(self, fn):
        self.adata = anndata.read_h5ad(fn)
        if self.source == "aws":
            self.adata.X = self.adata.raw.X
            self.adata.var = self.adata.raw.var
            del self.adata.raw
            self.adata.obsm = {}
            self.adata.varm = {}
            self.adata.uns = {}

    def _get_protocol_tms(self, x) -> str:
        return "smartseq2" if "smartseq2" in x else "10x"
