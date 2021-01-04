import anndata
import os
from typing import Union
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):
    """
    This data loader directly processes the two raw data files which can be obtained from the `download_website`
    and `download_website_meta` attributes of this class.

    :param path:
    :param meta_path:
    :param kwargs:
    """

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.organism = "human"
        self.id = "human_lung_2019_dropseq_braga_003_10.1038/s41591-019-0468-5"
        self.download = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130148/suppl/GSE130148%5Fraw%5Fcounts%2Ecsv%2Egz"
        self.download_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130148/suppl/GSE130148%5Fbarcodes%5Fcell%5Ftypes%2Etxt%2Egz"
        self.organ = "lung"
        self.sub_tissue = "parenchymal lung and distal airway specimens"
        self.author = 'Teichmann'
        self.year = 2019
        self.doi = "10.1038/s41591-019-0468-5"
        self.protocol = 'dropseq'
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = 'uninvolved areas of tumour resection material'
        self.var_symbol_col = 'index'
        self.obs_key_cellontology_original = 'celltype'

        self.class_maps = {
            "0": {
                'Fibroblast': 'Fibroblasts',
                'Type 2': 'AT2',
                'B cell': 'B cell lineage',
                'Macrophages': 'Macrophages',
                'NK cell': 'Innate lymphoid cells',
                'T cell': 'T cell lineage',
                'Ciliated': 'Multiciliated lineage',
                'Lymphatic': 'Lymphatic EC',
                'Type 1': 'AT1',
                'Transformed epithelium': '1_Epithelial',
                'Secretory': 'Secretory',
                'Endothelium': '1_Endothelial',
                'Mast cell': 'Mast cells',
            },
        }

    def _load(self, fn=None):
        if fn is None:
            fn = [
                os.path.join(self.path, "human", "lung", "GSE130148_raw_counts.csv.gz"),
                os.path.join(self.path, "human", "lung", "GSE130148_barcodes_cell_types.txt.gz"),
            ]
        self.adata = anndata.read_csv(fn[0]).T
        self.adata.obs = pd.read_csv(fn[1], sep='\t', index_col=0)

        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
