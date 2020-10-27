import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
import anndata
import pandas as pd


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
            **kwargs
    ):
        DatasetBase.__init__(self=self, path=path, meta_path=meta_path, **kwargs)
        self.species = "human"
        self.id = "human_lung_2019_dropseq_braga_003_10.1038/s41591-019-0468-5"
        self.download_website = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130148/suppl/GSE130148%5Fraw%5Fcounts%2Ecsv%2Egz"
        self.download_website_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130148/suppl/GSE130148%5Fbarcodes%5Fcell%5Ftypes%2Etxt%2Egz"
        self.organ = "lung"
        self.sub_tissue = "parenchymal lung and distal airway specimens"
        self.has_celltypes = True

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
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = [
                    os.path.join(self.path, "human/lung/GSE130148_raw_counts.csv.gz"),
                    os.path.join(self.path, "human/lung/GSE130148_barcodes_cell_types.txt.gz"),
                ]
            self.adata = anndata.read_csv(fn[0]).T
            self.adata.obs = pd.read_csv(fn[1], sep='\t', index_col=0)

        self.adata.uns[ADATA_IDS.author] = 'Teichmann'
        self.adata.uns[ADATA_IDS.year] = 2019
        self.adata.uns[ADATA_IDS.doi] = "10.1038/s41591-019-0468-5"
        self.adata.uns[ADATA_IDS.protocol] = 'dropseq'
        self.adata.uns[ADATA_IDS.organ] = self.organ
        self.adata.uns[ADATA_IDS.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS.animal] = "human"
        self.adata.uns[ADATA_IDS.id] = self.id
        self.adata.uns[ADATA_IDS.wget_download] = [self.download_website, self.download_website_meta]
        self.adata.uns[ADATA_IDS.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS.normalization] = 'raw'

        self.adata.obs[ADATA_IDS.cell_ontology_class] = self.adata.obs['celltype']
        self.set_unkown_class_id(ids=["1_Unicorns and artifacts"])
        self.adata.obs[ADATA_IDS.healthy] = True
        self.adata.obs['state_exact'] = 'uninvolved areas of tumour resection material'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS.gene_id_ensembl)
