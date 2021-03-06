import os
from typing import Union
from .external import DatasetBase
import anndata
import scipy.sparse


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data file which can be obtained from the `download_website` attribute of
    this class.

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
        self.id = "human_spleen_2019_10x_madissoon_001_10.1101/741405"
        self.download_website = "https://cellgeni.cog.sanger.ac.uk/tissue-stability/tissue-stability/spleen.cellxgene.h5ad"
        self.download_website_meta = None
        self.organ = "spleen"
        self.sub_tissue = "spleen"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                "B_Hypermutation": "B_Hypermutation",
                "B_T_doublet": "B_T_doublet",
                "B_follicular": "B_follicular",
                "B_mantle": "B_mantle",
                "CD34_progenitor": "CD34_progenitor",
                "DC_1": "DC_1",
                "DC_2": "DC_2",
                "DC_activated": "DC_activated",
                "DC_plasmacytoid": "DC_plasmacytoid",
                "ILC": "ILC",
                "Macrophage": "Macrophage",
                "Monocyte": "Monocyte",
                "NK_CD160pos": "NK_CD160pos",
                "NK_FCGR3Apos": "NK_FCGR3Apos",
                "NK_dividing": "NK_dividing",
                "Plasma_IgG": "Plasma_IgG",
                "Plasma_IgM": "Plasma_IgM",
                "Plasmablast": "Plasmablast",
                "Platelet": "Platelet",
                "T_CD4_conv": "T_CD4_conv",
                "T_CD4_fh": "T_CD4_fh",
                "T_CD4_naive": "T_CD4_naive",
                "T_CD4_reg": "T_CD4_reg",
                "T_CD8_CTL": "T_CD8_CTL",
                "T_CD8_MAIT": "T_CD8_MAIT",
                "T_CD8_activated": "T_CD8_activated",
                "T_CD8_gd": "T_CD8_gd",
                "T_cell_dividing": "Proliferating T cell",
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "spleen", "spleen.cellxgene.h5ad")
            self.adata = anndata.read(fn)
            self.adata.X = self.adata.X.multiply(scipy.sparse.csc_matrix(self.adata.obs['n_counts'].values[:, None]))\
                                       .multiply(1/10000)

        self.adata.uns["lab"] = "Meyer"
        self.adata.uns["year"] = 2019
        self.adata.uns["doi"] = "10.1101/741405"
        self.adata.uns["protocol"] = "10x"
        self.adata.uns["organ"] = self.organ
        self.adata.uns["subtissue"] = self.sub_tissue
        self.adata.uns["animal"] = "human"
        self.adata.uns["id"] = self.id
        self.adata.uns["wget_download"] = self.download_website
        self.adata.uns["has_celltypes"] = self.has_celltypes
        self.adata.uns["counts"] = 'raw'

        self.adata.obs["cell_ontology_class"] = self.adata.obs['Celltypes']
        self.set_unkown_class_id(ids=["Unknown"])
        self.adata.obs["healthy"] = True
        self.adata.obs["state_exact"] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col='gene_ids-HCATisStab7463846',
                                        new_index='ensembl')
