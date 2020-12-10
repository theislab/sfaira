import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
import anndata
import pandas as pd


class Dataset(DatasetBase):
    """
    The input files for this dataloader (GSE115469.csv.gz and GSE115469_labels.txt) were kindly provided to us by the
    authors of the publication. Please contact them directly to obtain the required
    files.

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
        self.id = "human_liver_2018_10x_macparland_001_10.1038/s41467-018-06318-7"
        self.download_website = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469"
        self.download_website_meta = 'private'
        self.organ = "liver"
        self.sub_tissue = "caudate lobe"
        self.annotated = True

        self.class_maps = {
            "0": {
                '1':'Hepatocyte 1',
                '2':'Alpha beta T cells',
                '3':'Hepatocyte 2',
                '4':'Inflammatory macrophages',
                '5':'Hepatocyte 3',
                '6':'Hepatocyte 4',
                '7':'Plasma cells',
                '8':'NK cell',
                '9':'Gamma delta T cells 1',
                '10':'Non inflammatory macrophages',
                '11':'Periportal LSECs',
                '12':'Central venous LSECs',
                '13':'Endothelial cell',
                '14':'Hepatocyte 5',
                '15':'Hepatocyte 6',
                '16':'Mature B cells',
                '17':'Cholangiocytes',
                '18':'Gamma delta T cells 2',
                '19':'Erythroid cells',
                '20':'Hepatic stellate cells'
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = [
                    os.path.join(self.path, "human", "liver", "GSE115469.csv.gz"),
                    os.path.join(self.path, "human", "liver", "GSE115469_labels.txt")
                ]
            self.adata = anndata.read_csv(fn[0]).T
            celltype_df = pd.read_csv(fn[1], sep='\t').set_index('CellName')
            self.adata.obs['celltype'] = [str(celltype_df.loc[i]['Cluster#']) for i in self.adata.obs.index]

        self.adata.uns[ADATA_IDS_SFAIRA.author] = 'McGilvray'
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2018
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = '10.1038/s41467-018-06318-7'
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['celltype']
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
