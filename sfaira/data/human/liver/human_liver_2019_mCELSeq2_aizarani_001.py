import anndata
import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
import pandas as pd


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
        self.id = "human_liver_2019_mCELSeq2_aizarani_001_10.1038/s41586-019-1373-2"
        self.download_website = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395%5FNormalhumanlivercellatlasdata%2Etxt%2Egz"
        self.download_website_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395%5Fclusterpartition%2Etxt%2Egz"
        self.organ = "liver"
        self.sub_tissue = "liver"
        self.annotated = True

        self.class_maps = {
            "0": {
                '1': 'NK, NKT and T cells',
                '2': 'Kupffer Cell',
                '3': 'NK, NKT and T cells',
                '4': 'Cholangiocytes',
                '5': 'NK, NKT and T cells',
                '6': 'Kupffer Cell',
                '7': 'Cholangiocytes',
                '8': 'B Cell',
                '9': 'Liver sinusoidal endothelial cells',
                '10': 'Macrovascular endothelial cells',
                '11': 'Hepatocyte',
                '12': 'NK, NKT and T cells',
                '13': 'Liver sinusoidal endothelial cells',
                '14': 'Hepatocyte',
                '15': 'Other endothelial cells',
                '16': 'Unknown',
                '17': 'Hepatocyte',
                '18': 'NK, NKT and T cells',
                '19': 'Unknown',
                '20': 'Liver sinusoidal endothelial cells',
                '21': 'Macrovascular endothelial cells',
                '22': 'B Cell',
                '23': 'Kupffer Cell',
                '24': 'Cholangiocytes',
                '25': 'Kupffer Cell',
                '26': 'Other endothelial cells',
                '27': 'Unknown',
                '28': 'NK, NKT and T cells',
                '29': 'Macrovascular endothelial cells',
                '30': 'Hepatocyte',
                '31': 'Kupffer Cell',
                '32': 'Liver sinusoidal endothelial cells',
                '33': 'Hepatic stellate cells',
                '34': 'B Cell',
                '35': 'Other endothelial cells',
                '36': 'Unknown',
                '37': 'Unknown',
                '38': 'B Cell',
                '39': 'Cholangiocytes'
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = [
                    os.path.join(self.path, "human", "liver", "GSE124395_Normalhumanlivercellatlasdata.txt.gz"),
                    os.path.join(self.path, "human", "liver", "GSE124395_clusterpartition.txt.gz")
                ]
            self.adata = anndata.AnnData(pd.read_csv(fn[0], sep='\t').T)
            celltype_df = pd.read_csv(fn[1], sep=' ')
            self.adata = self.adata[[i in celltype_df.index for i in self.adata.obs.index]].copy()
            self.adata.obs['CellType'] = [str(celltype_df.loc[i]['sct@cpart']) for i in self.adata.obs.index]

        self.adata.uns[ADATA_IDS_SFAIRA.author] = 'Gruen'
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = '10.1038/s41586-019-1373-2'
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = 'mCEL-Seq2'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs['cell_ontology_class'] = self.adata.obs['CellType']
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None)
