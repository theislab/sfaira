import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
import anndata
import pandas as pd


class Dataset(DatasetBase):
    """
    This data loader directly processes the raw data files which can be obtained from the `download_website` and
    `download_website_meta` attributes of this class.

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
        self.id = "human_kidney_2019_10xSn_lake_001_10.1038/s41467-019-10861-2"
        self.download_website = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121862/suppl/GSE121862%5FUCSD%2DWU%5FSingle%5FNuclei%5FCluster%5FAnnotated%5FRaw%5FUMI%5FMatrix%2Etsv%2Egz"
        self.download_website_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121862/suppl/GSE121862%5FUCSD%2DWU%5FSingle%5FNuclei%5FCluster%5FAnnotations%2Ecsv%2Egz"
        self.organ = "kidney"
        self.sub_tissue = "kidney"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'Collecting Duct - Intercalated Cells Type A (cortex)': 'Collecting Duct - Intercalated Cells Type A (cortex)',
                'Collecting Duct - Intercalated Cells Type A (medulla)': 'Collecting Duct - Intercalated Cells Type A (medulla)',
                'Collecting Duct - Intercalated Cells Type B': 'Collecting Duct - Intercalated Cells Type B',
                'Collecting Duct - PCs - Stressed Dissoc Subset': 'Collecting Duct - PCs - Stressed Dissoc Subset',
                'Collecting Duct - Principal Cells (cortex)': 'Collecting Duct - Principal Cells (cortex)',
                'Collecting Duct - Principal Cells (medulla)': 'Collecting Duct - Principal Cells (medulla)',
                'Connecting Tubule': 'Connecting tubule',
                'Decending Limb': 'Decending Limb',
                'Distal Convoluted Tubule': 'Distal Convoluted Tubule',
                'Endothelial Cells (unassigned)': 'Endothelial Cells (unassigned)',
                'Endothelial Cells - AEA & DVR ': 'Endothelial Cells - AEA & DVR',
                'Endothelial Cells - AVR': 'Endothelial Cells - AVR',
                'Endothelial Cells - glomerular capillaries': 'Endothelial Cells - glomerular capillaries',
                'Epithelial Cells (unassigned)': 'Epithelial Cells (unassigned)',
                'Immune Cells - Macrophages': 'Macrophage',
                'Interstitium': 'Interstitium',
                'Mesangial Cells': 'Mesangial Cells',
                'Podocytes': 'Podocyte',
                'Proximal Tubule Epithelial Cells (S1)': 'Proximal Tubule Epithelial Cells (S1)',
                'Proximal Tubule Epithelial Cells (S2)': 'Proximal Tubule Epithelial Cells (S2)',
                'Proximal Tubule Epithelial Cells (S3)': 'Proximal Tubule Epithelial Cells (S3)',
                'Proximal Tubule Epithelial Cells - Fibrinogen+ (S3 )': 'Proximal Tubule Epithelial Cells - Fibrinogen+ (S3)',
                'Proximal Tubule Epithelial Cells - Stress/Inflam': 'Proximal Tubule Epithelial Cells - Stress/Inflam',
                'Thick Ascending Limb': 'Thick ascending limb of Loop of Henle',
                'Thin ascending limb': 'Thin ascending limb',
                'Unknown - Novel PT CFH+ Subpopulation (S2)': 'Unknown - Novel PT CFH+ Subpopulation (S2)',
                'Vascular Smooth Muscle Cells and pericytes': 'Vascular Smooth Muscle Cells and pericytes',
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = [
                    os.path.join(self.path, "human", "kidney", "GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotated_Raw_UMI_Matrix.tsv.gz"),
                    os.path.join(self.path, "human", "kidney", "GSE121862_UCSD-WU_Single_Nuclei_Cluster_Annotations.csv.gz")
                ]
            self.adata = anndata.AnnData(pd.read_csv(fn[0], sep='\t').T)
            annot = pd.read_csv(fn[1], index_col=0, dtype='category')
            self.adata.obs['celltype'] = [annot.loc[i.split('_')[0][1:]]['Annotation'] for i in self.adata.obs.index]

        self.adata.uns[ADATA_IDS_SFAIRA.author] = 'Jain'
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = '10.1038/s41467-019-10861-2'
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = '10xSn'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.has_celltypes
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['celltype']
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS_SFAIRA.gene_id_ensembl)
