import anndata
import numpy as np
import os
import pandas as pd
import scipy.sparse
from typing import Union
import urllib.request
import zipfile

from .external import DatasetBase


class DatasetHcl(DatasetBase):
    """
    This is a dataloader template for human cell landscape data.
    """

    def __init__(
            self,
            path: Union[str, None],
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, **kwargs)

        self.obs_key_cellontology_class = "cell_ontology_class"
        self.obs_key_cellontology_original = "cell_ontology_class"
        self.obs_key_dev_stage = "dev_stage"
        self.obs_key_sex = "gender"
        self.obs_key_age = "age"

        self.author = 'Guo'
        self.doi = '10.1038/s41586-020-2157-4'
        self.organism = "human"
        self.download = 'https://ndownloader.figshare.com/files/17727365'
        self.download_meta = [
            'https://ndownloader.figshare.com/files/21758835',
            'https://ndownloader.figshare.com/files/22447898',
        ]
        self.healthy = True
        self.normalization = 'raw'
        self.protocol = 'microwell-seq'
        self.state_exact = 'healthy'
        self.year = 2020

        self.var_ensembl_col = "ensembl"
        self.var_symbol_col = "names"

    def _download(self):

        # download required files from human cell landscape publication data: https://figshare.com/articles/HCL_DGE_Data/7235471
        print(urllib.request.urlretrieve(
            'https://ndownloader.figshare.com/files/17727365',
            os.path.join(self.path, self.organism, self.organ, 'HCL_Fig1_adata.h5ad')
        ))
        print(urllib.request.urlretrieve(
            'https://ndownloader.figshare.com/files/21758835',
            os.path.join(self.path, self.organism, self.organ, 'HCL_Fig1_cell_Info.xlsx')
        ))
        print(urllib.request.urlretrieve(
            'https://ndownloader.figshare.com/files/22447898',
            os.path.join(self.path, self.organism, self.organ, 'annotation_rmbatch_data_revised417.zip')
        ))
        # extract the downloaded zip archive
        with zipfile.ZipFile(
                os.path.join(self.path, self.organism, self.organ, 'annotation_rmbatch_data_revised417.zip'),
                'r'
        ) as zip_ref:
            zip_ref.extractall(os.path.join(self.path, self.organism, self.organ))

    def _load_hcl(self, fn, sample_id: str):
        """
        Attempt to find file, cache entire HCL if file was not found.

        :param fn:
        :return:
        """
        adata = anndata.read(os.path.join(self.path, self.organism, self.organ, 'HCL_Fig1_adata.h5ad'))
        # convert to sparse matrix
        adata.X = scipy.sparse.csr_matrix(adata.X).copy()

        # harmonise annotations
        for col in ['batch', 'tissue']:
            adata.obs[col] = adata.obs[col].astype('str')
        adata.obs.index = adata.obs.index.str.replace('AdultJeJunum', 'AdultJejunum', regex=True).str.replace(
            'AdultGallBladder', 'AdultGallbladder', regex=True).str.replace(
            'FetalFemaleGonald', 'FetalFemaleGonad', regex=True)
        adata.obs.replace({'AdultJeJunum': 'AdultJejunum', 'AdultGallBladder': 'AdultGallbladder',
                           'FetalFemaleGonald': 'FetalFemaleGonad'}, regex=True, inplace=True)
        adata.obs.index = ["-".join(i.split('-')[:-1]) for i in adata.obs.index]

        # load celltype labels and harmonise them
        fig1_anno = pd.read_excel(
            os.path.join(self.path, self.organism, self.organ, 'HCL_Fig1_cell_Info.xlsx'),
            index_col='cellnames'
        )
        fig1_anno.index = fig1_anno.index.str.replace('AdultJeJunum', 'AdultJejunum', regex=True).str.replace(
            'AdultGallBladder', 'AdultGallbladder', regex=True).str.replace(
            'FetalFemaleGonald', 'FetalFemaleGonad', regex=True)

        # check that the order of cells and cell labels is the same
        assert np.all(fig1_anno.index == adata.obs.index)

        # add annotations to adata object and rename columns
        adata.obs = pd.concat([adata.obs, fig1_anno[['cluster', 'stage', 'donor', 'celltype']]], axis=1)
        adata.obs.columns = ['sample', 'tissue', 'n_genes', 'n_counts', 'cluster_global', 'stage', 'donor',
                             'celltype_global']

        # add sample-wise annotations to the full adata object
        df = pd.DataFrame(
            columns=['Cell_barcode', 'Sample', 'Batch', 'Cell_id', 'Cluster_id', 'Ages', 'Development_stage', 'Method',
                     'Gender', 'Source', 'Biomaterial', 'Name', 'ident', 'Celltype'])
        for f in os.listdir('annotation_rmbatch_data_revised417/'):
            df1 = pd.read_csv('annotation_rmbatch_data_revised417/' + f, encoding='unicode_escape')
            df = pd.concat([df, df1], sort=True)
        df = df.set_index('Cell_id')
        adata = adata[[i in df.index for i in adata.obs.index]].copy()
        a_idx = adata.obs.index.copy()
        adata.obs = pd.concat([adata.obs, df[['Ages', 'Celltype', 'Cluster_id', 'Gender', 'Method', 'Source']]], axis=1)
        assert np.all(a_idx == adata.obs.index)

        # remove mouse cells from the object
        adata = adata[adata.obs['Source'] != 'MCA2.0'].copy()

        # tidy up the column names of the obs annotations
        adata.obs.columns = ['sample', 'sub_tissue', 'n_genes', 'n_counts', 'cluster_global', 'dev_stage',
                             'donor', 'celltype_global', 'age', 'celltype_specific', 'cluster_specific', 'gender',
                             'protocol', 'source']

        # create some annotations that are used in sfaira
        adata.obs["healthy"] = True
        adata.obs["state_exact"] = 'healthy'
        adata.obs["cell_ontology_class"] = adata.obs["celltype_global"]
        adata.obs["cell_ontology_id"] = None
        adata.var = adata.var.reset_index().rename({'index': 'names'}, axis='columns')

        # create a tidy organ annotation which is then used in sfaira
        adata.obs['organ'] = adata.obs['sub_tissue'] \
            .str.replace("Adult", "") \
            .str.replace("Fetal", "") \
            .str.replace("Neonatal", "") \
            .str.replace("Transverse", "") \
            .str.replace("Sigmoid", "") \
            .str.replace("Ascending", "") \
            .str.replace("Cord", "") \
            .str.replace("Peripheral", "") \
            .str.replace("CD34P", "") \
            .str.replace("Cerebellum", "Brain") \
            .str.replace("TemporalLobe", "Brain") \
            .str.replace("BoneMarrow", "Bone") \
            .str.replace("Spinal", "SpinalCord") \
            .str.replace("Intestine", "Stomach") \
            .str.replace("Eyes", "Eye") \
            .str.lower()

        self.adata = adata[adata.obs['sample'] == sample_id].copy()
