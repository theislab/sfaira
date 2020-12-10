import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
import anndata
import pandas as pd


class Dataset(DatasetBase):
    """
    This data loader supports reading of the downloaded raw data file if `load_raw=True` is passed to self.load()
    To download the datafile required by this dataloader, use the link provided as the `download_website` attribute of
    this class and obtain cell type annotations ('hc_meta_data_stromal_with_donor.txt' and
    'uc_meta_data_stromal_with_donor.txt') directly from the authors of the paper. For (up
    to 100-fold faster) repeated data loading, please pass `load_raw=False` when calling the self.load() method. For
    this, you need to preprocess the raw files as below and place the resulting h5ad file in the data folder of this
    organ:

    import anndata
    import pandas as pd

    adata = anndata.read_loom('f8aa201c-4ff1-45a4-890e-840d63459ca2.homo_sapiens.loom')
    ctuc = pd.read_csv('uc_meta_data_stromal_with_donor.txt', sep='\t')
    cthealthy = pd.read_csv('hc_meta_data_stromal_with_donor.txt', sep='\t')

    adata = adata[adata.obs['emptydrops_is_cell'] == 't'].copy()
    adata = adata[adata.X.sum(axis=1).flatten() >= 250].copy()

    uc = adata[adata.obs['donor_organism.diseases.ontology_label'] == "ulcerative colitis (disease)"].copy()
    bcuc = [i.split('-')[0] for i in ctuc['Barcode']]
    seluc = []
    for i in uc.obs['barcode']:
        seluc.append((uc.obs['barcode'].str.count(i).sum() == 1) and i in bcuc)
    uc = uc[seluc].copy()
    ctuc.index = [i.split('-')[0] for i in ctuc['Barcode']]
    uc.obs['celltype'] = [ctuc.loc[i]['Cluster'] for i in uc.obs['barcode']]
    uc.var = uc.var.reset_index().rename(columns={'index': 'names'}).set_index('featurekey')

    healthy = adata[adata.obs['donor_organism.diseases.ontology_label'] == "normal"].copy()
    bchealthy = [i.split('-')[0] for i in cthealthy['Barcode']]
    selhealthy = []
    for i in healthy.obs['barcode']:
        selhealthy.append((healthy.obs['barcode'].str.count(i).sum() == 1) and i in bchealthy)
    healthy = healthy[selhealthy].copy()
    cthealthy.index = [i.split('-')[0] for i in cthealthy['Barcode']]
    healthy.obs['celltype'] = [cthealthy.loc[i]['Cluster'] for i in healthy.obs['barcode']]
    healthy.var = healthy.var.reset_index().rename(columns={'index': 'names'}).set_index('featurekey')

    adata = healthy.concatenate(uc)
    adata.write('kinchenetal.h5ad')

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
        self.id = "human_colon_2019_10x_kinchen_001_10.1016/j.cell.2018.08.067"
        self.download_website = "https://data.humancellatlas.org/project-assets/project-matrices/f8aa201c-4ff1-45a4-890e-840d63459ca2.homo_sapiens.loom"
        self.download_website_meta = 'private'
        self.organ = "colon"
        self.sub_tissue = "lamina propria of mucosa of colon"
        self.annotated = True

        self.class_maps = {
            "0": {
                "Endothelial 1": "Endothelial",
                "Endothelial 2": "Endothelial",
                "Glial": "Glial cells",
                "Myofibroblasts": "Myofibroblasts",
                "Pericyte 1": "Pericytes",
                "Pericyte 2": "Pericytes",
                "Pericytes": "Pericytes",
                "Plasma Cells": "Plasma Cells",
                "Smooth Muscle": "Smooth Muscle",
                "Stromal 1": "Stromal",
                "Stromal 2a": "Stromal",
                "Stromal 2b": "Stromal",
                "Stromal 3": "Stromal",
                "Stromal 4": "Stromal",
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw:
            if fn is None:
                fn = [
                    os.path.join(self.path, "human", "colon", "f8aa201c-4ff1-45a4-890e-840d63459ca2.homo_sapiens.loom"),
                    os.path.join(self.path, "human", "colon", "uc_meta_data_stromal_with_donor.txt"),
                    os.path.join(self.path, "human", "colon", "hc_meta_data_stromal_with_donor.txt")
                ]
            adata = anndata.read_loom(fn[0])
            ctuc = pd.read_csv(fn[1], sep='\t')
            cthealthy = pd.read_csv(fn[2], sep='\t')
            adata = adata[adata.obs['emptydrops_is_cell'] == 't'].copy()
            adata = adata[adata.X.sum(axis=1).flatten() >= 250].copy()
            uc = adata[adata.obs['donor_organism.diseases.ontology_label'] == "ulcerative colitis (disease)"].copy()
            bcuc = [i.split('-')[0] for i in ctuc['Barcode']]
            seluc = []
            for i in uc.obs['barcode']:
                seluc.append((uc.obs['barcode'].str.count(i).sum() == 1) and i in bcuc)
            uc = uc[seluc].copy()
            ctuc.index = [i.split('-')[0] for i in ctuc['Barcode']]
            uc.obs['celltype'] = [ctuc.loc[i]['Cluster'] for i in uc.obs['barcode']]
            uc.var = uc.var.reset_index().rename(columns={'index': 'names'}).set_index('featurekey')
            healthy = adata[adata.obs['donor_organism.diseases.ontology_label'] == "normal"].copy()
            bchealthy = [i.split('-')[0] for i in cthealthy['Barcode']]
            selhealthy = []
            for i in healthy.obs['barcode']:
                selhealthy.append((healthy.obs['barcode'].str.count(i).sum() == 1) and i in bchealthy)
            healthy = healthy[selhealthy].copy()
            cthealthy.index = [i.split('-')[0] for i in cthealthy['Barcode']]
            healthy.obs['celltype'] = [cthealthy.loc[i]['Cluster'] for i in healthy.obs['barcode']]
            healthy.var = healthy.var.reset_index().rename(columns={'index': 'names'}).set_index('featurekey')
            self.adata = healthy.concatenate(uc)

        else:
            if fn is None:
                fn = os.path.join(self.path, "human", "colon", "kinchenetal.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[ADATA_IDS_SFAIRA.author] = 'Simmons'
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = "10.1016/j.cell.2018.08.067"
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs['celltype']
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = [line == 'normal' for line in
                                     self.adata.obs['donor_organism.diseases.ontology_label']]
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = self.adata.obs['donor_organism.diseases.ontology_label'].astype('category')
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = self.adata.obs[ADATA_IDS_SFAIRA.state_exact]\
            .cat.rename_categories({'normal': 'healthy', 'ulcerative colitis (disease)': 'ulcerative colitis'})

        self._convert_and_set_var_names(symbol_col=ADATA_IDS_SFAIRA.gene_id_names, ensembl_col='Accession')
