import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
import anndata


class Dataset(DatasetBase):
    """
    The input file for this dataloader (fetal_liver_alladata_.h5ad) was kindly provided to us by the
    authors of the publication. Please contact them directly to obtain the required file.

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
        self.id = "human_liver_2019_10x_popescu_001_10.1038/s41586-019-1652-y"
        self.download_website = "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7407/"
        self.organ = "liver"
        self.sub_tissue = "liver"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'B cell': 'Mature B cells',
                'DC1': 'Dendritic cell 1',
                'DC2': 'Dendritic cell 2',
                'DC precursor': 'Dendritic cell precursor',
                'Early Erythroid': 'Early Erythroid',
                'Early lymphoid_T lymphocyte': 'Early lymphoid T lymphocyte',
                'Endothelial cell': 'Endothelial cell',
                'Fibroblast': 'Fibroblast',
                'HSC_MPP': 'HSC MPP',
                'Hepatocyte': 'Hepatocyte',
                'ILC precursor': 'ILC precursor',
                'Kupffer Cell': 'Kupffer Cell',
                'Late Erythroid': 'Late Erythroid',
                'MEMP': 'MEMP',
                'Mast cell': 'Mast cell',
                'Megakaryocyte': 'Megakaryocyte',
                'Mid Erythroid': 'Mid Erythroid',
                'Mono-Mac': 'Mono Macrophage',
                'Monocyte': 'Monocyte',
                'Monocyte precursor': 'Monocyte precursor',
                'NK': 'NK cell',
                'Neutrophil-myeloid progenitor': 'Neutrophil myeloid progenitor',
                'Pre pro B cell': 'Pre pro B cell',
                'VCAM1+ EI macrophage': 'VCAM1pos EI macrophage',
                'pDC precursor': 'pDendritic cell precursor',
                'pre-B cell': 'pre B cell',
                'pro-B cell': 'pro B cell'
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw or not self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human/liver/fetal_liver_alladata_.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[ADATA_IDS_SFAIRA.author] = 'Haniffa'
        self.adata.uns[ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[ADATA_IDS_SFAIRA.doi] = '10.1038/s41586-019-1652-y'
        self.adata.uns[ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS_SFAIRA.animal] = "human"
        self.adata.uns[ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[ADATA_IDS_SFAIRA.download] = self.download_website
        self.adata.uns[ADATA_IDS_SFAIRA.annotated] = self.has_celltypes
        self.adata.uns[ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[ADATA_IDS_SFAIRA.cell_ontology_class] = self.adata.obs["cell.labels"]
        self.adata.obs[ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS_SFAIRA.gene_id_ensembl)
