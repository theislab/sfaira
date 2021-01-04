import anndata
import os
from typing import Union

from sfaira.data import DatasetBase


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
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)
        self.organism = "loaders"
        self.id = "human_liver_2019_10x_popescu_001_10.1038/s41586-019-1652-y"
        self.download = "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7407/"
        self.download_meta = 'private'
        self.organ = "liver"
        self.sub_tissue = "liver"
        self.author = 'Haniffa'
        self.year = 2019
        self.doi = '10.1038/s41586-019-1652-y'
        self.protocol = '10x'
        self.normalization = 'raw'
        self.healthy = True
        self.state_exact = 'healthy'
        self.var_symbol_col = 'index'
        self.obs_key_cellontology_original = 'cell.labels'

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
        if fn is None:
            fn = os.path.join(self.path, "human", "liver", "fetal_liver_alladata_.h5ad")
        self.adata = anndata.read(fn)
