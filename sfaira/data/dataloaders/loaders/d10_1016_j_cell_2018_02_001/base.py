import anndata
import numpy as np
import pandas
from typing import Union
from sfaira.data import DatasetBase
import zipfile
import tarfile
import os


class Dataset_d10_1016_j_cell_2018_02_001(DatasetBase):
    """
    This is a dataloader template for mca data.
    """

    def __init__(
            self,
            path: Union[str, None],
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, **kwargs)

        self.download_url_data = "https://ndownloader.figshare.com/articles/5435866?private_link=865e694ad06d5857db4b"
        self.download_url_meta = None

        self.obs_key_cellontology_class = "Annotation"
        self.obs_key_cellontology_original = "Annotation"

        self.author = "Guo"
        self.doi = "10.1016/j.cell.2018.02.001"
        self.normalization = "raw"
        self.healthy = True
        self.organism = "mouse"
        self.protocol = "microwell-seq"
        self.state_exact = "healthy"
        self.year = 2018

        self.var_symbol_col = "index"

    def _load_generalized(self, fn, samplename):
        if fn is None:
            fn = os.path.join(self.path, "mouse", "temp_mouse_atlas")

        with zipfile.ZipFile(os.path.join(fn, '5435866.zip')) as archive:
            celltypes = pandas.read_csv(archive.open('MCA_CellAssignments.csv'), index_col=1)
            celltypes = celltypes.drop(["Unnamed: 0"], axis=1)

            with tarfile.open(fileobj=archive.open('MCA_500more_dge.tar.gz')) as tar:
                data = pandas.read_csv(tar.extractfile(f'500more_dge/{samplename}.txt.gz'),
                                       compression="gzip",
                                       sep=" ",
                                       header=0
                                       )

        self.adata = anndata.AnnData(data.T)
        self.adata = self.adata[np.array([x in celltypes.index for x in self.adata.obs_names])].copy()
        self.adata.obs = celltypes.loc[self.adata.obs_names, :]
