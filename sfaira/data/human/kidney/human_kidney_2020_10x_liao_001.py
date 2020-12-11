import anndata
import os
from typing import Union
from .external import DatasetBase
import pandas as pd
import scipy.io
import gzip
import tarfile


class Dataset(DatasetBase):
    """
    This data loader supports reading of the downloaded raw data file if `load_raw=True` is passed to self.load()
    To download the datafile required by this dataloader, use the link provided as the `download_website` attribute of
    this class. For (up to 100-fold faster) repeated data loading, please pass `load_raw=False` when calling the
    self.load() method. For this, you need to preprocess the raw files as below and place the resulting h5ad file in the
    data folder of this organ:

    import anndata
    import pandas as pd
    import scipy.io
    import gzip
    import tarfile
    adatas = []
    with tarfile.open("GSE131685_RAW.tar") as tar:
        for member in tar.getmembers():
            if '_matrix.mtx.gz' in member.name:
                name = '_'.join(member.name.split('_')[:-1])
                with gzip.open(tar.extractfile(member), 'rb') as mm:
                    X = scipy.io.mmread(mm).T.tocsr()
                obs = pd.read_csv(tar.extractfile(name+'_barcodes.tsv.gz'), compression='gzip', header=None, sep='\t', index_col=0)
                obs.index.name = None
                var = pd.read_csv(tar.extractfile(name+'_features.tsv.gz'), compression='gzip', header=None, sep='\t').iloc[:,:2]
                var.columns = ['ensembl', 'names']
                var.index = var['ensembl'].values
                adata = anndata.AnnData(X=X, obs=obs, var=var)
                adata.obs['sample'] = name
                adatas.append(adata)
    adata = adatas[0].concatenate(adatas[1:])
    del adata.obs['batch']
    adata.write('GSE131685.h5ad')

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
        self.id = "human_kidney_2020_10x_liao_001_10.1038/s41597-019-0351-8"
        self.download = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131685/suppl/GSE131685_RAW.tar"
        self.download_meta = None
        self.organ = "kidney"
        self.sub_tissue = "kidney"
        self.annotated = False

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw:
            if fn is None:
                fn = os.path.join(self.path, "human", "kidney", "GSE131685_RAW.tar")
            adatas = []
            with tarfile.open(fn) as tar:
                for member in tar.getmembers():
                    if '_matrix.mtx.gz' in member.name:
                        name = '_'.join(member.name.split('_')[:-1])
                        with gzip.open(tar.extractfile(member), 'rb') as mm:
                            X = scipy.io.mmread(mm).T.tocsr()
                        obs = pd.read_csv(tar.extractfile(name + '_barcodes.tsv.gz'), compression='gzip', header=None,
                                          sep='\t', index_col=0)
                        obs.index.name = None
                        var = pd.read_csv(tar.extractfile(name + '_features.tsv.gz'), compression='gzip', header=None,
                                          sep='\t').iloc[:, :2]
                        var.columns = ['ensembl', 'names']
                        var.index = var['ensembl'].values
                        self.adata = anndata.AnnData(X=X, obs=obs, var=var)
                        self.adata.obs['sample'] = name
                        adatas.append(self.adata)
            self.adata = adatas[0].concatenate(adatas[1:])
            del self.adata.obs['batch']

        else:
            if fn is None:
                fn = os.path.join(self.path, "human", "kidney", "GSE131685.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = 'Mo'
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = 2020
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = '10.1038/s41597-019-0351-8'
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = '10x'
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = self.species
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download
        self.adata.uns[self._ADATA_IDS_SFAIRA.download_meta] = self.download_meta
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[self._ADATA_IDS_SFAIRA.cell_ontology_class] = None
        self.adata.obs[self._ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[self._ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col="names", ensembl_col="ensembl")
