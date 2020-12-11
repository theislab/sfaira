import anndata
import os
from typing import Union
from .external import DatasetBase
import tarfile
import pandas as pd
import scipy.sparse


class Dataset(DatasetBase):
    """
    This data loader supports reading of the downloaded raw data files if `load_raw=True` is passed to self.load()
    To download the datafile required by this dataloader, use the link provided as the `download_website` attribute of
    this class. The required celltype annotations for the data were kindly provided to us by the authors of the paper.
    Please contact them directly to pbtain the required annotation files (donor1.annotation.txt and
    donor2.annotation.txt). For (up to 100-fold faster) repeated data loading, please pass `load_raw=False` when calling
    the self.load() method. For this, you need to preprocess the raw files as below and place the resulting h5ad file in
    the data folder of this organ:

    import anndata
    import tarfile
    import pandas as pd
    import scipy.sparse
    adatas = []
    with tarfile.open('GSE126030_RAW.tar') as tar:
        for member in tar.getmembers():
            df = pd.read_csv(tar.extractfile(member.name), compression='gzip', sep='\t')
            df.index = [i.split('.')[0] for i in df['Accession']]
            var = pd.concat([df.pop(x) for x in ['Gene', 'Accession']], 1)
            if df.columns[-1].startswith('Un'):
                df.drop(df.columns[len(df.columns)-1], axis=1, inplace=True)
            adata = anndata.AnnData(df.T)
            adata.var = var
            if "PP001" in member.name or "PP002" in member.name:
                adata.obs['donor'] = 'Donor1'
                adata.obs['organ'] = 'Lung'
            elif "PP003" in member.name or "PP004" in member.name:
                adata.obs['donor'] = 'Donor1'
                adata.obs['organ'] = 'Bone Marrow'
            elif "PP005" in member.name or "PP006" in member.name:
                adata.obs['donor'] = 'Donor1'
                adata.obs['organ'] = 'Lymph Node'
            elif "PP009" in member.name or "PP010" in member.name:
                adata.obs['donor'] = 'Donor2'
                adata.obs['organ'] = 'Lung'
            elif "PP011" in member.name or "PP012" in member.name:
                adata.obs['donor'] = 'Donor2'
                adata.obs['organ'] = 'Bone Marrow'
            elif "PP013" in member.name or "PP014" in member.name:
                adata.obs['donor'] = 'Donor2'
                adata.obs['organ'] = 'Lymph Node'
            else:
                continue
            adata.obs.index = member.name.split('_')[1].split('s')[0]+'nskept.'+adata.obs.index
            adatas.append(adata)
    adata = adatas[0].concatenate(adatas[1:], index_unique=None)
    adata.obs.drop('batch', axis=1, inplace=True)
    adata = adata[:,adata.X.sum(axis=0) > 0].copy()
    adata.obs['cell_ontology_class'] = 'Unknown'
    df1 = pd.read_csv('donor1.annotation.txt', sep='\t', index_col=0, header=None)
    df2 = pd.read_csv('donor2.annotation.txt', sep='\t', index_col=0, header=None)
    for i in df1.index:
        adata.obs['cell_ontology_class'].loc[i] = df1.loc[i][1]
    for i in df2.index:
        adata.obs['cell_ontology_class'].loc[i] = df2.loc[i][1]
    adata.X = scipy.sparse.csc_matrix(adata.X)
    adata.write('GSE126030.h5ad')

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
        self.id = "human_mixed_2019_10x_szabo_001_10.1038/s41467-019-12464-3"
        self.download = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126030/suppl/GSE126030_RAW.tar"
        self.download_meta = 'private'
        self.organ = "mixed"
        self.sub_tissue = "Bone Marrow, Lung, Lymph Node"
        self.annotated = True
        self.loaded = False

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw:
            if fn is None:
                fn = [
                    os.path.join(self.path, "human", "mixed", "GSE126030_RAW.tar"),
                    os.path.join(self.path, "human", "mixed", "donor1.annotation.txt"),
                    os.path.join(self.path, "human", "mixed", "donor2.annotation.txt"),
                ]
            adatas = []
            with tarfile.open(fn[0]) as tar:
                for member in tar.getmembers():
                    df = pd.read_csv(tar.extractfile(member.name), compression='gzip', sep='\t')
                    df.index = [i.split('.')[0] for i in df['Accession']]
                    var = pd.concat([df.pop(x) for x in ['Gene', 'Accession']], 1)
                    if df.columns[-1].startswith('Un'):
                        df.drop(df.columns[len(df.columns) - 1], axis=1, inplace=True)
                    self.adata = anndata.AnnData(df.T)
                    self.adata.var = var
                    if "PP001" in member.name or "PP002" in member.name:
                        self.adata.obs['donor'] = 'Donor1'
                        self.adata.obs['organ'] = 'Lung'
                    elif "PP003" in member.name or "PP004" in member.name:
                        self.adata.obs['donor'] = 'Donor1'
                        self.adata.obs['organ'] = 'Bone Marrow'
                    elif "PP005" in member.name or "PP006" in member.name:
                        self.adata.obs['donor'] = 'Donor1'
                        self.adata.obs['organ'] = 'Lymph Node'
                    elif "PP009" in member.name or "PP010" in member.name:
                        self.adata.obs['donor'] = 'Donor2'
                        self.adata.obs['organ'] = 'Lung'
                    elif "PP011" in member.name or "PP012" in member.name:
                        self.adata.obs['donor'] = 'Donor2'
                        self.adata.obs['organ'] = 'Bone Marrow'
                    elif "PP013" in member.name or "PP014" in member.name:
                        self.adata.obs['donor'] = 'Donor2'
                        self.adata.obs['organ'] = 'Lymph Node'
                    else:
                        continue
                    self.adata.obs.index = member.name.split('_')[1].split('s')[0] + 'nskept.' + self.adata.obs.index
                    adatas.append(self.adata)
            self.adata = adatas[0].concatenate(adatas[1:], index_unique=None)
            self.adata.obs.drop('batch', axis=1, inplace=True)
            self.adata = self.adata[:, self.adata.X.sum(axis=0) > 0].copy()
            self.adata.obs['cell_ontology_class'] = 'Unknown'
            df1 = pd.read_csv(fn[1], sep='\t', index_col=0, header=None)
            df2 = pd.read_csv(fn[2], sep='\t', index_col=0, header=None)
            for i in df1.index:
                self.adata.obs['cell_ontology_class'].loc[i] = df1.loc[i][1]
            for i in df2.index:
                self.adata.obs['cell_ontology_class'].loc[i] = df2.loc[i][1]
            self.adata.X = scipy.sparse.csc_matrix(self.adata.X)
        else:
            if fn is None:
                fn = os.path.join(self.path, "human", "mixed", "GSE126030.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[self._ADATA_IDS_SFAIRA.author] = "Sims"
        self.adata.uns[self._ADATA_IDS_SFAIRA.year] = 2019
        self.adata.uns[self._ADATA_IDS_SFAIRA.doi] = "10.1038/s41467-019-12464-3"
        self.adata.uns[self._ADATA_IDS_SFAIRA.protocol] = "10x"
        self.adata.uns[self._ADATA_IDS_SFAIRA.organ] = self.organ
        self.adata.uns[self._ADATA_IDS_SFAIRA.subtissue] = self.sub_tissue
        self.adata.uns[self._ADATA_IDS_SFAIRA.species] = "human"
        self.adata.uns[self._ADATA_IDS_SFAIRA.id] = self.id
        self.adata.uns[self._ADATA_IDS_SFAIRA.download] = self.download
        self.adata.uns[self._ADATA_IDS_SFAIRA.annotated] = self.annotated
        self.adata.uns[self._ADATA_IDS_SFAIRA.normalization] = 'raw'

        self.adata.obs[self._ADATA_IDS_SFAIRA.subtissue] = self.adata.obs["organ"]
        self.adata.obs[self._ADATA_IDS_SFAIRA.healthy] = True
        self.adata.obs[self._ADATA_IDS_SFAIRA.state_exact] = 'healthy'

        self._convert_and_set_var_names(symbol_col='Gene', ensembl_col='Accession')

        # If the subset_organs() method has been run before, subset to specified organs
        if "organsubset" in self.__dict__:
            self.adata = self.adata[self.adata.obs['organ'].isin(self.organsubset)]
        # If adata object is empty, set it to None
        if not len(self.adata):
            self.adata = None
        self.loaded = True

    @property
    def ncells(self):
        if "organsubset" in self.__dict__:
            if not self.loaded:
                self._load()
            if self.adata is None:
                return 0
            else:
                return self.adata.n_obs
        else:
            return super().ncells
