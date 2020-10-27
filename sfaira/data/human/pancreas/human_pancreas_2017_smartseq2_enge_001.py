import os
from typing import Union
from .external import DatasetBase
from .external import ADATA_IDS_SFAIRA
import anndata
import tarfile
import gzip
from io import StringIO
import anndata as ad
import pandas as pd
import scipy.sparse


class Dataset(DatasetBase):
    """
    This data loader supports reading of the downloaded raw data file if `load_raw=True` is passed to self.load()
    To download the datafile required by this dataloader, use the link provided as the `download_website` and
    `download_website_meta` attributes of this class. For (up to 100-fold faster) repeated data loading, please pass
    `load_raw=False` when calling the self.load() method. For this, you need to preprocess the raw files as below and
    place the resulting h5ad file in the data folder of this organ:

    import tarfile
    import os
    import gzip
    from io import StringIO
    import anndata as ad
    import pandas as pd
    import scipy.sparse
    dfs = []
    with tarfile.open("GSE81547_RAW.tar") as tar:
        for member in tar.getmembers():
            d = pd.read_csv(tar.extractfile(member), compression='gzip', header=None, sep='\t', index_col=0, names=[member.name.split("_")[0]])
            dfs.append(d)
    adata = ad.AnnData(pd.concat(dfs, axis=1).iloc[1:-6].T)
    adata.X = scipy.sparse.csc_matrix(adata.X)
    with gzip.open('GSE81547_series_matrix.txt.gz') as f:
        file_content = [i.decode("utf-8") for i in f.readlines()]
    inputstring = ''
    for line in file_content:
        if '"ID_REF"' in line:
            inputstring += line
        if '!Sample_title' in line:
            inputstring += line[1:]
        if '!Sample_characteristics_ch1\t"inferred_cell_type: alpha' in line:
            inputstring += line[1:]
    data = StringIO(inputstring)
    d = pd.read_csv(data, sep='\t').T
    d.columns=d.iloc[0]
    d.drop('Sample_title', inplace=True)
    d = d.reset_index().set_index('ID_REF')
    d.columns.name = None
    d.index.name = None
    adata.obs['celltype'] = [d.loc[i]['Sample_characteristics_ch1'].split(": ")[1] for i in adata.obs.index]
    adata.obs['patient'] = ["_".join(d.loc[i]['index'].split('_')[:2]) for i in adata.obs.index]
    adata.write('GSE81547.h5ad')

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
        self.id = "human_pancreas_2017_smartseq2_enge_001_10.1016/j.cell.2017.09.004"
        self.download_website = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81547/suppl/GSE81547_RAW.tar"
        self.download_website_meta = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81547/matrix/GSE81547_series_matrix.txt.gz"
        self.organ = "pancreas"
        self.sub_tissue = "islet of Langerhans"
        self.has_celltypes = True

        self.class_maps = {
            "0": {
                'alpha': 'Alpha cell',
                'acinar': 'Acinar cell',
                'ductal': 'Ductal cell',
                'beta': 'Beta cell',
                'unsure': 'Unknown',
                'delta': 'Delta cell',
                'mesenchymal': 'Mesenchymal Cell'
            },
        }

    def _load(self, fn=None):
        if fn is None and self.path is None:
            raise ValueError("provide either fn in load or path in constructor")

        if self._load_raw:
            if fn is None:
                fn = [
                    os.path.join(self.path, "human/pancreas/GSE81547_RAW.tar"),
                    os.path.join(self.path, "human/pancreas/GSE81547_series_matrix.txt.gz")
                ]
            dfs = []
            with tarfile.open(fn[0]) as tar:
                for member in tar.getmembers():
                    d = pd.read_csv(tar.extractfile(member), compression='gzip', header=None, sep='\t', index_col=0,
                                    names=[member.name.split("_")[0]])
                    dfs.append(d)
            self.adata = ad.AnnData(pd.concat(dfs, axis=1).iloc[1:-6].T)
            self.adata.X = scipy.sparse.csc_matrix(self.adata.X)
            with gzip.open(fn[1]) as f:
                file_content = [i.decode("utf-8") for i in f.readlines()]
            inputstring = ''
            for line in file_content:
                if '"ID_REF"' in line:
                    inputstring += line
                if '!Sample_title' in line:
                    inputstring += line[1:]
                if '!Sample_characteristics_ch1\t"inferred_cell_type: alpha' in line:
                    inputstring += line[1:]
            data = StringIO(inputstring)
            d = pd.read_csv(data, sep='\t').T
            d.columns = d.iloc[0]
            d.drop('Sample_title', inplace=True)
            d = d.reset_index().set_index('ID_REF')
            d.columns.name = None
            d.index.name = None
            self.adata.obs['celltype'] = [d.loc[i]['Sample_characteristics_ch1'].split(": ")[1] for i in self.adata.obs.index]
            self.adata.obs['patient'] = ["_".join(d.loc[i]['index'].split('_')[:2]) for i in self.adata.obs.index]

        else:
            if fn is None:
                fn = os.path.join(self.path, "human/pancreas/GSE81547.h5ad")
            self.adata = anndata.read(fn)

        self.adata.uns[ADATA_IDS.author] = "Quake"
        self.adata.uns[ADATA_IDS.year] = 2017
        self.adata.uns[ADATA_IDS.doi] = "10.1016/j.cell.2017.09.004"
        self.adata.uns[ADATA_IDS.protocol] = 'Smartseq2'
        self.adata.uns[ADATA_IDS.organ] = self.organ
        self.adata.uns[ADATA_IDS.subtissue] = self.sub_tissue
        self.adata.uns[ADATA_IDS.animal] = "human"
        self.adata.uns[ADATA_IDS.id] = self.id
        self.adata.uns[ADATA_IDS.wget_download] = self.download_website
        self.adata.uns[ADATA_IDS.has_celltypes] = self.has_celltypes
        self.adata.uns[ADATA_IDS.normalization] = 'raw'

        self.adata.obs[ADATA_IDS.healthy] = True
        self.adata.obs[ADATA_IDS.state_exact] = "healthy"

        self.adata.obs[ADATA_IDS.cell_ontology_class] = self.adata.obs['celltype']
        self._convert_and_set_var_names(symbol_col='index', ensembl_col=None, new_index=ADATA_IDS.gene_id_ensembl)
