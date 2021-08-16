import anndata
import os
import pandas as pd

from sfaira.data import DatasetBase


class Dataset(DatasetBase):

    """
    TODO extra meta data in obs2

    age: columns "Age" contains integer entries and Unknown
    diseases: column "Diagnosis" contains entries NSIP, cHP, Control, IPF, ILD, Sarcoidosis
        column Tobacco contains entries Y,N
    ethnicity: column "Ethnicity" contains entries African_American, Caucasian, Hispanic, Unknown
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.download_url_data = [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5Fmatrix%2Emtx%2Egz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5Fgenes%2Etsv%2Egz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5Fbarcodes%2Etsv%2Egz"
        ]
        self.download_url_meta = [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl/GSE135893%5FIPF%5Fmetadata%2Ecsv%2Egz",
            "https://advances.sciencemag.org/highwire/filestream/234522/field_highwire_adjunct_files/2/aba1972_Table_S2.csv",
        ]

        self.author = "Habermann"
        self.doi_journal = "10.1126/sciadv.aba1972"
        self.doi_preprint = "10.1101/753806"
        self.normalization = "raw"
        self.organ = "lung parenchyma"
        self.organism = "human"
        self.assay_sc_obs_key = "Chemistry"
        self.year = 2020
        self.sample_source = "primary_tissue"
        self.sex_obs_key = "Gender"
        self.tech_sample_obs_key = "Sample_Name"

        self.gene_id_symbols_var_key = "index"

        self.cell_type_obs_key = "celltype"
        self.state_exact_obs_key = "Diagnosis"

        self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "GSE135893_matrix.mtx.gz"),
        os.path.join(data_dir, "GSE135893_genes.tsv.gz"),
        os.path.join(data_dir, "GSE135893_barcodes.tsv.gz"),
        os.path.join(data_dir, "GSE135893_IPF_metadata.csv.gz"),
        os.path.join(data_dir, "aba1972_Table_S2.csv"),
    ]
    adata = anndata.read_mtx(fn[0]).T
    adata.var = pd.read_csv(fn[1], index_col=0, header=None, names=["ids"])
    adata.obs = pd.read_csv(fn[2], index_col=0, header=None, names=["barcodes"])
    obs = pd.read_csv(fn[3], index_col=0)
    obs2 = pd.read_csv(fn[4], index_col=0)
    obs["Chemistry"] = [{"3_prime_V2": "10x 3' v2", "5_prime": "10x 5' v1"}[obs2.loc[x, "Chemistry"]]
                        for x in obs["orig.ident"].values]
    obs["Gender"] = [{"F": "female", "M": "male", "Unknown": "unknown"}[obs2.loc[x, "Gender"]]
                     for x in obs["orig.ident"].values]
    adata = adata[obs.index.tolist(), :].copy()
    adata.obs = obs

    return adata
