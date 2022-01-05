import anndata
import os
import scipy.sparse
import pandas as pd
import gzip
import scipy.io

assay_sc_map = {
    "10X 3'": "10x 3' transcription profiling",
    "10X 5'": "10x 5' transcription profiling",
}

disease_map = {
    "positive": "COVID-19",
    "negative": "healthy",
}

organ_map = {
    "B cells sorted from frozen PBMC (MACS, STEMCELL 19054)": "blood",
    "CD19+ B cell sorted from fresh PBMC (FACS)": "blood",
    "CD19+ B cell sorted from fresh PBMC (MACS)": "blood",
    "CD3+ T cell and CD19+ B cell sorted from fresh PBMC (FACS)": "blood",
    "CD3+ T cell sorted from fresh PBMC (FACS)": "blood",
    "fresh BALF": "lung",
    "fresh PBMC": "blood",
    "fresh PFMC": "lung",
    "fresh Sputum": "lung",
    "frozen PBMC": "blood",
}

sex_map = {
    "F": "female",
    "M": "male",
    "unknown": "unknown",
    "nan": "unknown",
}


def load(data_dir, **kwargs):
    fn = [os.path.join(data_dir, "GSE158055_covid19_counts.mtx.gz"),
          os.path.join(data_dir, "GSE158055_covid19_barcodes.tsv.gz"),
          os.path.join(data_dir, "GSE158055_covid19_features.tsv.gz"),
          os.path.join(data_dir, "GSE158055_cell_annotation.csv.gz"),
          os.path.join(data_dir, "GSE158055_sample_metadata.xlsx")]

    with gzip.open(fn[0], 'rb') as mm:
        X = scipy.io.mmread(mm).T.tocsr()
    obs = pd.read_csv(fn[1], header=None, sep="\t", index_col=0)
    obs.index.name = None
    var = pd.read_csv(fn[2], header=None, sep="\t", names=['names'])
    var.index = var['names'].values

    adata = anndata.AnnData(X=X, obs=obs, var=var)

    annotation = pd.read_csv(fn[3])
    metadata = pd.read_excel(fn[4], skiprows=20, skipfooter=37)
    metadata.columns = metadata.columns.str.split('characteristics: ').str[-1]
    tmp = annotation.merge(metadata, left_on='sampleID', right_on='Sample name').set_index('cellName')

    adata.obs = adata.obs.join(tmp).astype(str)

    adata.obs["assay_sc"] = adata.obs["Single cell sequencing platform"].map(assay_sc_map)
    adata.obs["disease"] = adata.obs["SARS-CoV-2"].map(disease_map)
    adata.obs["organ"] = adata.obs["Sample type"].map(organ_map)
    adata.obs["Sex"] = adata.obs["Sex"].map(sex_map)

    return adata
