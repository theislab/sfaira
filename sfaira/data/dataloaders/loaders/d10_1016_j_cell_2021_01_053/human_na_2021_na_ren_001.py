import anndata
import os
import scipy.sparse

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


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, sample_fn)
    adata = anndata.read(fn)
    adata.X = scipy.sparse.csr_matrix(adata.X)

    adata.obs["assay_sc"] = adata.obs["Single cell sequencing platform"].map(assay_sc_map)
    adata.obs["disease"] = adata.obs["SARS-CoV-2"].map(disease_map)
    adata.obs["organ"] = adata.obs["Sample type"].map(organ_map)
    adata.obs["Sex"] = adata.obs["Sex"].map(sex_map)

    return adata
