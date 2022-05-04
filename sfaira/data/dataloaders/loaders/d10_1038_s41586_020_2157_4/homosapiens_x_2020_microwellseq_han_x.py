import anndata
import numpy as np
import os
import pandas as pd
import scipy.sparse
import zipfile


def load(data_dir, **kwargs):
    adata = anndata.read(os.path.join(data_dir, "HCL_Fig1_adata.h5ad"))
    # convert to sparse matrix
    adata.X = scipy.sparse.csr_matrix(adata.X).copy()

    # harmonise annotations
    for col in ["batch", "tissue"]:
        adata.obs[col] = adata.obs[col].astype("str")
    adata.obs.index = adata.obs.index.str.replace("AdultJeJunum", "AdultJejunum", regex=True).str.replace(
        "AdultGallBladder", "AdultGallbladder", regex=True).str.replace(
        "FetalFemaleGonald", "FetalFemaleGonad", regex=True)
    adata.obs.replace({"AdultJeJunum": "AdultJejunum", "AdultGallBladder": "AdultGallbladder",
                       "FetalFemaleGonald": "FetalFemaleGonad"}, regex=True, inplace=True)
    adata.obs.index = ["-".join(i.split("-")[:-1]) for i in adata.obs.index]

    # load celltype labels and harmonise them
    # This pandas code should work with pandas 1.2 but it does not and yields an empty data frame:
    fig1_anno = pd.read_excel(
        os.path.join(data_dir, "HCL_Fig1_cell_Info.xlsx"),
        index_col="cellnames",
        engine="xlrd",  # ToDo: Update when pandas xlsx reading with openpyxl is fixed: yields empty tables
    )
    fig1_anno.index = fig1_anno.index.str.replace("AdultJeJunum", "AdultJejunum", regex=True).str.replace(
        "AdultGallBladder", "AdultGallbladder", regex=True).str.replace(
        "FetalFemaleGonald", "FetalFemaleGonad", regex=True)

    # check that the order of cells and cell labels is the same
    assert np.all(fig1_anno.index == adata.obs.index)

    # add annotations to adata object and rename columns
    adata.obs = pd.concat([adata.obs, fig1_anno[["cluster", "stage", "donor", "celltype"]]], axis=1)
    adata.obs.columns = ["sample", "tissue", "n_genes", "n_counts", "cluster_global", "stage", "donor",
                         "celltype_global"]

    # add sample-wise annotations to the full adata object
    df = pd.DataFrame(
        columns=["Cell_barcode", "Sample", "Batch", "Cell_id", "Cluster_id", "Ages", "Development_stage", "Method",
                 "Gender", "Source", "Biomaterial", "Name", "ident", "Celltype"])
    archive = zipfile.ZipFile(os.path.join(data_dir, "annotation_rmbatch_data_revised417.zip"))
    for f in archive.namelist():
        df1 = pd.read_csv(archive.open(f), encoding="unicode_escape")
        df = pd.concat([df, df1], sort=True)
    df = df.set_index("Cell_id")
    adata = adata[[i in df.index for i in adata.obs.index]].copy()
    a_idx = adata.obs.index.copy()
    adata.obs = pd.concat([adata.obs, df[
        ["Ages", "Celltype", "Cluster_id", "Gender", "Method", "Source"]
    ]], axis=1)
    assert np.all(a_idx == adata.obs.index)

    # remove musmusculus cells from the object  # ToDo: add this back in as musmusculus data sets?
    adata = adata[adata.obs["Source"] != "MCA2.0"].copy()

    # tidy up the column names of the obs annotations
    adata.obs.columns = [
        "sample", "sub_tissue", "n_genes", "n_counts", "cluster_global", "dev_stage", "donor", "celltype_global",
        "age", "celltype_specific", "cluster_specific", "sex", "assay_sc", "source"]
    # Remove new line characters from cell type:
    adata.obs["celltype_specific"] = [
        x.replace("\n", "").rstrip()
        for x in adata.obs["celltype_specific"].values
    ]
    adata.obs["organ"] = adata.obs["sample"].values
    adata.obs["sex"] = adata.obs["sex"].values
    # TODO are the more exact developmental stages in dev_stage?
    adata.obs["dev_stage"] = adata.obs["sample"].values

    return adata
