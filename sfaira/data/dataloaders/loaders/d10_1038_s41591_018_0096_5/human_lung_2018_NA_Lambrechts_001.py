import anndata
import tarfile
import scipy.io
import pandas as pd


def load(data_dir, **kwargs):
    fn = [
        data_dir + "387ZM2s",
        data_dir + "3dGyjGr",
    ]  # file names (these are tar.gz files without extension)
    # load counts, gene names and barcode names
    with tarfile.open(fn[0], mode="r:gz") as archive:  # mode for tar.gz files
        counts = scipy.io.mmread(
            archive.extractfile("./export/LC_counts/matrix.mtx")
        )  # path relative to folder in which .tar.gz file is located (so path is after extraction)
        gene_names = pd.read_csv(
            archive.extractfile("./export/LC_counts/genes.tsv"), sep="\t", header=None
        )
        barcodes = pd.read_csv(
            archive.extractfile("./export/LC_counts/barcodes.tsv"),
            sep="\t",
            header=None,
        )
    adata = anndata.AnnData(counts.T)
    adata.obs_names = barcodes[0]
    adata.var_names = gene_names[0]
    # import metadata
    meta = pd.read_csv(fn[1], compression="gzip", sep=",", index_col=0)
    # add metadata
    adata.obs = meta.loc[adata.obs.index, :]
    # extract sample names:
    adata.obs["sample"] = [bc.split("_")[0] for bc in adata.obs.index]
    # add single cell assay info (based on paper):
    patient_to_assay_dict = dict()
    for patient in [1, 2]:
        patient_to_assay_dict[patient] = "10X 3' v1 sequencing"
    for patient in range(3, 9):
        patient_to_assay_dict[patient] = "10X 3' v2 sequencing"
    adata.obs['sc_platform'] = adata.obs.PatientNumber.map(patient_to_assay_dict)
    return adata
