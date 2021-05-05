import os
import pandas
import shutil
import zipfile


def load(data_dir, **kwargs):
    age_dict = {
        17: "17th week post-fertilization human stage",
        18: "18th week post-fertilization human stage",
    }

    import anndata2ri
    from rpy2.robjects import r
    fn = os.path.join(data_dir, "sc_dev_cortex_geschwind.zip")
    fn_tmp = os.path.join(os.path.expanduser("~"), "sfaira_tmp")
    if not os.path.exists(fn_tmp):
        os.makedirs(fn_tmp)
    with zipfile.ZipFile(fn, 'r') as zip_ref:
        zip_ref.extractall(fn_tmp)
    anndata2ri.activate()  # TODO: remove global activation of anndata2ri and use localconverter once it's fixed
    adata = r(
        f"library(Seurat)\n"
        f"load('{os.path.join(fn_tmp, 'raw_counts_mat.rdata')}')\n"
        f"new_obj = CreateSeuratObject(raw_counts_mat)\n"
        f"as.SingleCellExperiment(new_obj)\n"
    )
    obs = pandas.read_csv(os.path.join(fn_tmp, "cell_metadata.csv"), index_col=0)
    adata = adata[obs.index.tolist()].copy()
    adata.obs = obs
    shutil.rmtree(fn_tmp)
    adata.obs['DevStage'] = [age_dict[i] for i in adata.obs['Gestation_week']]
    return adata
