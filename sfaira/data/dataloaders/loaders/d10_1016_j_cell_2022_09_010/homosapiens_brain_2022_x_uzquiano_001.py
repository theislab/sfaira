import anndata as ad
import os
import scipy.io
import pandas as pd
import numpy as np


def load(data_dir, sample_fn, **kwargs):
    p_ex = os.path.join(data_dir, "SCP1756", "expression")
    p_cl = os.path.join(data_dir, "SCP1756", "cluster")
    p_me = os.path.join(data_dir, "SCP1756", "metadata", "meta_all.txt")

    loadpaths = {
        "RNA_23days": [
            os.path.join(p_ex, "expression_23days.txt"),
            os.path.join(p_ex, "62041b09771a5b0d8cd54adb", "NormExpression_23days_barcodes.txt"),
            os.path.join(p_ex, "62041b09771a5b0d8cd54adb", "NormExpression_23days_genes.txt"),
            os.path.join(p_cl, "umap_23days.txt")
        ],
        "RNA_1month": [
            os.path.join(p_ex, "expression_1month.txt"),
            os.path.join(p_ex, "62041b09771a5b0d8cd54ad8", "NormExpression_1month_barcodes.txt"),
            os.path.join(p_ex, "62041b09771a5b0d8cd54ad8", "NormExpression_1month_genes.txt"),
            os.path.join(p_cl, "umap_1month.txt")
        ],
        "RNA_1.5month": [
            os.path.join(p_ex, "620176f7771a5b0dcfc1a0c8", "expression_1.5month.txt"),
            os.path.join(p_ex, "620176f7771a5b0dcfc1a0c8", "expression_1.5month_barcodes.txt"),
            os.path.join(p_ex, "62041b09771a5b0d8cd54ad3", "NormExpression_1.5month_genes.txt"),
            os.path.join(p_cl, "umap_1.5month.txt")
        ],
        "RNA_2month": [
            os.path.join(p_ex, "expression_2month.txt"),
            os.path.join(p_ex, "62042e17771a5b0db1d548c8", "NormExpression_2month_barcodes.txt"),
            os.path.join(p_ex, "62042e17771a5b0db1d548c8", "NormExpression_2month_genes.txt"),
            os.path.join(p_cl, "umap_2month.txt")
        ],
        "RNA_3month": [
            os.path.join(p_ex, "expression_3month.txt"),
            os.path.join(p_ex, "62041b09771a5b0d8cd54ae0", "NormExpression_3month_barcodes.txt"),
            os.path.join(p_ex, "62041b09771a5b0d8cd54ae0", "NormExpression_3month_genes.txt"),
            os.path.join(p_cl, "umap_3month.txt")
        ],
        "RNA_4month": [
            os.path.join(p_ex, "620178db771a5b0dcfc1a0d1", "expression_4month.txt"),
            os.path.join(p_ex, "620178db771a5b0dcfc1a0d1", "NormExpression_4month_barcodes.txt"),
            os.path.join(p_ex, "62042e17771a5b0db1d548d0", "NormExpression_4month_genes.txt"),
            os.path.join(p_cl, "umap_4month.txt")
        ],
        "RNA_5month": [
            os.path.join(p_ex, "62017a6e771a5b0dcfc1a0da", "expression_5month.txt"),
            os.path.join(p_ex, "62017a6e771a5b0dcfc1a0da", "NormExpression_5month_barcodes.txt"),
            os.path.join(p_ex, "62017a6e771a5b0dcfc1a0da", "NormExpression_5month_genes.txt"),
            os.path.join(p_cl, "umap_5month.txt")
        ],
        "RNA_6month": [
            os.path.join(p_ex, "62042e18771a5b0db1d548ed", "expression_6month.txt"),
            os.path.join(p_ex, "62042e18771a5b0db1d548ed", "expression_6month_barcodes.txt"),
            os.path.join(p_ex, "62042e18771a5b0db1d548ed", "expression_6month_genes.txt"),
            os.path.join(p_cl, "umap_6month.txt")
        ],
        "ATAC_1month": [
            os.path.join(p_ex, "expression_1month-ATAC.txt"),
            "from_meta_atac_1",
            os.path.join(p_ex, "62042e18771a5b0db1d548ea", "expression_6month-ATAC_genes.txt"),
            os.path.join(p_cl, "umap_1month-ATAC.txt")
        ],
        "ATAC_3month": [
            os.path.join(p_ex, "expression_3month-ATAC.txt"),
            "from_meta_atac_3",
            os.path.join(p_ex, "62042e18771a5b0db1d548ea", "expression_6month-ATAC_genes.txt"),
            os.path.join(p_cl, "umap_3month-ATAC.txt")
        ],
        "ATAC_6month": [
            os.path.join(p_ex, "62042e18771a5b0db1d548ea", "expression_6month-ATAC.txt"),
            os.path.join(p_ex, "62042e18771a5b0db1d548ea", "expression_6month-ATAC_barcodes.txt"),
            os.path.join(p_ex, "62042e18771a5b0db1d548ea", "expression_6month-ATAC_genes.txt"),
            os.path.join(p_cl, "umap_6month-ATAC.txt")
        ],
        "SHARE_RNA": [
            os.path.join(p_ex, "62019734771a5b0dcfc1a1de", "expression_SHARE-RNA.txt"),
            os.path.join(p_ex, "62019734771a5b0dcfc1a1de", "NormExpression_SHARE-RNA_barcodes.txt"),
            os.path.join(p_ex, "62042e18771a5b0db1d548e7", "NormExpression_SHARE-RNA_genes.txt"),
            os.path.join(p_cl, "umap_SHARE-RNA.txt")
        ],
        "SHARE_ATAC": [
            os.path.join(p_ex, "expression_SHARE-ATAC.txt"),
            os.path.join(p_ex, "62042e18771a5b0db1d548e4", "NormExpression_SHARE-ATAC_barcodes.txt"),
            os.path.join(p_ex, "62042e18771a5b0db1d548e4", "NormExpression_SHARE-ATAC_genes.txt"),
            os.path.join(p_cl, "umap_SHARE-ATAC.txt")
        ],
        "QUAD-GM_3.5month": [
            os.path.join(p_ex, "expression_Quad-GM.txt"),
            os.path.join(p_ex, "62042e18771a5b0db1d548de", "NormExpression_Quad-GM_barcodes.txt"),
            os.path.join(p_ex, "62042e18771a5b0db1d548de", "NormExpression_Quad-GM_genes.txt"),
            os.path.join(p_cl, "umap_Quad-GM.txt")
        ],
        "QUAD-H66_3.5month": [
            os.path.join(p_ex, "6201a53b771a5b0dcfc1a200", "expression_Quad-H66.txt"),
            os.path.join(p_ex, "6201a53b771a5b0dcfc1a200", "NormExpression_Quad-H66_barcodes.txt"),
            os.path.join(p_ex, "6201a53b771a5b0dcfc1a200", "NormExpression_Quad-H66_genes.txt"),
            os.path.join(p_cl, "umap_Quad-H66.txt")
        ],
        "FETAL": [
            os.path.join(p_ex, "62bf51ba57b44b9aa42292c4", "expression_HumanFetal.txt"),
            os.path.join(p_ex, "62bf51ba57b44b9aa42292c4", "expression_HumanFetal_barcodes.txt"),
            os.path.join(p_ex, "62bf51ba57b44b9aa42292c4", "expression_HumanFetal_genes.txt"),
            os.path.join(p_cl, "umap_HumanFetal.txt")
        ],
    }

    age_dict = {
        "23days": "23",
        "1month": "30",
        "1.5month": "45",
        "2month": "60",
        "3month": "90",
        "3.5month": "105",
        "4month": "120",
        "5month": "150",
        "6month": "180",
    }

    index_dict = {
        "from_meta_atac_1": "1m-ATAC",
        "from_meta_atac_3": "3m-ATAC",
    }

    meta = pd.read_csv(p_me, sep="\t", index_col=0, low_memory=False).tail(-1)
    x = scipy.io.mmread(loadpaths[sample_fn][0]).T.tocsr().astype(np.float32)
    var = pd.read_csv(loadpaths[sample_fn][2], sep="\t", index_col=0, low_memory=False, header=None, names=[None])
    if loadpaths[sample_fn][1].startswith("from_meta"):  # for a few samples there is no barcode information, so we extract it from the meta file
        obs = meta.loc[[index_dict[loadpaths[sample_fn][1]] in i for i in meta.index]]
    else:  # in all other cases, we just load the barcodes file
        obs = pd.read_csv(loadpaths[sample_fn][1], sep="\t", index_col=0, low_memory=False, header=None, names=[None])

    clusters = pd.read_csv(loadpaths[sample_fn][3], sep="\t", index_col=0, low_memory=False).tail(-1)
    clusters = clusters.loc[obs.index]

    meta = meta.loc[obs.index]
    if "SHARE" in sample_fn:
        meta["age"] = meta["biosample_id"].str.split("-").str[1]
        meta["organoid_age_days"] = meta["biosample_id"].str.split("-").str[1].replace({"d23": 23, "1m": 30, "2m": 60, "3m": 90}).values
    elif sample_fn == "FETAL":
        meta["age"] = meta["biosample_id"].str.split("-").str[2].values
        meta["organoid_age_days"] = meta["biosample_id"].str.split("-").str[2].str[3:].values.astype(int) * 7
    else:
        meta["age"] = sample_fn.split("_")[1]
        meta["organoid_age_days"] = age_dict[sample_fn.split("_")[-1]]
    meta["organoid_age_days"] = meta["organoid_age_days"] .astype(str)
    meta = pd.concat((meta, clusters[["CellType", "Cluster"]]), axis=1)

    adata = ad.AnnData(
        X=x,
        obs=meta,
        var=var,
        obsm={"X_umap": clusters[["X", "Y"]].values}
    )

    return adata
