import anndata as ad
import os
import scipy.io
import pandas as pd
import numpy as np


def load(data_dir, sample_fn, **kwargs):
    p_ex = os.path.join(data_dir, "SCP1129", "expression")
    p_cl = os.path.join(data_dir, "SCP1129", "cluster")
    p_me = os.path.join(data_dir, "SCP1129", "metadata", "meta_all.txt")

    loadpaths = {
        "GM_wt_harmonized": [
            os.path.join(p_ex, "5f5fd2ef771a5b6c993955bd", "gene_sorted-expression_GM_wt_harmonized.txt"),
            os.path.join(p_ex, "5f5fd2ef771a5b6c993955bd", "expression_GM_wt_harmonized_barcodes.txt"),
            os.path.join(p_ex, "5f5fd2ef771a5b6c993955bd", "expression_GM_wt_harmonized_genes.txt"),
            os.path.join(p_cl, "tsne_GM_wt_harmonized.txt")
        ],
        "Mito_wt_harmonized": [
            os.path.join(p_ex, "6101c9ad771a5b6c29d85cfe", "gene_sorted-expression_Mito_wt_harmonized.txt"),
            os.path.join(p_ex, "6101c9ad771a5b6c29d85cfe", "expression_Mito_wt_harmonized_barcodes.txt"),
            os.path.join(p_ex, "6101c9ad771a5b6c29d85cfe", "expression_Mito_wt_harmonized_genes.txt"),
            os.path.join(p_cl, "tsne_Mito_wt_harmonized.txt")
        ],
        "CHD8_GM_3mr1": [
            os.path.join(p_ex, "6101c977771a5b6c06d85df9", "gene_sorted-expression_CHD8_GM_3mr1.txt"),
            "from_tsne",
            # using genes files from different sample as this one is not available and the're all the same
            os.path.join(p_ex, "6101bf1e771a5b6c06d85dc0", "expression_CHD8_GM_3mr2_genes.txt"),
            os.path.join(p_cl, "tsne_CHD8_GM_3mr1.txt")
        ],
        "CHD8_GM_3mr2": [
            os.path.join(p_ex, "6101bf1e771a5b6c06d85dc0", "gene_sorted-expression_CHD8_GM_3mr2.txt"),
            os.path.join(p_ex, "6101bf1e771a5b6c06d85dc0", "expression_CHD8_GM_3mr2_barcodes.txt"),
            os.path.join(p_ex, "6101bf1e771a5b6c06d85dc0", "expression_CHD8_GM_3mr2_genes.txt"),
            os.path.join(p_cl, "tsne_CHD8_GM_3mr2.txt")
        ],
        "CHD8_H1_3m": [
            os.path.join(p_ex, "6101c39f771a5b6c06d85dcf", "gene_sorted-expression_CHD8_H1_3m.txt"),
            os.path.join(p_ex, "6101c39f771a5b6c06d85dcf", "expression_CHD8_H1_3m_barcodes.txt"),
            os.path.join(p_ex, "6101c39f771a5b6c06d85dcf", "expression_CHD8_H1_3m_genes.txt"),
            os.path.join(p_cl, "tsne_CHD8_H1_3m.txt")
        ],
        "CHD8_HUES66_3mr1": [
            os.path.join(p_ex, "5f694c76771a5b0de0ce9ccf", "gene_sorted-expression_HUES66_3.5mon_r1.txt"),
            os.path.join(p_ex, "5f694c76771a5b0de0ce9ccf", "expression_CHD8_HUES66_3mr1_barcodes.txt"),
            os.path.join(p_ex, "5f694c76771a5b0de0ce9ccf", "expression_HUES66_3.5mon_r1_genes.txt"),
            os.path.join(p_cl, "tsne_CHD8_HUES66_3mr1.txt")
        ],
        "CHD8_HUES66_3mr2": [
            os.path.join(p_ex, "5f69533f771a5b0de0ce9ce2", "gene_sorted-expression_HUES66_3.5mon_r2.txt"),
            os.path.join(p_ex, "5f69533f771a5b0de0ce9ce2", "expression_CHD8_HUES66_3mr2_barcodes.txt"),
            os.path.join(p_ex, "5f69533f771a5b0de0ce9ce2", "expression_HUES66_3.5mon_r2_genes.txt"),
            os.path.join(p_cl, "tsne_CHD8_HUES66_3mr2.txt")
        ],
        "CHD8_HUES66_6m": [
            os.path.join(p_ex, "6101c407771a5b6c06d85dd8", "gene_sorted-expression_CHD8_HUES66_6m.txt"),
            os.path.join(p_ex, "6101c407771a5b6c06d85dd8", "expression_CHD8_HUES66_6m_barcodes.txt"),
            os.path.join(p_ex, "6101c407771a5b6c06d85dd8", "expression_CHD8_HUES66_6m_genes.txt"),
            os.path.join(p_cl, "tsne_CHD8_HUES66_6m.txt")
        ],
        "SUV_PGP1_1m": [
            os.path.join(p_ex, "6101ca7d771a5b6c06d85e10", "gene_sorted-expression_SUV_PGP1_1m.txt"),
            os.path.join(p_ex, "6101ca7d771a5b6c06d85e10", "expression_SUV_PGP1_1m_barcodes.txt"),
            os.path.join(p_ex, "6101ca7d771a5b6c06d85e10", "expression_SUV_PGP1_1m_genes.txt"),
            None
        ],
        "SUV_Mito210_1m28d": [
            os.path.join(p_ex, "5f6924ce771a5b0dbe65da39", "gene_sorted-expression_SUV_28div.txt"),
            os.path.join(p_ex, "5f6924ce771a5b0dbe65da39", "expression_SUV_Mito210_1m28d_barcodes.txt"),
            os.path.join(p_ex, "5f6924ce771a5b0dbe65da39", "expression_SUV_28div_genes.txt"),
            os.path.join(p_cl, "tsne_SUV_Mito210_1m28d.txt")
        ],
        "SUV_Mito210_1m35d": [
            os.path.join(p_ex, "5f693648771a5b0dbe65da60", "gene_sorted-expression_SUV_35div.txt"),
            os.path.join(p_ex, "5f693648771a5b0dbe65da60", "expression_SUV_Mito210_1m35d_barcodes.txt"),
            os.path.join(p_ex, "5f693648771a5b0dbe65da60", "expression_SUV_35div_genes.txt"),
            os.path.join(p_cl, "tsne_SUV_Mito210_1m35d.txt")
        ],
        "SUV_Mito210_3mr1": [
            os.path.join(p_ex, "5f6939b3771a5b0dbe65da6b", "gene_sorted-expression_SUV_3mon.txt"),
            os.path.join(p_ex, "5f6939b3771a5b0dbe65da6b", "expression_SUV_Mito210_3mr1_barcodes.txt"),
            os.path.join(p_ex, "5f6939b3771a5b0dbe65da6b", "expression_SUV_3mon_genes.txt"),
            os.path.join(p_cl, "tsne_SUV_Mito210_3mr1.txt")
        ],
        "SUV_Mito210_3mr2": [
            os.path.join(p_ex, "6101c49b771a5b6c29d85cfa", "gene_sorted-expression_SUV_mito210_3mr2.txt"),
            os.path.join(p_ex, "6101c49b771a5b6c29d85cfa", "expression_SUV_mito210_3mr2_barcodes.txt"),
            os.path.join(p_ex, "6101c49b771a5b6c29d85cfa", "expression_SUV_mito210_3mr2_genes.txt"),
            os.path.join(p_cl, "tsne_SUV_mito210_3mr2.txt")
        ],
        "SUV_Mito294_1m": [
            os.path.join(p_ex, "6101c9fb771a5b6c06d85e07", "gene_sorted-expression_SUV_Mito294_1m.txt"),
            os.path.join(p_ex, "6101c9fb771a5b6c06d85e07", "expression_SUV_Mito294_1m_barcodes.txt"),
            os.path.join(p_ex, "6101c9fb771a5b6c06d85e07", "expression_SUV_Mito294_1m_genes.txt"),
            os.path.join(p_cl, "tsne_SUV_Mito294_1m.txt")
        ],
        "ARID1B_Mito210_1mr1": [
            os.path.join(p_ex, "61019a03771a5b6c06d85cf3", "gene_sorted-expression_ARID1B_Mito210_1mr1.txt"),
            os.path.join(p_ex, "61019a03771a5b6c06d85cf3", "expression_ARID1B_Mito210_1mr1_barcodes.txt"),
            os.path.join(p_ex, "61019a03771a5b6c06d85cf3", "expression_ARID1B_Mito210_1mr1_genes.txt"),
            os.path.join(p_cl, "tsne_ARID1B_Mito210_1mr1.txt")
        ],
        "ARID1B_Mito210_1mr2": [
            os.path.join(p_ex, "6101bb4d771a5b6c06d85d9a", "gene_sorted-expression_ARID1B_Mito210_1mr2.txt"),
            os.path.join(p_ex, "6101bb4d771a5b6c06d85d9a", "expression_ARID1B_Mito210_1mr2_barcodes.txt"),
            os.path.join(p_ex, "6101bb4d771a5b6c06d85d9a", "expression_ARID1B_Mito210_1mr2_genes.txt"),
            os.path.join(p_cl, "tsne_ARID1B_Mito210_1mr2.txt")
        ],
        "ARID1B_Mito210_3m": [
            os.path.join(p_ex, "gene_sorted-expression_ARID1B_Mito210_3m.txt"),
            "from_tsne",
            # using genes files from different sample as this one is not available and the're all the same
            os.path.join(p_ex, "6101bb4d771a5b6c06d85d9a", "expression_ARID1B_Mito210_1mr2_genes.txt"),
            os.path.join(p_cl, "tsne_ARID1B_Mito210_3m.txt")
        ],
        "ARID1B_Mito294_1m": [
            os.path.join(p_ex, "6101bdda771a5b6c29d85cf2", "gene_sorted-expression_ARID1B_Mito294_1m.txt"),
            os.path.join(p_ex, "6101bdda771a5b6c29d85cf2", "expression_ARID1B_Mito294_1m_barcodes.txt"),
            os.path.join(p_ex, "6101bdda771a5b6c29d85cf2", "expression_ARID1B_Mito294_1m_genes.txt"),
            os.path.join(p_cl, "tsne_ARID1B_Mito294_1m.txt")
        ],
    }

    days_dict = {
        "Mito_wt_harmonized": "178",
        "CHD8_GM_3mr1": "108",
        "CHD8_GM_3mr2": "108",
        "CHD8_H1_3m": "105",
        "CHD8_HUES66_3mr1": "109",
        "CHD8_HUES66_3mr2": "107",
        "CHD8_HUES66_6m": "190",
        "SUV_PGP1_1m": "36",
        "SUV_Mito210_1m28d": "28",
        "SUV_Mito210_1m35d": "35",
        "SUV_Mito210_3mr1": "92",
        "SUV_Mito210_3mr2": "90",
        "SUV_Mito294_1m": "38",
        "ARID1B_Mito210_1mr1": "35",
        "ARID1B_Mito210_1mr2": "35",
        "ARID1B_Mito210_3m": "90",
        "ARID1B_Mito294_1m": "35",
    }

    x = scipy.io.mmread(loadpaths[sample_fn][0]).T.tocsr().astype(np.float32)
    var = pd.read_csv(loadpaths[sample_fn][2], sep="\t", index_col=0, low_memory=False, header=None, names=[None])
    if loadpaths[sample_fn][1] == "from_tsne":  # for a few samples there is no barcode information, so we extract it from the tsne file index
        obs = pd.read_csv(loadpaths[sample_fn][3], sep="\t", index_col=0, low_memory=False).tail(-1).loc[:, []]
        obs.index.name = None
    else:  # in all other cases, we just load the barcodes file
        obs = pd.read_csv(loadpaths[sample_fn][1], sep="\t", index_col=0, low_memory=False, header=None, names=[None])
    meta = pd.read_csv(p_me, sep="\t", index_col=0, low_memory=False).tail(-1).loc[obs.index]

    if loadpaths[sample_fn][3] is not None:
        clusters = pd.read_csv(loadpaths[sample_fn][3], sep="\t", index_col=0, low_memory=False).tail(-1)
        obsm = {"X_tsne": clusters.loc[obs.index, ["X", "Y"]].values.astype(np.float32)}
    else:
        obsm = {}

    # convert lognormalized counts back to raw counts
    x = np.round(np.expm1(x).multiply(meta["nUMI"].values.astype(np.float32)[:, None]).multiply(1e-6))
    adata = ad.AnnData(
        X=x,
        obs=meta,
        var=var,
        obsm=obsm
    )

    # exclude samples which were generated in other publications
    if sample_fn == "GM_wt_harmonized":
        adata = adata[adata.obs["dataset"].isin(["GM_1mon", "GM_3mon"])].copy()
    elif sample_fn == "Mito_wt_harmonized":
        adata = adata[adata.obs["biosample_id"].isin(["Mito210map_7", "Mito210map_8", "Mito210map_9"])].copy()
    elif sample_fn == "CHD8_HUES66_3mr1":
        adata = adata[adata.obs["Condition"] == "mut"].copy()

    # add age in days annotation
    if sample_fn == "GM_wt_harmonized":
        adata.obs["organoid_age_days"] = adata.obs["biosample_id"].replace({"GM_1mon": "32", "GM_3mon": "98"})
    else:
        adata.obs["organoid_age_days"] = days_dict[sample_fn]

    # convert "Condition" obs column to be more meaningful
    adata.obs["Condition"] = adata.obs["Condition"].replace({"wt": "Wildtype", "mut": f"{sample_fn.split('_')[0]}_mutant"})

    return adata
