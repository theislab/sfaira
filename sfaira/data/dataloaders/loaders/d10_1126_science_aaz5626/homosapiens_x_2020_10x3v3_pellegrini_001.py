import anndata
import os
import scipy.sparse
import pandas as pd
import numpy as np


def load(data_dir, **kwargs):

    em = pd.read_csv(os.path.join(data_dir, "exprMatrix.tsv.gz"), sep="\t", index_col=0).T
    umap_coords = pd.read_csv(os.path.join(data_dir, "Seurat_umap.coords.tsv.gz"), sep="\t", index_col=0, header=None)

    X = scipy.sparse.csr_matrix(em.values, dtype=np.float32)
    var = pd.DataFrame(data={}, index=[i.split("|")[0] for i in em.columns])
    obs = pd.read_csv(os.path.join(data_dir, "meta.tsv"), sep="\t", index_col=0)
    obsm = {"X_umap": umap_coords.loc[obs.index].values}

    assay_diff_dict = {
        "Tel": "Lancaster, 2014 (doi: 10.1038/nprot.2014.158); Lancaster, 2017 (doi: 10.1038/nbt.3906)",
        "ChP1": "Lancaster, 2014 (doi: 10.1038/nprot.2014.158); Lancaster, 2017 (doi: 10.1038/nbt.3906);"
        "added CHIR and BMP4 for ChP patterning on day 10 (Pellegrini, 2020 (doi: 10.1126/science.aaz5626))",
        "ChP2": "Lancaster, 2014 (doi: 10.1038/nprot.2014.158); Lancaster, 2017 (doi: 10.1038/nbt.3906);"
        "added CHIR and BMP4 for ChP patterning on day 10 (Pellegrini, 2020 (doi: 10.1126/science.aaz5626))",
        "ChP3": "Lancaster, 2014 (doi: 10.1038/nprot.2014.158); Lancaster, 2017 (doi: 10.1038/nbt.3906);"
        "added CHIR and BMP4 for ChP patterning on day 10 (Pellegrini, 2020 (doi: 10.1126/science.aaz5626))"}

    assay_type_diff_dict = {
        "Tel": "unguided",
        "ChP1": "guided",
        "ChP2": "guided",
        "ChP3": "guided"}

    organ_dict = {
        "Tel": "telencephalon",
        "ChP1": "choroid plexus",
        "ChP2": "choroid plexus",
        "ChP3": "choroid plexus"}

    cell_line_dict = {
        "Tel": "WA09",
        "ChP1": "WA01",
        "ChP2": "WA01",
        "ChP3": "WA01"}

    organoid_age_dict = {
        "Tel": "55",
        "ChP1": "27",
        "ChP2": "46",
        "ChP3": "53"}

    obs["assay_diff"] = [assay_diff_dict[i] for i in obs["orig.ident"]]
    obs["assay_type_diff"] = [assay_type_diff_dict[i] for i in obs["orig.ident"]]
    obs["organ"] = [organ_dict[i] for i in obs["orig.ident"]]
    obs["cell_line"] = [cell_line_dict[i] for i in obs["orig.ident"]]
    obs["organoid_age_days"] = [organoid_age_dict[i] for i in obs["orig.ident"]]

    adata = anndata.AnnData(X=X, obs=obs, var=var, obsm=obsm)

    return adata
