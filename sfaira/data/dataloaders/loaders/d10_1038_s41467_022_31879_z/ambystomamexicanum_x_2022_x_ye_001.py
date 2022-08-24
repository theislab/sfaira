import anndata
import numpy as np
import os
import scipy.sparse

import pandas as pd

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    sample_dir = sample_fn.split(":")[0]
    sample_fn = sample_fn.split(":")[1]
    if sample_dir == "15732027":
        dir_counts = os.path.join(data_dir, sample_dir, sample_fn)
        dir_counts = buffered_decompress(dir_counts)
        # Fine names in folders (these are not directly inferable from zip name:
        fn_counts_raw = [x for x in os.listdir(dir_counts) if x.endswith("dge.csv")][0]
        fn_genes = [x for x in os.listdir(dir_counts) if x.endswith("dge_gene.csv")][0]
        fn_counts = os.path.join(dir_counts, fn_counts_raw)
        fn_genes = os.path.join(dir_counts, fn_genes)
        # Need to use name from inside of zip here:
        fn_meta = os.path.join(data_dir, sample_dir,
                               sample_fn.split("_")[0] + "_" + fn_counts_raw.replace("dge.csv", "annotation.csv"))
        tab_x = pd.read_csv(fn_counts, header=0, index_col=False).T
        tab_genes = pd.read_csv(fn_genes, header=0, index_col=0)
        genes = ["-".join(x.split("-")[2:]) for x in tab_genes["x"].values]
        if sample_fn == "Neotenic_cloaca_dge.zip":
            # Table is structured differently:
            tab_meta = pd.read_csv(fn_meta, header=0, index_col=False)
            tab_meta.columns = ["Cluster"]
        else:
            tab_meta = pd.read_csv(fn_meta, header=0, index_col=0)
        # Take out sub cluster information:
        tab_meta["cell_type"] = [x.split("_") for x in tab_meta["Cluster"].values]
    else:
        fn_counts = os.path.join(data_dir, sample_dir, sample_fn)
        fn_meta = os.path.join(data_dir, sample_dir, sample_fn.replace("dge.zip", "Annotation.csv"))
        tab_x = pd.read_csv(fn_counts, header=0, index_col=0).T
        genes = ["_".join(x.split("_")[2:]) for x in tab_x.columns]
        tab_meta = pd.read_csv(fn_meta, header=0, index_col=0)
        tab_meta.columns = ["cell_type"]
    if sample_fn != "Neotenic_cloaca_dge.zip":
        assert np.all(tab_meta.index == tab_meta["cell_type"].index)
    adata = anndata.AnnData(scipy.sparse.csr_matrix(tab_x.values), var=pd.DataFrame({}, index=genes), obs=tab_meta)
    return adata
