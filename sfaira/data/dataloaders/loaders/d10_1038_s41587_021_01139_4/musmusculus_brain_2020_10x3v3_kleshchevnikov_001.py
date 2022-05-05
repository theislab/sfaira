import os
import pandas as pd

import scanpy

from sfaira.data.dataloaders.utils import buffered_decompress


def load(data_dir, sample_fn, **kwargs):
    fn_archive = os.path.join(data_dir, "E-MTAB-11115.processed.5.zip")
    fn_meta_archive = os.path.join(data_dir, "E-MTAB-11115.processed.4.zip")

    fn = os.path.join(buffered_decompress(fn_archive), sample_fn)
    fn_meta = os.path.join(buffered_decompress(fn_meta_archive), "cell_annotation.csv")

    adata = scanpy.read_10x_h5(fn)
    sample = sample_fn.split("_")[0]
    meta = pd.read_csv(fn_meta)
    meta = meta.loc[meta["sample"].values == sample, :]
    meta.index = [x.split("_")[1] for x in meta["Cell ID"].values]
    adata.obs["cell_type"] = [meta.loc[x, "annotation_1"] if x in meta.index else "unknown" for x in adata.obs_names]

    return adata
