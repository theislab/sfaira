import numpy as np
import os
import pytest

from sfaira.data import DatasetSuperGroup, DistributedStore
from sfaira.data import Universe

MOUSE_GENOME_ANNOTATION = "Mus_musculus.GRCm38.102"

dir_data = "../test_data"
dir_meta = "../test_data/meta"


"""
TODO tests from here on down require cached data for mouse lung
"""


def test_dsg_store_config():
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.streamline_features(remove_gene_version=True, match_to_reference={"mouse": MOUSE_GENOME_ANNOTATION},
                           subset_genes_to_type="protein_coding")
    ds.streamline_metadata(schema="sfaira", uns_to_obs=False, clean_obs=True, clean_var=True, clean_uns=True,
                           clean_obs_names=True)
    store_path = os.path.join(dir_data, "store")
    config_path = os.path.join(store_path, "lung")
    ds.write_distributed_store(dir_cache=store_path, store="h5ad", dense=True)
    store = DistributedStore(cache_path=store_path)
    store.subset(attr_key="organ", values=["lung"])
    store.write_config(fn=config_path)
    store2 = DistributedStore(cache_path=store_path)
    store2.load_config(fn=config_path)
    assert np.all(store.indices == store2.indices)
