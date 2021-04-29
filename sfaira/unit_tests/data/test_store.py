import numpy as np
import os
import pytest

from sfaira.data import DistributedStore
from sfaira.data import Universe

MOUSE_GENOME_ANNOTATION = "Mus_musculus.GRCm38.102"

dir_data = "../test_data"
dir_meta = "../test_data/meta"


"""
TODO tests from here on down require cached data for mouse lung
"""


def test_store_config():
    """
    Test that data set config files can be set, written and recovered.
    """
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
    store.subset(attr_key="assay_sc", values=["10x sequencing"])
    store.subset_cells(attr_key="assay_sc", values=["10x sequencing"])
    store.write_config(fn=config_path)
    store2 = DistributedStore(cache_path=store_path)
    store2.load_config(fn=config_path)
    assert np.all(store.indices.keys() == store2.indices.keys())
    assert np.all([np.all(store.indices[k] == store2.indices[k]) for k in store.indices.keys()])


def test_store_type_targets():
    """
    Test that target leave nodes can be set, written and recovered.
    """
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.streamline_features(remove_gene_version=True, match_to_reference={"mouse": MOUSE_GENOME_ANNOTATION},
                           subset_genes_to_type="protein_coding")
    ds.streamline_metadata(schema="sfaira", uns_to_obs=False, clean_obs=True, clean_var=True, clean_uns=True,
                           clean_obs_names=True)
    store_path = os.path.join(dir_data, "store")
    target_path = os.path.join(store_path, "lung")
    ds.write_distributed_store(dir_cache=store_path, store="h5ad", dense=True)
    store = DistributedStore(cache_path=store_path)
    observed_nodes = np.unique(np.concatenate([
        x.obs[store._adata_ids_sfaira.cell_ontology_class]
        for x in store.adatas.values()
    ])).tolist()
    leaves_all = store.celltypes_universe.onto_cl.leaves
    effective_leaves = store.celltypes_universe.onto_cl.get_effective_leaves(x=observed_nodes)
    store.celltypes_universe.onto_cl.leaves = effective_leaves
    leaves1 = store.celltypes_universe.onto_cl.leaves
    store.celltypes_universe.write_target_universe(fn=target_path, x=effective_leaves)
    store2 = DistributedStore(cache_path=store_path)
    store2.celltypes_universe.load_target_universe(fn=target_path)
    leaves2 = store2.celltypes_universe.onto_cl.leaves
    assert len(leaves_all) > len(leaves1)
    assert len(set(leaves1).union(set(leaves2))) == len(leaves1)
    assert np.all([x in leaves1 for x in leaves2])
