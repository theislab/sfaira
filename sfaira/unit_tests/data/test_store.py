import numpy as np
import os
import pytest
from typing import List

from sfaira.data import DistributedStore
from sfaira.versions.genomes import GenomeContainer
from sfaira.unit_tests.utils import cached_store_writing


MOUSE_GENOME_ANNOTATION = "Mus_musculus.GRCm38.102"

dir_data = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data")
dir_meta = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data/meta")


"""
TODO tests from here on down require cached data for mouse lung
"""


def test_config():
    """
    Test that data set config files can be set, written and recovered.
    """
    store_path = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=MOUSE_GENOME_ANNOTATION)
    config_path = os.path.join(store_path, "lung")
    store = DistributedStore(cache_path=store_path)
    store.subset(attr_key="assay_sc", values=["10x sequencing"])
    store.subset_cells(attr_key="assay_sc", values=["10x sequencing"])
    store.write_config(fn=config_path)
    store2 = DistributedStore(cache_path=store_path)
    store2.load_config(fn=config_path)
    assert np.all(store.indices.keys() == store2.indices.keys())
    assert np.all([np.all(store.indices[k] == store2.indices[k]) for k in store.indices.keys()])


def test_type_targets():
    """
    Test that target leave nodes can be set, written and recovered.
    """
    store_path = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=MOUSE_GENOME_ANNOTATION)
    target_path = os.path.join(store_path, "lung")
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


@pytest.mark.parametrize("batch_size", [1, 10])
@pytest.mark.parametrize("obs_keys", [[], ["cell_ontology_class"]])
@pytest.mark.parametrize("continuous_batches", [True])
@pytest.mark.parametrize("assembly", [None, MOUSE_GENOME_ANNOTATION])
@pytest.mark.parametrize("subset", [{}, {"biotype": "protein_coding"}])
def test_generator_shapes(batch_size: int, obs_keys: List[str], continuous_batches: bool, assembly: str, subset: dict):
    """
    Test generators queries do not throw errors and that output shapes are correct.
    """
    if assembly is not None:
        gc = GenomeContainer(assembly=assembly)
        gc.subset(**subset)
    else:
        gc = None
    store_path = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=MOUSE_GENOME_ANNOTATION)
    store = DistributedStore(cache_path=store_path)
    store.genome_container = gc
    g = store.generator(
        batch_size=batch_size,
        obs_keys=obs_keys,
        continuous_batches=continuous_batches,
    )
    x, obs = next(g)
    assert x.shape[0] == batch_size, (x.shape, batch_size)
    assert obs.shape[0] == batch_size, (obs.shape, batch_size)
    assert x.shape[1] == store.n_vars, (x.shape, store.n_vars)
    assert obs.shape[1] == len(obs_keys), (x.shape, obs_keys)
    if assembly is not None:
        assert x.shape[1] == gc.n_var, (x.shape, gc.n_var)
