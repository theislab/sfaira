import numpy as np
import os
import pytest
import time
from typing import List

from sfaira.data import load_store
from sfaira.versions.genomes import GenomeContainer
from sfaira.unit_tests.utils import cached_store_writing


MOUSE_GENOME_ANNOTATION = "Mus_musculus.GRCm38.102"
HUMAN_GENOME_ANNOTATION = "Homo_sapiens.GRCh38.102"

dir_data = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data")
dir_meta = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data", "meta")


"""
TODO tests from here on down require cached data for mouse lung
"""


@pytest.mark.parametrize("store_format", ["h5ad", "zarr"])
def test_fatal(store_format: str):
    """
    Test if basic methods abort.
    """
    store_path = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=MOUSE_GENOME_ANNOTATION,
                                      store_format=store_format)
    store = load_store(cache_path=store_path, store_format=store_format)
    store.subset(attr_key="organism", values=["mouse"])
    store.subset(attr_key="assay_sc", values=["10x sequencing"])
    _ = store.n_obs
    _ = store.n_vars
    _ = store.var_names
    _ = store.shape
    _ = store.obs
    _ = store.indices
    _ = store.genome_container
    _ = store.n_counts(idx=[1, 3])


@pytest.mark.parametrize("store_format", ["h5ad", "zarr"])
def test_config(store_format: str):
    """
    Test that data set config files can be set, written and recovered.
    """
    store_path = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=MOUSE_GENOME_ANNOTATION,
                                      store_format=store_format)
    config_path = os.path.join(store_path, "config_lung")
    store = load_store(cache_path=store_path, store_format=store_format)
    store.subset(attr_key="organism", values=["mouse"])
    store.subset(attr_key="assay_sc", values=["10x sequencing"])
    store.write_config(fn=config_path)
    store2 = load_store(cache_path=store_path, store_format=store_format)
    store2.load_config(fn=config_path + ".pickle")
    assert np.all(store.indices.keys() == store2.indices.keys())
    assert np.all([np.all(store.indices[k] == store2.indices[k]) for k in store.indices.keys()])


@pytest.mark.parametrize("store_format", ["h5ad", "zarr"])
@pytest.mark.parametrize("idx", [np.array([2, 1020, 3, 20000, 20100]),
                                 np.concatenate([np.arange(150, 200), np.array([1, 100, 2003, 33])])])
@pytest.mark.parametrize("batch_size", [1, 7])
@pytest.mark.parametrize("obs_keys", [[], ["cell_ontology_class"]])
@pytest.mark.parametrize("gc", [(None, {}), (MOUSE_GENOME_ANNOTATION, {"biotype": "protein_coding"})])
@pytest.mark.parametrize("randomized_batch_access", [True, False])
def test_generator_shapes(store_format: str, idx, batch_size: int, obs_keys: List[str], gc: tuple,
                          randomized_batch_access: bool):
    """
    Test generators queries do not throw errors and that output shapes are correct.
    """
    assembly, subset = gc
    store_path = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=MOUSE_GENOME_ANNOTATION,
                                      store_format=store_format)
    store = load_store(cache_path=store_path, store_format=store_format)
    store.subset(attr_key="organism", values=["mouse"])
    if assembly is not None:
        gc = GenomeContainer(assembly=assembly)
        gc.subset(**subset)
        store.genome_container = gc
    g = store.generator(
        idx=idx,
        batch_size=batch_size,
        obs_keys=obs_keys,
        randomized_batch_access=randomized_batch_access,
    )
    nobs = len(idx) if idx is not None else store.n_obs
    batch_sizes = []
    t0 = time.time()
    for i, z in enumerate(g()):
        x_i, obs_i = z
        assert x_i.shape[0] == obs_i.shape[0]
        if i == 0:
            x = x_i
            obs = obs_i
        batch_sizes.append(x_i.shape[0])
    tdelta = time.time() - t0
    print(f"time for iterating over generator:"
          f" {tdelta}s for {np.sum(batch_sizes)} cells in {len(batch_sizes)} batches,"
          f" {tdelta / len(batch_sizes)}s per batch.")
    assert x.shape[1] == store.n_vars, (x.shape, store.n_vars)
    assert obs.shape[1] == len(obs_keys), (obs.shape, obs_keys)
    assert np.sum(batch_sizes) == nobs, (batch_sizes, nobs)
    if assembly is not None:
        assert x.shape[1] == gc.n_var, (x.shape, gc.n_var)
