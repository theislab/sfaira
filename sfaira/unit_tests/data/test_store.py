import anndata
import dask.array
import numpy as np
import os
import pytest
import scipy.sparse
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
Tests from here on down require cached data for mouse lung
"""


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
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


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
@pytest.mark.parametrize("dataset", ["mouse_lung_2019_10xsequencing_pisco_022_10.1101/661728"])
def test_data(store_format: str, dataset: str):
    """
    Test if the data exposed by the store are the same as in the original Dataset instance after streamlining.
    """
    store_path, ds = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=MOUSE_GENOME_ANNOTATION,
                                          store_format=store_format, return_ds=True)
    store = load_store(cache_path=store_path, store_format=store_format)
    dataset_key_reduced = dataset.split("_10.")[0]
    store.subset(attr_key="id", values=[dataset_key_reduced])
    adata_store = store.adata_by_key[dataset]
    adata_ds = ds.datasets[dataset].adata
    # Check .X
    x_store = adata_store.X
    x_ds = adata_ds.X.todense()
    if isinstance(x_store, dask.array.Array):
        x_store = x_store.compute()
    if isinstance(x_store, anndata._core.sparse_dataset.SparseDataset):
        # Need to load sparse matrix into memory if it comes from a backed anndata object.
        x_store = x_store[:, :]
    if isinstance(x_store, scipy.sparse.csr_matrix):
        x_store = x_store.todense()
    # Check that non-zero elements are the same:
    assert np.all(np.where(x_store > 0)[0] == np.where(x_ds > 0)[0])
    assert np.all(np.where(x_store > 0)[1] == np.where(x_ds > 0)[1])
    assert np.all(x_store - x_ds == 0.)
    assert x_store.dtype == x_ds.dtype
    # Note: Do not run test on sum across entire object if dtype is float32 as this can result in test failures because
    # of float overflows.
    # Check .obs
    obs_store = adata_store.obs
    obs_ds = adata_ds.obs
    assert np.all(obs_store.columns == obs_ds.columns), (obs_store.columns, obs_ds.columns)
    for k, v in obs_store.items():
        assert np.all(np.asarray(v.values.tolist()) == np.asarray(obs_ds[k].values.tolist()))
    # Check .var
    var_store = adata_store.var
    var_ds = adata_ds.var
    assert np.all(var_store.columns == var_ds.columns), (var_store.columns, var_ds.columns)
    for k, v in var_store.items():
        assert np.all(np.asarray(v.values.tolist()) == np.asarray(var_ds[k].values.tolist()))
    # Check .uns
    uns_store = adata_store.uns
    uns_ds = adata_ds.uns
    assert np.all(uns_store.keys() == uns_ds.keys()), (uns_store.keys(), uns_ds.keys())
    for k, v in uns_store.items():
        assert np.all(v == uns_ds[k])


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
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


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
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
    x = None
    obs = None
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
