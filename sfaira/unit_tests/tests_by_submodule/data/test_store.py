import anndata
import dask.array
import h5py
import numpy as np
import os
import pytest
import scipy.sparse
from typing import List

from sfaira.data import load_store
from sfaira.versions.genomes.genomes import GenomeContainer

from sfaira.unit_tests.data_for_tests.loaders import ASSEMBLY_MOUSE, prepare_dsg, prepare_store


def _get_single_store(store_format: str):
    store_path = prepare_store(store_format=store_format)
    stores = load_store(cache_path=store_path, store_format=store_format)
    stores.subset(attr_key="organism", values=["mouse"])
    store = stores.stores["mouse"]
    return store


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
def test_fatal(store_format: str):
    """
    Test if basic methods abort.
    """
    store_path = prepare_store(store_format=store_format)
    stores = load_store(cache_path=store_path, store_format=store_format)
    stores.subset(attr_key="organism", values=["mouse"])
    store = stores.stores["mouse"]
    # Test both single and multi-store:
    for x in [store, stores]:
        _ = x.n_obs
        _ = x.n_vars
        _ = x.var_names
        _ = x.shape
        _ = x.obs
        _ = x.indices
        _ = x.genome_container


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
@pytest.mark.parametrize("as_sparse", [True, False])
def test_x_slice(store_format: str, as_sparse: bool):
    """
    Test if basic methods abort.
    """
    store = _get_single_store(store_format=store_format)
    data = store.X_slice(idx=np.arange(0, 5), as_sparse=as_sparse)
    assert data.shape[0] == 5
    if as_sparse:
        assert isinstance(data, scipy.sparse.csr_matrix)
    else:
        assert isinstance(data, np.ndarray)


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
@pytest.mark.parametrize("as_sparse", [True, False])
def test_adata_slice(store_format: str, as_sparse: bool):
    """
    Test if basic methods abort.
    """
    store = _get_single_store(store_format=store_format)
    data = store.adata_slice(idx=np.arange(0, 5), as_sparse=as_sparse)
    assert data.shape[0] == 5
    assert isinstance(data, anndata.AnnData)
    if as_sparse:
        assert isinstance(data.X, scipy.sparse.csr_matrix)
    else:
        assert isinstance(data.X, np.ndarray)


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
def test_data(store_format: str):
    """
    Test if the data exposed by the store are the same as in the original Dataset instance after streamlining.
    """
    # Run standard streamlining workflow on dsg and compare to object relayed via store.
    # Prepare dsg.
    dsg = prepare_dsg(load=True)
    # Prepare store.
    # Rewriting store to avoid mismatch of randomly generated data in cache and store.
    store_path = prepare_store(store_format=store_format, rewrite=False, rewrite_store=True)
    store = load_store(cache_path=store_path, store_format=store_format)
    store.subset(attr_key="doi_journal", values=["no_doi_mock1"])
    dataset_id = store.adata_by_key[list(store.indices.keys())[0]].uns["id"]
    adata_store = store.adata_by_key[dataset_id]
    x_store = store.data_by_key[dataset_id]
    adata_ds = dsg.datasets[dataset_id].adata
    x_ds = adata_ds.X.todense()
    if isinstance(x_store, dask.array.Array):
        x_store = x_store.compute()
    if isinstance(x_store, h5py.Dataset):
        # Need to load sparse matrix into memory if it comes from a backed anndata object.
        x_store = x_store[:, :]
    if isinstance(x_store, anndata._core.sparse_dataset.SparseDataset):
        # Need to load sparse matrix into memory if it comes from a backed anndata object.
        x_store = x_store[:, :]
    if isinstance(x_store, scipy.sparse.csr_matrix):
        x_store = x_store.todense()
    if isinstance(x_ds, anndata._core.sparse_dataset.SparseDataset):
        # Need to load sparse matrix into memory if it comes from a backed anndata object.
        x_ds = x_ds[:, :]
    if isinstance(x_ds, scipy.sparse.csr_matrix):
        x_ds = x_ds.todense()
    # Check that non-zero elements are the same:
    assert x_store.shape[0] == x_ds.shape[0]
    assert x_store.shape[1] == x_ds.shape[1]
    assert np.all(np.where(x_store > 0)[0] == np.where(x_ds > 0)[0]), (np.sum(x_store > 0), np.sum(x_ds > 0))
    assert np.all(np.where(x_store > 0)[1] == np.where(x_ds > 0)[1]), (np.sum(x_store > 0), np.sum(x_ds > 0))
    assert np.all(x_store - x_ds == 0.), (np.sum(x_store), np.sum(x_ds))
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
    store_path = prepare_store(store_format=store_format)
    config_path = os.path.join(store_path, "config_lung")
    store = load_store(cache_path=store_path, store_format=store_format)
    store.subset(attr_key="organism", values=["mouse"])
    store.subset(attr_key="assay_sc", values=["10x technology"])
    store.write_config(fn=config_path)
    store2 = load_store(cache_path=store_path, store_format=store_format)
    store2.load_config(fn=config_path + ".pickle")
    assert np.all(store.indices.keys() == store2.indices.keys())
    assert np.all([np.all(store.indices[k] == store2.indices[k])
                   for k in store.indices.keys()])


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
@pytest.mark.parametrize("idx", [np.arange(1, 10),
                                 np.concatenate([np.arange(30, 50), np.array([1, 4, 98])])])
@pytest.mark.parametrize("batch_size", [1, ])
@pytest.mark.parametrize("obs_keys", [["cell_type"]])
@pytest.mark.parametrize("randomized_batch_access", [True, False])
def test_generator_shapes(store_format: str, idx, batch_size: int, obs_keys: List[str], randomized_batch_access: bool):
    """
    Test generators queries do not throw errors and that output shapes are correct.
    """
    # Need to re-write because specific obs_keys are required:
    store_path = prepare_store(store_format=store_format, rewrite_store=True)
    store = load_store(cache_path=store_path, store_format=store_format)
    store.subset(attr_key="organism", values=["mouse"])
    gc = GenomeContainer(assembly=ASSEMBLY_MOUSE)
    gc.subset(**{"biotype": "protein_coding"})
    store.genome_container = gc
    g = store.generator(
        idx={"mouse": idx},
        batch_size=batch_size,
        obs_keys=obs_keys,
        randomized_batch_access=randomized_batch_access,
    )
    g = g.iterator
    nobs = len(idx) if idx is not None else store.n_obs
    batch_sizes = []
    x = None
    obs = None
    counter = 0
    for i, z in enumerate(g()):
        counter += 1
        x_i, obs_i = z
        if len(x_i.shape) == 1:
            # x is flattened if batch size is 1:
            assert batch_size == 1
            x_i = np.expand_dims(x_i, axis=0)
        assert x_i.shape[0] == obs_i.shape[0]
        if i == 0:
            x = x_i
            obs = obs_i
        batch_sizes.append(x_i.shape[0])
    assert counter > 0
    assert x.shape[1] == store.n_vars["mouse"], (x.shape, store.n_vars["mouse"])
    assert obs.shape[1] == len(obs_keys), (obs.shape, obs_keys)
    assert np.sum(batch_sizes) == nobs, (batch_sizes, nobs)
    assert x.shape[1] == gc.n_var, (x.shape, gc.n_var)
