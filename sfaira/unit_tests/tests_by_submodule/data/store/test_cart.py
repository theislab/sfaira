import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse
from typing import List

from sfaira.consts import AdataIdsSfaira
from sfaira.unit_tests.tests_by_submodule.data.store.utils import _get_cart


@pytest.mark.parametrize("feature_space", ["single", "multi"])
@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
def test_properties(store_format: str, feature_space: str):
    """Tests if properties of carts are available and of the right format."""
    idx = np.arange(0, 5)
    kwargs = {"idx": {"Mus musculus": idx}, "obs_keys": ["cell_type"], "randomized_batch_access": False,
              "return_dense": False}
    cart = _get_cart(store_format=store_format, feature_space=feature_space, **kwargs)

    # Meta data:
    assert cart.n_obs > len(idx)
    assert cart.n_obs_selected == len(idx)
    obs_idx = cart.obs_idx
    if feature_space != "single":
        obs_idx = obs_idx["Mus musculus"]
    if store_format == "dao":
        cart.move_to_memory()
    assert np.all(obs_idx == idx)
    # Data arrays:
    arr = [
        cart.adata,
        cart.x,
        cart.obs,
        cart.var,
    ]
    if feature_space != "single":
        arr = [x["Mus musculus"] for x in arr]
    assert np.all([x.shape[0] == len(idx) for x in arr[:3]])
    assert isinstance(arr[0], anndata.AnnData)
    assert scipy.sparse.issparse(arr[0].X)
    assert scipy.sparse.issparse(arr[1])
    assert isinstance(arr[2], pd.DataFrame)
    assert isinstance(arr[3], pd.DataFrame)
    # Other
    _ = cart.n_var


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
@pytest.mark.parametrize("idx", [np.arange(1, 10),
                                 np.concatenate([np.arange(30, 50), np.array([1, 4, 98])])])
@pytest.mark.parametrize("obs_keys", [["cell_type"]])
@pytest.mark.parametrize("randomized_batch_access", [True, False])
def test_data(store_format: str, idx, obs_keys: List[str], randomized_batch_access: bool):
    """
    Test generators queries do not throw errors and that output shapes are correct.
    """
    # Need to re-write because specific obs_keys are required:
    kwargs = {"idx": {"Mus musculus": idx}, "obs_keys": obs_keys, "randomized_batch_access": randomized_batch_access}
    cart = _get_cart(store_format=store_format, feature_space="multi", **kwargs)
    it = cart.iterator
    batch_sizes = []
    x = None
    counter = 0
    for i, z in enumerate(it()):
        counter += 1
        x_i, = z[0]
        if len(x_i.shape) == 1:
            x_i = np.expand_dims(x_i, axis=0)
        if i == 0:
            x = x_i
        batch_sizes.append(x_i.shape[0])
    assert counter > 0
    assert x.shape[1] == cart.n_var["Mus musculus"], (x.shape, cart.n_var)
    assert np.sum(batch_sizes) == len(idx), (batch_sizes, len(idx))


@pytest.mark.parametrize("store_format", ["h5ad", "dao"])
@pytest.mark.parametrize("idx", [None, np.array([1, 4, 98])])
def test_schedule_blocked(store_format: str, idx):
    """
    Test if schedule "blocks" works.
    """
    block_col = AdataIdsSfaira().cell_type
    obs_keys = [block_col]

    def map_fn(x, obs):
        return (obs, ),

    kwargs = {"idx": {"Mus musculus": idx}, "obs_keys": obs_keys, "randomized_batch_access": False,
              "batch_schedule": "blocks", "grouping": block_col, "batch_size": 0}
    cart = _get_cart(store_format=store_format, feature_space="single", map_fn=map_fn, **kwargs)
    it = cart.iterator
    block_vals = cart.obs[block_col].values
    batches = []
    for i, z in enumerate(it()):
        obs_i, = z[0]
        # Check that this batch is one single block:
        assert len(np.unique(obs_i[block_col].values)) == 1, obs_i
        batches.append(obs_i)
    assert len(batches) > 0
    # Check that one batch was emitted per block:
    assert len(np.unique(block_vals)) == len(batches), batches
    # Check that the total number of observations across blocks is correct:
    if idx is None:
        target_nobs = cart.n_obs
    else:
        target_nobs = len(idx)
    assert np.sum([y.shape[0] for y in batches]) == target_nobs


@pytest.mark.parametrize(
    "adaptor", ["python", "tensorflow", "torch", "torch-loader", "torch-iter", "torch-iter-loader"]
)
@pytest.mark.parametrize("shuffle_buffer_size", [0, 1000])
def test_adaptors(adaptor: str, shuffle_buffer_size: int):
    """
    Test if framework-specific generator adpators yield batches.
    """
    idx = np.arange(0, 10)

    def map_fn(x_, obs_):
        """
        Note: Need to convert to numpy in output because torch does not accept dask.
        """
        return (np.asarray(x_[:, :2]),),

    kwargs = {"idx": {"Mus musculus": idx}, "obs_keys": [], "randomized_batch_access": False, "retrieval_batch_size": 2,
              "map_fn": map_fn}
    cart = _get_cart(store_format="dao", feature_space="single", **kwargs)

    if adaptor == "python":
        kwargs = {}
    elif adaptor == "tensorflow":
        import tensorflow as tf

        kwargs = {"output_signature": (
            tf.TensorSpec(shape=(2,), dtype=tf.float32),
        )}
    elif adaptor in ["torch", "torch-loader", "torch-iter-loader", "torch-iter"]:
        kwargs = {}
    else:
        assert False
    it = cart.adaptor(generator_type=adaptor, shuffle_buffer=shuffle_buffer_size, **kwargs)
    if adaptor == "tensorflow":
        it = iter(it.range(2))
    if adaptor in ["torch", "torch-iter"]:
        from torch.utils.data import DataLoader
        it = list(DataLoader(it))
        it = iter(it)
    if adaptor in ["torch-loader", "torch-iter-loader"]:
        import torch
        it = iter(list(it))
    _ = next(it)
