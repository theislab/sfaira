import dask.array
import pandas as pd
import numpy as np
import scipy.sparse
import sparse
import torch
from typing import Dict, List, Tuple, Union


def flatten_single_obs(xy, is_torch: bool):
    # Flatten batch dim for torch.Dataset [not necessary for IteratableDataset]
    # TODO this might be inefficient, might need different solution.
    # TODO has consequences in setup_cache as well
    if isinstance(xy, dict):
        # Assume is dictionary of arrays.
        if xy[list(xy.keys())[0]].shape[0] == 1 and len(xy[list(xy.keys())[0]].shape) >= 2:
            xy = dict([(k, v.squeeze(axis=0)) for k, v in xy.items()])
        if not is_torch:
            xy = dict([(k, torch.from_numpy(v)) for k, v in xy.items()])
    elif isinstance(xy[0], np.ndarray):
        # Assume is tuple of arrays.
        if xy[0].shape[0] == 1 and len(xy[0].shape) >= 2:
            xy = tuple(z.squeeze(axis=0) for z in xy)
        if not is_torch:
            xy = tuple(torch.from_numpy(z) for z in xy)
    else:
        # Assume is tuple of tuples.
        if xy[0][0].shape[0] == 1 and len(xy[0][0].shape) >= 2:
            xy = tuple(tuple(zz.squeeze(axis=0) for zz in z) for z in xy)
        if not is_torch:
            xy = tuple(tuple(torch.from_numpy(zz) for zz in z) for z in xy)
    return xy


class IndexDataset(torch.utils.data.Dataset):

    """
    Only yields index of selected observations back as tensor, can be used as proxy data set if data is already full
    copied to GPU, for example.
    """

    def __init__(self, number_steps: int, **kwargs):
        super(IndexDataset, self).__init__()
        self._number_steps = number_steps

    def __len__(self):
        return self._number_steps

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        if isinstance(idx, int):
            idx = [idx]
        xy = ((torch.from_numpy(np.asarray([idx])),),)
        return xy


class SfairaDataset(torch.utils.data.Dataset):

    _shapes: List[int]
    cached_data: Union[None, Union[Dict[str, torch.Tensor], Tuple[torch.Tensor], Tuple[Tuple[torch.Tensor]]]]
    cache_element_structure: Union[None, Union[List[str], int, List[int]]]
    cache_structure: Union[None, str]
    use_cache: bool

    """
    Note on caching:

    Caches can be iteratively build during the first epoch, this does however constrain the data pipeline to be serial
    during the first epoch which results in an awkward interface to torch lightning if parallelization is desired
    in subsequent epoch. Therefore, we implemented a one-off pre-loading of the cache which has the disadvantage of an
    overhead at initialization. We may provide sequential caching during first epoch in the future.
    See also: https://discuss.pytorch.org/t/best-practice-to-cache-the-entire-dataset-during-first-epoch/19608
    """

    def __init__(
            self,
            map_fn,
            x: Union[np.ndarray, scipy.sparse.spmatrix, dask.array.Array],
            obs: pd.DataFrame,
            obsm: Dict[str, Union[np.ndarray, scipy.sparse.spmatrix, dask.array.Array]],
            use_cache: bool = False,
            **kwargs):
        """

        :param map_fn:
        :param obs:
        :param obsm:
        :param x:
        :param dask_to_memory: Whether to convert X to a sparse matrix in memory if it is a dask array.
            This implies that x is also moved from disk to memory if it is currently a lazy dask representation of
            a zarr array.
        :param use_cache: Whether to cache transformed data. Use with caution: Copy of full data will be in memory!
        :param kwargs:
        """
        super(SfairaDataset, self).__init__()
        self.map_fn = map_fn
        self.obs = obs
        self.obsm = obsm
        self.emit_obsm = len(obsm.keys()) > 0
        if isinstance(x, scipy.sparse.spmatrix) or isinstance(x, sparse.spmatrix):
            self._using_sparse = True
        else:
            self._using_sparse = False
        self.x = x
        assert self.x.shape[0] == self.obs.shape[0], (self.x.shape, self.obs.shape)
        self._len = self.x.shape[0]
        # Caching:
        self.setup_cache(use_cache)

    def setup_cache(self, use_cache):
        """
        Sets cache up.

        The cache has the same structure as a sample from the data set but each tensor has the full length of the
        dataset as observation dimension. This implies that no stacking needs to be performed during query, only
        integer-based query of a dataset-shaped tensor.

        :param use_cache: Whether to use cache. If true, sets .cached_data: loads full transformed data.
        """
        self.use_cache = use_cache
        if self.use_cache:
            xy = [self.__getitem_raw(idx=[i]) for i in range(self._len)]
            # Expand observation dimension via stack if __getitem_raw collapsed the observation axis:
            if isinstance(xy[0], dict):
                # Assume is dictionary of tensors.
                self.cache_element_structure = list(xy[0].keys())  # tensor keys
                xy = dict([(k, torch.cat([xy[n][k] for n in range(self._len)], dim=0))
                           if len(xy[0][k].shape) > 1 else
                           (k, torch.stack([xy[n][k] for n in range(self._len)], dim=0))
                           for k in self.cache_element_structure])
                self.cache_structure = "dict_tensor"
            elif isinstance(xy[0][0], np.ndarray):
                # Assume is tuple of tensors.
                self.cache_element_structure = len(xy[0])  # length of data tuple
                xy = tuple(torch.cat([xy[n][i] for n in range(self._len)], dim=0) if len(xy[0][i].shape) > 1
                           else torch.stack([xy[n][i] for n in range(self._len)], dim=0)
                           for i in range(self.cache_element_structure))
                self.cache_structure = "tuple_tensor"
            else:
                # Assume is tuple of tuples of tensors.
                self.cache_element_structure = [len(z) for z in xy[0]]
                xy = tuple(tuple(torch.cat([xy[n][i][j] for n in range(self._len)], dim=0) if len(xy[0][i][j].shape) > 1
                                 else torch.stack([xy[n][i][j] for n in range(self._len)], dim=0)
                                 for j in range(xi))
                           for i, xi in enumerate(self.cache_element_structure))
                self.cache_structure = "tuple_tuple_tensor"
            self.cached_data = xy
        else:
            self.cached_data = None

    def __len__(self):
        return self._len

    def __getitem_cache(self, idx):
        # Flatten batch dim for torch.Dataset [not necessary for IteratableDataset]
        if len(idx) == 1:
            idx = idx[0]
        if self.cache_structure == "dict_tensor":
            xy = dict([(k, self.cached_data[k][idx]) for k in self.cache_element_structure])
        elif self.cache_structure == "tuple_tensor":
            xy = tuple(self.cached_data[i][idx] for i in range(self.cache_element_structure))
        elif self.cache_structure == "tuple_tuple_tensor":
            xy = tuple(tuple(self.cached_data[i][j][idx] for j in range(xi))
                       for i, xi in enumerate(self.cache_element_structure))
        else:
            assert False
        return xy

    def __getitem_raw(self, idx):
        x = self.x[idx, :]
        if self._using_sparse:
            x = x.toarray()
        obs = self.obs.iloc[idx, :]
        if self.emit_obsm:
            obsm = dict([(k, v[idx, :]) for k, v in self.obsm.items()])
            data_tuple = (x, obs, obsm)
        else:
            data_tuple = (x, obs)
        xy = self.map_fn(*data_tuple)
        xy = flatten_single_obs(xy, is_torch=False)
        return xy

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        if isinstance(idx, int):
            idx = [idx]
        if self.use_cache:
            xy = self.__getitem_cache(idx=idx)
        else:
            xy = self.__getitem_raw(idx=idx)
        return xy


class SfairaMultiDataset(torch.utils.data.Dataset):

    datasets: List[torch.utils.data.Dataset]
    n_obs: List[int]
    n_obs_cumulative_lb: List[int]
    n_obs_cumulative_ub: List[int]
    map_fn_merge: callable

    def __init__(self, datasets: Dict[str, torch.utils.data.Dataset], map_fn_merge: Union[None, callable], **kwargs):
        """

        :param datasets:
        :param kwargs:
        """
        assert map_fn_merge is not None, "set map_fn_merge"
        super(SfairaMultiDataset, self).__init__()
        self.dataset_keys = list(datasets.keys())
        self.datasets = [datasets[k] for k in self.dataset_keys]
        self.map_fn_merge = map_fn_merge
        self.n_obs = [len(x) for x in self.datasets]
        self.n_obs_cumulative_lb = [0] + np.cumsum(self.n_obs).tolist()[:-1]
        self.n_obs_cumulative_ub = np.cumsum(self.n_obs).tolist()

    def __len__(self):
        return sum(self.n_obs)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        if isinstance(idx, int):
            idx = [idx]
        # Map indices to individual data sets:
        idx = [[int(i - lb) for i in idx if lb <= i < ub]
               for lb, ub in zip(self.n_obs_cumulative_lb, self.n_obs_cumulative_ub)]
        # Note: the calls to the individual datasets already yield torch tensors, not numpy arrays!
        xy = self.map_fn_merge([(k, x[i]) for k, x, i in zip(self.dataset_keys, self.datasets, idx) if len(i) > 0])
        xy = flatten_single_obs(xy, is_torch=True)
        return xy


class SfairaIterableDataset(torch.utils.data.IterableDataset):

    def __init__(self, iterator_fun, **kwargs):
        super(SfairaIterableDataset, self).__init__(**kwargs)
        self.iterator_fun = iterator_fun

    def __iter__(self):
        return self.iterator_fun()
