import dask.array
import pandas as pd
import numpy as np
import scipy.sparse
import sparse
import torch
from typing import Dict, List, Tuple, Union


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
    cached_data: Union[None, Tuple[Tuple[torch.Tensor]]]
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
            xy = [self.__getitem_raw(idx=i) for i in range(self._len)]
            self._shapes = [len(z) for z in xy[0]]  # length of each data tuple, e.g. number of x and y tensors.
            # Expand observation dimension via stack if __getitem_raw collapsed the observation axis:
            xy = tuple(tuple(torch.cat([xy[n][i][j] for n in range(self._len)], dim=0) if len(xy[0][i][j].shape) > 1
                             else torch.stack([xy[n][i][j] for n in range(self._len)], dim=0)
                             for j in range(xi))
                       for i, xi in enumerate(self._shapes))
            self.cached_data = xy
        else:
            self.cached_data = None

    def __len__(self):
        return self._len

    def __getitem_cache(self, idx):
        # Flatten batch dim for torch.Dataset [not necessary for IteratableDataset]
        if len(idx) == 1:
            idx = idx[0]
        xy = tuple(tuple(self.cached_data[i][j][idx] for j in range(xi)) for i, xi in enumerate(self._shapes))
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
        # Flatten batch dim for torch.Dataset [not necessary for IteratableDataset]
        # TODO this might be inefficient, might need different solution.
        # TODO has consequences in setup_cache as well
        if xy[0][0].shape[0] == 1 and len(xy[0][0].shape) >= 2:
            xy = tuple(tuple(zz.squeeze(axis=0) for zz in z) for z in xy)
        xy = tuple(tuple(torch.from_numpy(zz) for zz in z) for z in xy)
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


class SfairaIterableDataset(torch.utils.data.IterableDataset):

    def __init__(self, iterator_fun, **kwargs):
        super(SfairaIterableDataset, self).__init__(**kwargs)
        self.iterator_fun = iterator_fun

    def __iter__(self):
        return self.iterator_fun()
