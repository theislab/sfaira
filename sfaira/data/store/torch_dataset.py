import dask.array
import pandas as pd
import numpy as np
import scipy.sparse
import sparse
import torch
from typing import Union


class SfairaDataset(torch.utils.data.Dataset):

    def __init__(
            self,
            map_fn,
            x: Union[np.ndarray, scipy.sparse.spmatrix, dask.array.Array],
            obs: pd.DataFrame,
            **kwargs):
        """

        :param map_fn:
        :param obs:
        :param x:
        :param dask_to_memory: Whether to convert X to a sparse matrix in memory if it is a dask array.
            This implies that x is also moved from disk to memory if it is currently a lazy dask representation of
            a zarr array.
        :param kwargs:
        """
        super(SfairaDataset, self).__init__()
        self.map_fn = map_fn
        self.obs = obs
        if isinstance(x, scipy.sparse.spmatrix) or isinstance(x, sparse.spmatrix):
            self._using_sparse = True
        else:
            self._using_sparse = False
        self.x = x
        assert self.x.shape[0] == self.obs.shape[0], (self.x.shape, self.obs.shape)
        self._len = self.x.shape[0]

    def __len__(self):
        return self._len

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()
        if isinstance(idx, int):
            idx = [idx]
        x = self.x[idx, :]
        if self._using_sparse:
            x = x.toarray()
        obs = self.obs.iloc[idx, :]
        xy = self.map_fn(x, obs)
        # Flatten batch dim for torch.Dataset [not necessary for IteratableDataset]
        # TODO this might be inefficient, might need different solution.
        if xy[0][0].shape[0] == 1 and len(xy[0][0].shape) >= 2:
            xy = tuple(tuple(zz.squeeze(axis=0) for zz in z) for z in xy)
        xy = tuple(tuple(torch.from_numpy(zz) for zz in z) for z in xy)
        return xy


class SfairaIterableDataset(torch.utils.data.IterableDataset):

    def __init__(self, iterator_fun, **kwargs):
        super(SfairaIterableDataset, self).__init__(**kwargs)
        self.iterator_fun = iterator_fun

    def __iter__(self):
        return self.iterator_fun()
