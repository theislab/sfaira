import abc
from functools import partial
from typing import Dict, List, Tuple, Union

import anndata
import dask.array
import numpy as np
import pandas as pd
import scipy.sparse

from sfaira.data.store.batch_schedule import BATCH_SCHEDULE, BatchDesignBase
from sfaira.data.store.carts.base import CartBase
from sfaira.data.store.carts.dataset_utils import _ShuffleBuffer, _DatasetIteratorRepeater
from sfaira.data.store.carts.utils import split_batch


class CartSingle(CartBase):

    """
    Cart for a DistributedStoreSingleFeatureSpace().
    """

    _obs_idx: Union[np.ndarray, None]
    _batch_schedule_name: str
    batch_size: int
    map_fn: callable
    map_fn_base: callable
    obs_keys: List[str]
    obsm: dict
    schedule: BatchDesignBase
    var: pd.DataFrame
    var_idx: Union[None, np.ndarray]

    def __init__(self, obs_idx, obs_keys, var, var_idx=None, batch_schedule="base", batch_size=1, map_fn=None,
                 map_fn_base=None, obsm={}, **kwargs):
        """

        :param batch_schedule: A valid batch schedule name or a class that inherits from BatchDesignBase.

            - "base": sfaira.data.store.batch_schedule.BatchDesignBasic
            - "balanced": sfaira.data.store.batch_schedule.BatchDesignBalanced
            - "blocks": sfaira.data.store.batch_schedule.BatchDesignBlocks
            - "full":  sfaira.data.store.batch_schedule.BatchDesignFull
            - class: batch_schedule needs to be a class (not instance), subclassing BatchDesignBase.

        :param batch_size: Emission batch size. Must be 1 or 0.
        :param map_fn: Map function to apply to output tuple of raw generator. Each draw i from the generator is then:
            `yield map_fn(x[i, var_idx], obs[i, obs_keys])`. This map function is designed for torch Datasets, ie
            supply this if "torch" or "torch-loader" is used as an adaptor.
        :param map_fn: Map function to apply to output tuple of raw generator. Each draw i from the generator is then:
            `yield map_fn(x[i, var_idx], obs[i, obs_keys])`. This map function is designed for torch IteratableDatasets,
            ie supply this if "torch-iter" or "torch-loader-iter" is used as an adaptor.
        :param obs_idx: np.ndarray: The observations to emit.
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adata_by_key.
        :parm obsm: Empty dict or dict with additional observation-indexed arrays that are in memory.
        :param var_idx: The features to emit.
        :parm split_to_obs: Whether to split tensors to observation-wise slices at the emission stage of the generator.
        """
        self._obs_idx = None
        self._batch_schedule_name = batch_schedule
        self.batch_size = batch_size
        self.map_fn = map_fn
        self.map_fn_base = map_fn_base
        self.obs_keys = obs_keys
        self.var = var
        self.var_idx = var_idx
        if isinstance(batch_schedule, str):
            if batch_schedule == "blocks":
                # Extract group column values from .obs.
                assert "grouping" in kwargs.keys()
                group_vals = self._obs_full[kwargs["grouping"]].values
                kwargs["grouping"] = group_vals
            batch_schedule = BATCH_SCHEDULE[batch_schedule]
        else:
            if not isinstance(batch_schedule, BatchDesignBase):
                raise ValueError("Either supply a name of a valid batch schedule or a class inherting from "
                                 f"BatchDesignBase to the constructor of GeneratorSingle. Found {type(batch_schedule)}.")
        self.schedule = batch_schedule(**kwargs)
        self.obs_idx = obs_idx  # This needs to be set after .schedule.
        self.obsm = obsm

    def adaptor(
            self,
            generator_type: str,
            dataset_kwargs: dict = None,
            shuffle_buffer: int = 0,
            repeat: int = 1,
            **kwargs
    ):
        """
        See documentation of self._adaptor().
        """
        iter_kwargs = {"shuffle_buffer": shuffle_buffer, "repeat": repeat}
        return self._adaptor(generator_type=generator_type, dataset_kwargs=dataset_kwargs, iter_kwargs=iter_kwargs,
                             **kwargs)

    def adaptor_torch(self, dataset_kwargs, loader, **kwargs):
        from torch.utils.data import DataLoader
        # Only import this module if torch is used to avoid strict torch dependency:
        from sfaira.data.store.torch_dataset import SfairaDataset

        g = SfairaDataset(map_fn=self.map_fn_base, obs=self.obs, obsm=self.obsm, x=self.x, **dataset_kwargs)
        if loader:
            g = DataLoader(g, **kwargs)
        return g

    def adaptor_torch_iter(self, loader, repeat, shuffle_buffer, **kwargs):
        from torch.utils.data import DataLoader
        # Only import this module if torch is used to avoid strict torch dependency:
        from sfaira.data.store.torch_dataset import SfairaIterableDataset

        g = SfairaIterableDataset(iterator_fun=partial(self.iterator, repeat=repeat, shuffle_buffer=shuffle_buffer))
        if loader:
            g = DataLoader(g, **kwargs)
        return g

    @property
    def _emit_obsm(self):
        return len(self.obsm.keys()) > 0

    @property
    def _obs_full(self):
        """
        Full meta data matrix (cells x meta data) that is emitted in batches by .iterator().
        """
        raise NotImplementedError()

    def _validate_idx(self, idx: Union[np.ndarray, list]) -> np.ndarray:
        """
        Validate global index vector.
        """
        if len(idx) > 0:
            assert np.max(idx) < self.n_obs, f"maximum of supplied index vector {np.max(idx)} exceeds number of " \
                                             f"modelled observations {self.n_obs}"
            assert len(idx) == len(np.unique(idx)), f"repeated indices in idx: {len(idx) - len(np.unique(idx))}"
            if isinstance(idx, np.ndarray):
                assert len(idx.shape) == 1, idx.shape
                assert idx.dtype == int
            else:
                assert isinstance(idx, list)
                assert isinstance(idx[0], int) or isinstance(idx[0], int)
        idx = np.asarray(idx)
        return idx

    @property
    def _x_full(self):
        """
        Full data matrix (cells x features) that is emitted in batches by .iterator().
        """
        raise NotImplementedError()

    @property
    def adata(self) -> anndata.AnnData:
        """
        Assembles a slice of this cart based on .obs_idx as an anndata instance.
        """
        return anndata.AnnData(X=self.x, obs=self.obs, var=self.var)

    @property
    def n_obs(self) -> int:
        """Total number of observations in cart."""
        raise NotImplementedError()

    @property
    def n_obs_selected(self) -> int:
        """Total number of selected observations in cart."""
        return len(self.obs_idx)

    @property
    def obs_idx(self):
        """
        Integer observation indices to select from cart. These will be emitted if data is queried.
        """
        return self._obs_idx

    @obs_idx.setter
    def obs_idx(self, x):
        if x is None:
            x = np.arange(0, self.n_obs)
        else:
            x = self._validate_idx(x)
            x = np.sort(x)
        # Only reset if they are actually different:
        if self._obs_idx is None or \
                (self._obs_idx is not None and len(x) != len(self._obs_idx)) or \
                np.any(x != self._obs_idx):
            self._obs_idx = x
            self.schedule.idx = x

    @abc.abstractmethod
    def _iterator(self):
        pass

    def iterator(self, repeat: int = 1, shuffle_buffer: int = 0):
        """
        Iterator over data matrix and meta data table, yields batches of data points.
        """
        if shuffle_buffer > 2 and self.batch_size == 1:
            iterator = _ShuffleBuffer(self._iterator, shuffle_buffer).iterator
        else:
            iterator = self._iterator

        return _DatasetIteratorRepeater(iterator, n_repeats=repeat).iterator()


class CartAnndata(CartSingle):

    """
    Cart for a DistributedStoreAnndata().
    """

    adata_dict: Dict[str, anndata._core.views.ArrayView]
    return_dense: bool
    single_object: bool

    def __init__(self, adata_dict, return_dense=False, **kwargs):
        self.return_dense = return_dense
        self.single_object = len(adata_dict.keys()) == 1
        self.adata_dict = adata_dict
        super(CartAnndata, self).__init__(**kwargs)

    @property
    def _obs_full(self):
        """
        Full meta data matrix (cells x meta data) that is emitted in batches by .iterator().
        """
        return pd.concat([
            v.obs[self.obs_keys]
            for v in self.adata_dict.values()
        ], axis=0, join="inner", ignore_index=True, copy=False)

    @property
    def _x_full(self):
        """
        Full data matrix (cells x features) that is emitted in batches by .iterator().
        """
        x = scipy.sparse.vstack([
            self._parse_array(v.X, return_dense=False)
            for v in self.adata_dict.values()
        ])
        if self.return_dense and isinstance(x, scipy.sparse.spmatrix):
            x = x.todense()
        return x

    def _iterator(self):
        """
        Iterator over data matrix and meta data table, yields batches of data points.
        """
        # Speed up access to single object by skipping index overlap operations:
        for idx_i in self.schedule.design:
            if len(idx_i) > 0:
                # Match adata objects that overlap to batch:
                idx_i_dict = self._obs_idx_dict_query(idx=idx_i)
                if self.batch_size == 1:
                    # Emit each data set separately and avoid concatenation into larger chunks for emission.
                    for k, v in idx_i_dict.items():
                        # I) Prepare data matrix.
                        x = self.adata_dict[k].X[v, :]
                        x = self._parse_array(x=x, return_dense=self.return_dense)
                        # Prepare .obs.
                        obs = self.adata_dict[k].obs[self.obs_keys].iloc[v, :]
                        if self._emit_obsm:
                            obsm = dict([(k, v[idx_i, :]) for k, v in self.obsm.items()])
                            map_fn_args = (x, obs, obsm)
                        else:
                            map_fn_args = (x, obs)
                        data_tuple = self.map_fn(*map_fn_args)
                        for data_tuple_i in split_batch(x=data_tuple):
                            yield data_tuple_i
                else:
                    # Concatenates slices first before returning. Note that this is likely slower than emitting by
                    # observation in most scenarios.
                    # I) Prepare data matrix.
                    x = [
                        self._parse_array(self.adata_dict[k].X[v, :], return_dense=self.return_dense)
                        for k, v in idx_i_dict.items()
                    ]
                    is_dense = isinstance(x[0], np.ndarray)
                    # Concatenate blocks in observation dimension:
                    if len(x) > 1:
                        if is_dense:
                            x = np.concatenate(x, axis=0)
                        else:
                            x = scipy.sparse.vstack(x)
                    else:
                        x = x[0]
                    # Prepare .obs.
                    obs = pd.concat([
                        self.adata_dict[k].obs[self.obs_keys].iloc[v, :]
                        for k, v in idx_i_dict.items()
                    ], axis=0, join="inner", ignore_index=True, copy=False)
                    if self._emit_obsm:
                        obsm = dict([(k, v[idx_i, :]) for k, v in self.obsm.items()])
                        map_fn_args = (x, obs, obsm)
                    else:
                        map_fn_args = (x, obs)
                    data_tuple = self.map_fn(*map_fn_args)
                    yield data_tuple

    def move_to_memory(self):
        """
        No action.
        """
        pass

    @property
    def n_obs(self) -> int:
        """Total number of observations in cart."""
        return int(np.sum([v.n_obs for v in self.adata_dict.values()]))

    @property
    def n_var(self) -> int:
        """Total number of features defined for return in cart."""
        if self.var_idx is None:
            return self.adata_dict[list(self.adata_dict.keys())[0]].n_var
        else:
            return len(self.var_idx)

    @property
    def obs(self):
        """
        Selected meta data matrix (cells x meta data) that is emitted in batches by .iterator().
        """
        idx_dict = self._obs_idx_dict_query(idx=self.obs_idx)
        return pd.concat([
            self.adata_dict[k].obs[self.obs_keys].iloc[v, :]
            for k, v in idx_dict.items()
        ], axis=0, join="inner", ignore_index=True, copy=False)

    @property
    def x(self):
        """
        Selected data matrix (cells x features) that is emitted in batches by .iterator().
        """
        # Assumes that .X are scipy.sparse.csr or coo
        idx_dict = self._obs_idx_dict_query(idx=self.obs_idx)
        x = scipy.sparse.vstack([
            self._parse_array(self.adata_dict[k].X[v, :], return_dense=False)
            for k, v in idx_dict.items()
        ])
        if self.return_dense and isinstance(x, scipy.sparse.spmatrix):
            x = x.todense()
        return x

    # Methods that are specific to this child class:

    @property
    def _idx_dict_global(self) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """
        Indexing helper for selecting cells across a dictionary of anndata instances.

        Increasing indices across data sets which can be concatenated into a single index vector with unique entries
        for cells.

        E.g.: For two data sets of 10 cells each, the return value would be

            {A: ([0..9], [0..9]), B: ([0..9], [10..19])}.
        """
        counter = 0
        indices = {}
        for k, v in self.adata_dict.items():
            indices[k] = (np.arange(0, v.n_obs), np.arange(counter, counter + v.n_obs))
            counter += len(v)
        return indices

    def _obs_idx_dict_query(self, idx: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Dictionary of indices of selected observations by data set key.

        :param idx: Array (global indexing) of selected observations.
        """
        if self.single_object:
            idx_dict = dict([(k, np.sort(idx)) for k in self.adata_dict.keys()])
        else:
            idx_set = set(idx)
            # Return data set-wise index if global index is in target set.
            idx_dict = dict([
                (k, np.sort([x1 for x1, x2 in zip(v1, v2) if x2 in idx_set]))
                for k, (v1, v2) in self._idx_dict_global.items()
            ])
            # Only retain non-empty.
            idx_dict = dict([(k, v) for k, v in idx_dict.items() if len(v) > 0])
        return idx_dict

    def _parse_array(self, x, return_dense):
        # Move from ArrayView to numpy if backed and dense:
        if (
                isinstance(x, anndata._core.views.ArrayView) or
                isinstance(x, anndata._core.views.SparseCSRView) or
                isinstance(x, anndata._core.views.SparseCSCView)
        ):
            x = x.toarray()
        if return_dense:
            x = np.asarray(x.todense()) if isinstance(x, scipy.sparse.spmatrix) else x
        elif isinstance(x, np.ndarray):
            x = scipy.sparse.csr_matrix(x) if isinstance(x, np.ndarray) else x
        if self.var_idx is not None:
            x = x[:, self.var_idx]
        return x


class CartDask(CartSingle):

    """
    Cart for a DistributedStoreDao().
    """

    _x: dask.array
    _obs: pd.DataFrame

    def __init__(self, x: dask.array, obs: pd.DataFrame, obs_keys: List[str], **kwargs):
        self._x = x
        self._obs = obs[obs_keys]
        super(CartDask, self).__init__(obs_keys=obs_keys, **kwargs)
        self.schedule.batchsplits = self._x.chunks[0]  # align batchsplits with partitions of the underlying dask array

    @property
    def _obs_full(self):
        """
        Full meta data matrix (cells x meta data) that is emitted in batches by .iterator().
        """
        return self._obs

    @property
    def _x_full(self):
        """
        Full data matrix (cells x features) that is emitted in batches by .iterator().
        """
        return self._x

    def _iterator(self):
        """
        Iterator over data matrix and meta data table, yields batches of data points.
        """
        # Can all data sets corresponding to one organism as a single array because they share the second dimension
        # and dask keeps expression data and obs out of memory.
        # Trigger design once so that its random elements are different every time this iterator is called (ie
        # ever epoch for example):
        for batch_idxs in self.schedule.design:
            if len(batch_idxs) > 0:
                x_i = self._x[batch_idxs, :]
                if self.var_idx is not None:
                    x_i = x_i[:, self.var_idx]
                obs_i = self._obs.iloc[batch_idxs, :]
                if self._emit_obsm:
                    obsm_i = dict([(k, v[batch_idxs, :]) for k, v in self.obsm.items()])
                    map_fn_args = (x_i, obs_i, obsm_i)
                else:
                    map_fn_args = (x_i, obs_i)
                data_tuple = self.map_fn(*map_fn_args)
                if self.batch_size == 1:
                    for data_tuple_i in split_batch(x=data_tuple):
                        yield data_tuple_i
                else:
                    yield data_tuple

    def move_to_memory(self):
        """
        Persist underlying dask array into memory in sparse.CSR format.
        """
        self._x = self._x.map_blocks(scipy.sparse.csr_matrix).persist()

    @property
    def adata(self) -> anndata.AnnData:
        """
        Assembles a slice of this cart based on .obs_idx as an anndata instance.

        TODO this can be removed once anndata can deal with dask array in x
        """
        return anndata.AnnData(X=scipy.sparse.csr_matrix(self.x.compute()), obs=self.obs, var=self.var)

    @property
    def n_obs(self) -> int:
        """Total number of observations in cart."""
        return self._x.shape[0]

    @property
    def n_var(self) -> int:
        """Total number of features defined for return in cart."""
        if self.var_idx is None:
            return self._x.shape[1]
        else:
            return len(self.var_idx)

    @property
    def obs(self):
        """
        Selected meta data matrix (cells x meta data) that is emitted in batches by .iterator().
        """
        return self._obs.iloc[self.schedule.idx, :]

    @property
    def x(self):
        """
        Selected data matrix (cells x features) that is emitted in batches by .iterator() as csr matrix.

        Note: this is a dask array.
        """
        if self.var_idx is None:
            return self._x[self.schedule.idx, :]
        else:
            return self._x[self.schedule.idx, :][:, self.var_idx]
