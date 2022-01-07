import abc

import anndata
import dask.array
import numpy as np
import pandas as pd
import scipy.sparse
from typing import Dict, List, Union

from sfaira.data.store.batch_schedule import BATCH_SCHEDULE, BatchDesignBase


def split_batch(x):
    """
    Splits retrieval batch into consumption batches of length 1.

    Often, end-user consumption batches would be observation-wise, ie yield a first dimension of length 1.

    :param x: Data tuple of length 1 or 2: (input,) or (input, output,), where both input and output are also
         a tuple, but of batch-dimensioned tensors.
    """
    batch_dim = x[0][0].shape[0]
    for i in range(batch_dim):
        output = []
        for y in x:
            if isinstance(y, tuple):
                output.append(tuple([z[i, :] for z in y]))
            else:
                output.append(y[i, :])
        yield tuple(output)


class GeneratorBase:
    """
    A generator is a shallow class that is instantiated on a pointer to a data set in a Store instance.

    This class exposes an iterator generator through `.iterator`.
    The iterator can often directly be used without a class like this one around it.
    However, that often implies that the namespace with the pointer to the data set is destroyed after iterator
    function declaration, which means that the pointer needs to be redefined for every full pass over the iterator.
    The class around this property maintains the namespace that includes the pointer and its instance can be used to
    avoid redefining the pointer every time the generator runs out and is re-called.
    For this run time advantage to not be compromised by memory load and class initialisation run time cost induced by
    actually copying data objects, it is important that the data object stored in this class is indeed a pointer.
    This is the case for:

        - lazily loaded dask arrays
        - anndata.Anndata view

    which have their own classes below.
    """

    schedule: BatchDesignBase

    @property
    def iterator(self) -> iter:
        raise NotImplementedError()

    @property
    def obs_idx(self):
        raise NotImplementedError()

    @property
    def n_batches(self) -> int:
        return self.schedule.n_batches

    def adaptor(self, generator_type: str, dataset_kwargs: dict = {}, **kwargs):
        """
        The adaptor turns a python base generator into a different iteratable object, defined by generator_type.

        :param generator_type: Type of output iteratable.
            - python base generator (no change to `.generator`)
            - tensorflow dataset: This dataset is defined on a python iterator.
            - pytorch: We distinguish torch.data.Dataset and torch.data.DataLoader ontop of either.
                The Dataset vs DataLoader distinction is made by the "" suffix for Dataset or "-loader" suffix for +
                dataloader. The distinction between Dataset and IteratableDataset defines if the object is defined
                directly on a dask array or based on a python iterator on a dask array. Note that the python iterator
                can implement favorable remote access schemata but the torch.data.Dataset generally causes less trouble
                in out-of-the-box usage.
                    - torch.data.Dataset: "torch" prefix, ie "torch" or "torch-loader"
                    - torch.data.IteratableDataset: "torch-iter" prefix, ie "torch-iter" or "torch-iter-loader"
        :param dataset_kwargs:
        :returns: Modified iteratable (see generator_type).
        """
        if generator_type == "python":
            g = self.iterator()
        elif generator_type == "tensorflow":
            import tensorflow as tf

            g = tf.data.Dataset.from_generator(generator=self.iterator, **kwargs)
        elif generator_type in ["torch", "torch-loader"]:
            import torch
            # Only import this module if torch is used to avoid strict torch dependency:
            from .torch_dataset import SfairaDataset

            g = SfairaDataset(map_fn=self.map_fn, obs=self.obs.iloc[self.obs_idx, :], x=self.x[self.obs_idx, :], **dataset_kwargs)
            if generator_type == "torch-loader":
                g = torch.utils.data.DataLoader(g, **kwargs)
        elif generator_type in ["torch-iter", "torch-iter-loader"]:
            import torch
            # Only import this module if torch is used to avoid strict torch dependency:
            from .torch_dataset import SfairaIterableDataset

            g = SfairaIterableDataset(iterator_fun=self.iterator)
            if generator_type == "torch-iter-loader":
                g = torch.utils.data.DataLoader(g, **kwargs)
        else:
            raise ValueError(f"{generator_type} not recognized")
        return g

    @property
    def x(self):
        pass

    @property
    def obs(self):
        pass


class GeneratorSingle(GeneratorBase):

    batch_size: int
    _obs_idx: Union[np.ndarray, None]
    obs_keys: List[str]
    var_idx: np.ndarray

    def __init__(self, batch_schedule, batch_size, map_fn, obs_idx, obs_keys, var_idx, **kwargs):
        """

        :param batch_schedule: str or class.
            - "basic"
            - "balanced"
            - class: batch_schedule needs to be a class (not instance), subclassing BatchDesignBase.
        :param batch_size: Emission batch size. Must be 1.
        :param map_fn: Map function to apply to output tuple of raw generator. Each draw i from the generator is then:
            `yield map_fn(x[i, var_idx], obs[i, obs_keys])`
        :param obs_idx: np.ndarray: The cells to emit.
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adata_by_key.
        :param var_idx: The features to emit.
        :parm split_to_obs: Whether to split tensors to observation-wise slices at the emission stage of the generator.
        """
        self.var_idx = var_idx
        self._obs_idx = None
        self.batch_schedule = batch_schedule
        self.batch_size = batch_size
        self.map_fn = map_fn
        if isinstance(batch_schedule, str):
            batch_schedule = BATCH_SCHEDULE[batch_schedule]
        self.schedule = batch_schedule(**kwargs)
        self.obs_idx = obs_idx
        self.obs_keys = obs_keys

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
                assert idx.dtype == np.int
            else:
                assert isinstance(idx, list)
                assert isinstance(idx[0], int) or isinstance(idx[0], np.int)
        idx = np.asarray(idx)
        return idx

    @property
    def obs_idx(self):
        return self._obs_idx

    @obs_idx.setter
    def obs_idx(self, x):
        """Allows emission of different iterator on same generator instance (using same dask array)."""
        if x is None:
            x = np.arange(0, self.n_obs)
        else:
            x = self._validate_idx(x)
            x = np.sort(x)
        # Only reset if they are actually different:
        if self._obs_idx is not None and len(x) != len(self._obs_idx) or np.any(x != self._obs_idx):
            self._obs_idx = x
            self.schedule.idx = x

    @property
    def n_obs(self) -> int:
        raise NotImplementedError()


class GeneratorAnndata(GeneratorSingle):

    adata_dict: Dict[str, anndata._core.views.ArrayView]
    return_dense: bool
    single_object: bool

    def __init__(self, adata_dict, idx_dict_global, return_dense, **kwargs):
        self.return_dense = return_dense
        self.single_object = len(adata_dict.keys()) == 1
        self.idx_dict_global = idx_dict_global
        self.adata_dict = adata_dict
        super(GeneratorAnndata, self).__init__(**kwargs)

    @property
    def n_obs(self) -> int:
        return int(np.sum([v.n_obs for v in self.adata_dict.values()]))

    @property
    def iterator(self) -> iter:
        # Speed up access to single object by skipping index overlap operations:

        def g():
            _, obs_idx, batch_bounds = self.schedule.design
            for s, e in batch_bounds:
                idx_i = obs_idx[s:e]
                # Match adata objects that overlap to batch:
                if self.single_object:
                    idx_i_dict = dict([(k, np.sort(idx_i)) for k in self.adata_dict.keys()])
                else:
                    idx_i_set = set(idx_i)
                    # Return data set-wise index if global index is in target set.
                    idx_i_dict = dict([
                        (k, np.sort([x2 for x1, x2 in zip(v1, v2) if x1 in idx_i_set]))
                        for k, (v1, v2) in self.idx_dict_global.items()
                    ])
                    # Only retain non-empty.
                    idx_i_dict = dict([(k, v) for k, v in idx_i_dict.items() if len(v) > 0])
                if self.batch_size == 1:
                    # Emit each data set separately and avoid concatenation into larger chunks for emission.
                    for k, v in idx_i_dict.items():
                        # I) Prepare data matrix.
                        x = self.adata_dict[k].X[v, :]
                        # Move from ArrayView to numpy if backed and dense:
                        if isinstance(x, anndata._core.views.ArrayView):
                            x = x.toarray()
                        if isinstance(x, anndata._core.views.SparseCSRView) or \
                                isinstance(x, anndata._core.views.SparseCSCView):
                            x = x.toarray()
                        # Do dense conversion now so that col-wise indexing is not slow, often, dense conversion
                        # would be done later anyway.
                        if self.return_dense:
                            x = np.asarray(x.todense()) if isinstance(x, scipy.sparse.spmatrix) else x
                        if self.var_idx is not None:
                            x = x[:, self.var_idx]
                        # Prepare .obs.
                        obs = self.adata_dict[k].obs[self.obs_keys].iloc[v, :]
                        data_tuple = self.map_fn(x, obs)
                        for data_tuple_i in split_batch(x=data_tuple):
                            yield data_tuple_i
                else:
                    # Concatenates slices first before returning. Note that this is likely slower than emitting by
                    # observation in most scenarios.
                    # I) Prepare data matrix.
                    x = [
                        self.adata_dict[k].X[v, :]
                        for k, v in idx_i_dict.items()
                    ]
                    # Move from ArrayView to numpy if backed and dense:
                    x = [
                        xx.toarray()
                        if (isinstance(xx, anndata._core.views.ArrayView) or
                            isinstance(xx, anndata._core.views.SparseCSRView) or
                            isinstance(xx, anndata._core.views.SparseCSCView))
                        else xx
                        for xx in x
                    ]
                    # Do dense conversion now so that col-wise indexing is not slow, often, dense conversion
                    # would be done later anyway.
                    if self.return_dense:
                        x = [np.asarray(xx.todense()) if isinstance(xx, scipy.sparse.spmatrix) else xx for xx in x]
                        is_dense = True
                    else:
                        is_dense = isinstance(x[0], np.ndarray)
                    # Concatenate blocks in observation dimension:
                    if len(x) > 1:
                        if is_dense:
                            x = np.concatenate(x, axis=0)
                        else:
                            x = scipy.sparse.vstack(x)
                    else:
                        x = x[0]
                    if self.var_idx is not None:
                        x = x[:, self.var_idx]
                    # Prepare .obs.
                    obs = pd.concat([
                        self.adata_dict[k].obs[self.obs_keys].iloc[v, :]
                        for k, v in idx_i_dict.items()
                    ], axis=0, join="inner", ignore_index=True, copy=False)
                    data_tuple = self.map_fn(x, obs)
                    yield data_tuple

        return g

    @property
    def obs_idx_dict(self):
        return dict([(k, v2) for k, (_, v2) in self.idx_dict_global.items()])

    @property
    def x(self):
        # Assumes that .X are scipy.sparse.csr or coo
        if self.var_idx is None:
            return scipy.sparse.vstack([
                self.adata_dict[k].X[v, :].copy()
                for k, v in self.obs_idx_dict.items()
            ])
        else:
            return scipy.sparse.vstack([
                self.adata_dict[k].X[v, :][:, self.var_idx].copy()
                for k, v in self.obs_idx_dict.items()
            ])

    @property
    def obs(self):
        return pd.concat([
            self.adata_dict[k].obs[self.obs_keys].iloc[v, :]
            for k, v in self.obs_idx_dict.items()
        ], axis=0, join="inner", ignore_index=True, copy=False)


class GeneratorDask(GeneratorSingle):

    """
    In addition to the full data array, x, this class maintains a slice _x_slice which is indexed by the iterator.
    Access to the slice can be optimised with dask and is therefore desirable.
    """

    _x: dask.array.Array
    _x_slice: Union[dask.array.Array, None]
    _obs: pd.DataFrame
    _obs_slice: Union[pd.DataFrame, None]

    def __init__(self, x, obs, obs_keys, var_idx, **kwargs):
        if var_idx is not None:
            x = x[:, var_idx]
        self._x = x
        self._obs = obs[obs_keys]
        # Redefine index so that .loc indexing can be used instead of .iloc indexing:
        self._obs.index = np.arange(0, obs.shape[0])
        self._x_slice = None
        self._obs_slice = None
        super(GeneratorDask, self).__init__(obs_keys=obs_keys, var_idx=var_idx, **kwargs)

    @property
    def n_obs(self) -> int:
        return self._x.shape[0]

    @property
    def obs_idx(self):
        return self._obs_idx

    @obs_idx.setter
    def obs_idx(self, x):
        """
        Allows emission of different iterator on same generator instance (using same dask array).
        In addition to base method: allows for optimisation of dask array for batch draws.
        """
        if x is None:
            x = np.arange(0, self.n_obs)
        else:
            x = self._validate_idx(x)
            x = np.sort(x)
        # Only reset if they are actually different:
        if (self._obs_idx is not None and len(x) != len(self._obs_idx)) or np.any(x != self._obs_idx):
            self._obs_idx = x
            self.schedule.idx = x
            self._x_slice = dask.optimize(self._x[self._obs_idx, :])[0]
            self._obs_slice = self._obs.loc[self._obs.index[self._obs_idx], :]  # TODO better than iloc?
            # Redefine index so that .loc indexing can be used instead of .iloc indexing:
            self._obs_slice.index = np.arange(0, self._obs_slice.shape[0])

    @property
    def iterator(self) -> iter:
        # Can all data sets corresponding to one organism as a single array because they share the second dimension
        # and dask keeps expression data and obs out of memory.

        def g():
            obs_idx_slice, _, batch_bounds = self.schedule.design
            x_temp = self._x_slice
            obs_temp = self._obs_slice
            for s, e in batch_bounds:
                x_i = x_temp[obs_idx_slice[s:e], :]
                # Exploit fact that index of obs is just increasing list of integers, so we can use the .loc[]
                # indexing instead of .iloc[]:
                obs_i = obs_temp.loc[obs_temp.index[obs_idx_slice[s:e]], :]
                data_tuple = self.map_fn(x_i, obs_i)
                if self.batch_size == 1:
                    for data_tuple_i in split_batch(x=data_tuple):
                        yield data_tuple_i
                else:
                    yield data_tuple

        return g

    @property
    def x(self):
        return self._x_slice

    @property
    def obs(self):
        return self._obs_slice


class GeneratorMulti(GeneratorBase):

    """
    Generator that wraps multiple other generators.
    """

    generators: Dict[str, GeneratorSingle]
    intercalated: bool

    def __init__(self, generators: Dict[str, GeneratorSingle], intercalated: bool = False):
        """

        :param generators: The generators to combine.
        :param intercalated: Whether to intercalate batches from both generators. Intercalates at frequency such that
            all generators terminate roughly at the same time.
            time.
        """
        self.generators = generators
        self.intercalated = intercalated
        self._ratios = None

    @property
    def ratios(self):
        """
        Define relative drawing frequencies from iterators for intercalation.

        Note that these can be float and will be randomly rounded during intercalation.
        """
        if self._ratios is None:
            gen_lens = np.array([v.n_batches for v in self.generators.values()])  # eg. [10, 15, 20]
            # Compute ratios of sampling so that one batch is drawn from the smallest generator per intercalation cycle.
            # See also self.iterator.
            freq = gen_lens / np.min(gen_lens)  # eg. [1.0, 1.5, 2.0]
            self._ratios = freq
        return self._ratios

    @property
    def obs_idx(self):
        return dict([(k, v.obs_idx) for k, v in self.generators.items()])

    @obs_idx.setter
    def obs_idx(self, x):
        """Allows emission of different iterator on same generator instance (using same dask array)."""
        if x is None:
            x = dict([(k, None) for k in self.generators.keys()])
        for k in self.generators.keys():
            assert k in x.keys(), (x.keys(), self.generators.keys())
            self.generators[k].obs_idx = x[k]
        self._ratios = None  # Reset ratios.

    @property
    def iterator(self) -> iter:

        if self.intercalated:
            ratios = self.ratios.copy()
            print(f"GENERATOR: intercalating generators at ratios {ratios}")

            def g():
                # Document which generators are still yielding batches:
                yielding = np.ones((ratios.shape[0],)) == 1.
                iterators = [v.iterator() for v in self.generators.values()]
                while np.any(yielding):
                    # Loop over one iterator length adjusted cycle of emissions.
                    for i, (gi, n) in enumerate(zip(iterators, ratios)):
                        n_rounded = int(n) + np.random.binomial(n=1, p=n - int(n))
                        for _ in range(n_rounded):
                            try:
                                x = next(gi)
                                yield x
                            except StopIteration:
                                yielding[i] = False
        else:
            def g():
                for gi in self.generators.values():
                    for x in gi.iterator():
                        yield x

        return g

    @property
    def n_batches(self) -> int:
        return np.sum([v.n_batches for v in self.generators.values()])

    @property
    def x(self):
        return dict([(k, v.x) for k, v in self.generators])

    @property
    def obs(self):
        return dict([(k, v.obs) for k, v in self.generators])
