import anndata
import dask
import numpy as np
import pandas as pd
import scipy.sparse
import time
from typing import Dict, List, Tuple


def split_batch(x, obs):
    """
    Splits retrieval batch into consumption batches of length 1.

    Often, end-user consumption batches would be observation-wise, ie yield a first dimension of length 1.
    """
    for i in range(x.shape[0]):
        yield x[[i], :], obs.iloc[[i], :]


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

    @property
    def iterator(self) -> iter:
        raise NotImplementedError()

    @property
    def n_batches(self) -> int:
        raise NotImplementedError()

    def adaptor(self, generator_type: str, **kwargs):
        """
        The adaptor turns a python base generator into a different iteratable object, defined by generator_type.

        :param generator_type: Type of output iteratable.
            - python base generator (no change to `.generator`)
            - tensorflow dataset
            - pytorch dataset
        :returns: Modified iteratable (see generator_type).
        """
        if generator_type == "python":
            g = self.iterator
        elif generator_type == "tensorflow":
            import tensorflow as tf
            g = tf.data.Dataset.from_generator(generator=self.iterator, **kwargs)
        else:
            raise ValueError(f"{generator_type} not recognized")
        return g


class GeneratorSingle(GeneratorBase):

    batch_start_ends: List[Tuple[int, int]]
    obs_idx: np.ndarray
    obs_keys: List[str]
    var_idx: np.ndarray

    def __init__(self, obs_idx, batch_start_ends, var_idx, obs_keys, map_fn):
        """
        :param obs_idx: np.ndarray: The cells to emit.
        :parm batch_starts_ends: List[Tuple[int, int]: Batch start and end indices.
        :param var_idx: The features to emit.
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adata_by_key.
        :param map_fn: Map function to apply to output tuple of raw generator. Each draw i from the generator is then:
            `yield map_fn(x[i, var_idx], obs[i, obs_keys])`
        """
        self.obs_idx = obs_idx
        self.batch_starts_ends = batch_start_ends
        self.map_fn = map_fn
        self.obs_keys = obs_keys
        self.var_idx = var_idx

    @property
    def n_batches(self) -> int:
        return len(self.batch_starts_ends)


class GeneratorAnndata(GeneratorSingle):

    adata_dict: anndata._core.views.ArrayView
    return_dense: bool
    single_object: bool

    def __init__(self, adata_dict, idx_dict_global, return_dense, obs_idx, batch_start_ends, var_idx, obs_keys, map_fn):
        super(GeneratorAnndata, self).__init__(obs_idx=obs_idx, batch_start_ends=batch_start_ends, var_idx=var_idx,
                                               obs_keys=obs_keys, map_fn=map_fn)
        self.return_dense = return_dense
        self.single_object = len(adata_dict.keys()) == 1
        self.idx_dict_global = idx_dict_global
        t0 = time.time()
        self.adata_dict = adata_dict
        print(f"init adata: {time.time() - t0}")

    @property
    def iterator(self) -> iter:
        # Speed up access to single object by skipping index overlap operations:

        def g():
            for s, e in self.batch_starts_ends:
                idx_i = self.obs_idx[s:e]
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
                # I) Prepare data matrix.
                x = [
                    self.adata_dict[k].X[v, :]
                    for k, v in idx_i_dict.items()
                ]
                # Move from ArrayView to numpy if backed and dense:
                x = [
                    xx.toarray() if isinstance(xx, anndata._core.views.ArrayView) else xx
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
                for x_i, obs_i in split_batch(x=x, obs=obs):
                    if self.map_fn is None:
                        yield x_i, obs_i
                    else:
                        output = self.map_fn(x_i, obs_i)
                        if output is not None:
                            yield output

        return g


class GeneratorDask(GeneratorSingle):

    x: dask.array
    obs: pd.DataFrame

    def __init__(self, x, obs, obs_idx, batch_start_ends, var_idx, obs_keys, map_fn):
        super(GeneratorDask, self).__init__(obs_idx=np.sort(obs_idx), batch_start_ends=batch_start_ends,
                                            var_idx=var_idx, obs_keys=obs_keys, map_fn=map_fn)
        t0 = time.time()
        self.x = x
        print(f"init X: {time.time() - t0}")
        t0 = time.time()
        self.obs = obs[self.obs_keys]
        # Redefine index so that .loc indexing can be used instead of .iloc indexing:
        obs.index = np.arange(0, obs.shape[0])
        print(f"init obs: {time.time() - t0}")

    @property
    def iterator(self) -> iter:
        # Can all data sets corresponding to one organism as a single array because they share the second dimension
        # and dask keeps expression data and obs out of memory.

        def g():
            x_temp = self.x[self.obs_idx, :]
            obs_temp = self.obs.loc[self.obs.index[self.obs_idx], :]  # TODO better than iloc?
            for s, e in self.batch_starts_ends:
                x_i = x_temp[s:e, :]
                if self.var_idx is not None:
                    x_i = x_i[:, self.var_idx]
                # Exploit fact that index of obs is just increasing list of integers, so we can use the .loc[]
                # indexing instead of .iloc[]:
                obs_i = obs_temp.loc[obs_temp.index[s:e], :]
                for x_ii, obs_ii in split_batch(x=x_i, obs=obs_i):
                    if self.map_fn is None:
                        yield x_ii, obs_ii
                    else:
                        output = self.map_fn(x_ii, obs_ii)
                        if output is not None:
                            yield output

        return g


class GeneratorMulti(GeneratorBase):

    generators: Dict[str, GeneratorSingle]
    intercalated: bool

    def __init__(self, generators: Dict[str, GeneratorSingle], intercalated: bool = False):
        self.generators = generators
        self.intercalated = intercalated
        # Define relative drawing frequencies from iterators for intercalation.
        gen_lens = np.array([v.n_batches for v in self.generators.values()])
        self.ratios = np.asarray(np.round(np.max(gen_lens) / np.asarray(gen_lens), 0), dtype="int64")

    @property
    def iterator(self) -> iter:

        if self.intercalated:
            def g():
                # Document which generators are still yielding batches:
                yielding = np.ones((self.ratios.shape[0],)) == 1.
                iterators = [v.iterator() for v in self.generators.values()]
                while np.any(yielding):
                    # Loop over one iterator length adjusted cycle of emissions.
                    for i, (g, n) in enumerate(zip(iterators, self.ratios)):
                        for _ in range(n):
                            try:
                                x = next(g)
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
