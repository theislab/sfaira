from functools import partial
from typing import Dict, Union

import anndata
import numpy as np

from sfaira.data.store.carts.base import CartBase
from sfaira.data.store.carts.dataset_utils import _ShuffleBuffer, _DatasetIteratorRepeater
from sfaira.data.store.carts.single import CartSingle


class CartMulti(CartBase):

    """
    Cart for a DistributedStoreMultipleFeatureSpaceBase().
    """

    carts: Dict[str, CartSingle]
    intercalated: bool
    map_fn_merge: Union[None, callable]

    def __init__(self, carts: Dict[str, CartSingle], intercalated: bool = False, map_fn_merge=None):
        """

        :param carts: The generators to combine.
        :param intercalated: Whether to intercalate batches from both generators. Intercalates at frequency such that
            all generators terminate roughly at the same time.
            time.
        """
        self.carts = carts
        self.intercalated = intercalated
        self.map_fn_merge = map_fn_merge
        self._iterator_frequencies = None

    def adaptor(
            self,
            generator_type: str,
            dataset_kwargs: dict = None,
            shuffle_buffer: int = 0,
            shuffle_buffer_multi: int = 0,
            repeat: int = 1,
            **kwargs
    ):
        """
        See documentation of self._adaptor().
        """
        iter_kwargs = {"shuffle_buffer": shuffle_buffer, "shuffle_buffer_multi": shuffle_buffer_multi, "repeat": repeat}
        return self._adaptor(generator_type=generator_type, dataset_kwargs=dataset_kwargs, iter_kwargs=iter_kwargs,
                             **kwargs)

    def adaptor_torch(self, dataset_kwargs, loader, **kwargs):
        from torch.utils.data import DataLoader
        # Only import this module if torch is used to avoid strict torch dependency:
        from sfaira.data.store.torch_dataset import SfairaDataset, SfairaMultiDataset

        g = dict([(k, SfairaDataset(map_fn=c.map_fn_base, obs=c.obs, obsm=c.obsm, x=c.x, **dataset_kwargs))
                  for k, c in self.carts.items() if len(c.schedule.idx) > 0])
        # Wrap data set into a multi dataset that emits a list of batch tensors with one element for each dataset.
        # Note that in contrast to self.iterator which is used in iteratable-style access,
        # this list of tensor lists needs to be decomposed by map_fn to yield batch-dimensional tensors.
        g = SfairaMultiDataset(datasets=g, map_fn_merge=self.map_fn_merge)
        if loader:
            g = DataLoader(g, **kwargs)
        return g

    def adaptor_torch_iter(self, loader, repeat, shuffle_buffer, shuffle_buffer_multi, batch_size=None,
                           interleaved: bool = False, **kwargs):
        if interleaved:
            # Note: shuffle_buffer_multi is ignored here as this data sream batches per cart, which means that
            # observations do not need to be shuffled across streams from individual carts.
            from torch.utils.data import DataLoader
            # Only import this module if torch is used to avoid strict torch dependency:
            from sfaira.data.store.torch_dataset import InterleavedIterableDataset
            assert batch_size is not None, "supply batch_size for interleaved multi loader"

            datasets = [v.adaptor_torch_iter(loader=False, repeat=repeat, shuffle_buffer=shuffle_buffer)
                        for v in self.carts.values()]
            weights = [v.n_obs_selected for v in self.carts.values()]
            # Note on batch handling: the batch size is set in the dataset, the loader accepts the batches from the
            # data, this is enforced by giving batch_size=None to the loader.
            g = InterleavedIterableDataset(datasets=datasets, weights=weights, batch_size=batch_size)
            if loader:
                g = DataLoader(g, batch_size=None, **kwargs)
        else:
            from torch.utils.data import DataLoader
            # Only import this module if torch is used to avoid strict torch dependency:
            from sfaira.data.store.torch_dataset import SfairaIterableDataset

            g = SfairaIterableDataset(iterator_fun=partial(self.iterator, repeat=repeat, shuffle_buffer=shuffle_buffer,
                                                           shuffle_buffer_multi=shuffle_buffer_multi))
            if loader:
                g = DataLoader(g, **kwargs)
        return g

    @property
    def adata(self) -> Dict[str, anndata.AnnData]:
        """
        Assembles a dictionary of slices of this cart based on .obs_idx as anndata instances per organism.
        """
        return dict([(k, v.adata) for k, v in self.carts.items()])

    def iterator(self, repeat: int = 1, shuffle_buffer: int = 0, shuffle_buffer_multi: int = 0):
        """
        Iterator over data matrix and meta data table, yields batches of data points.

        Note: uses same shuffle buffer size across organisms and within organism, these are separate buffers though!
        """
        keep_repeating = True
        num_repetitions = 0
        if self.intercalated:

            #def _iterator():
            iterator_frequencies = self.iterator_frequencies.tolist().copy()
            iterators = [v.iterator(repeat=1, shuffle_buffer=shuffle_buffer) for v in self.carts.values()]
            while keep_repeating:
                while len(iterators) > 0:
                    # Sample iterator with frequencies so that in expectation, the frequency of samples from each
                    # iterator is uniform over an epoch.
                    itertor_idx = np.where(np.random.multinomial(n=1, pvals=iterator_frequencies))[0][0]
                    try:
                        x = next(iterators[itertor_idx])
                        yield x
                    except StopIteration:
                        # Remove iterator from list to sample. Once all are removed, the loop terminates.
                        del iterators[itertor_idx]
                        del iterator_frequencies[itertor_idx]
                    num_repetitions += 1
                    keep_repeating = (num_repetitions < repeat) or (repeat <= 0)

        else:
            if not repeat == 1 or len(self.carts.keys()) == 1:
                raise ValueError("using non-intercalated iterator with more than one cart and multiple repeats,"
                                 "likely not desired.")

            #def _iterator():
            while keep_repeating:
                for gi in self.carts.values():
                    for x in gi.iterator(): # repeat=repeat, shuffle_buffer=shuffle_buffer):
                        yield x
                num_repetitions += 1
                keep_repeating = (num_repetitions < repeat) or (repeat <= 0)

        #if shuffle_buffer_multi > 2 and np.all([x.batch_size == 1 for x in self.carts.values()]):
        #    iterator = _ShuffleBuffer(generator=_iterator, buffer_size=shuffle_buffer_multi).iterator
        #else:
        #    iterator = _iterator
        ## Note: do not repeat this overall iterator as the individual carts are already repeated.
        #return _DatasetIteratorRepeater(iterator, n_repeats=1).iterator()

    @property
    def n_batches(self) -> int:
        return int(np.sum([v.n_batches for v in self.carts.values()]))

    @property
    def n_obs(self) -> int:
        """Total number of observations in cart."""
        return int(np.sum([v.n_obs for v in self.carts.values()]))

    @property
    def n_obs_selected(self) -> int:
        """Total number of selected observations in cart."""
        return int(np.sum([v.n_obs_selected for v in self.carts.values()]))

    @property
    def n_var(self) -> Dict[str, int]:
        """Total number of features defined for return in cart."""
        return dict([(k, v.n_var) for k, v in self.carts.items()])

    @property
    def obs(self):
        return dict([(k, v.obs) for k, v in self.carts.items()])

    @property
    def obs_idx(self):
        """
        Dictionary of integer observation indices to select from cart. These will be emitted if data is queried.
        """
        return dict([(k, v.obs_idx) for k, v in self.carts.items()])

    @obs_idx.setter
    def obs_idx(self, x):
        if x is None:
            x = dict([(k, None) for k in self.carts.keys()])
        for k in self.carts.keys():
            assert k in x.keys(), (x.keys(), self.carts.keys())
            self.carts[k].obs_idx = x[k]
        self._iterator_frequencies = None  # Reset ratios.

    def move_to_memory(self):
        """
        Persist underlying arrays into memory in sparse.COO format.
        """
        for v in self.carts.values():
            v.move_to_memory()

    @property
    def iterator_frequencies(self):
        """
        Define relative drawing frequencies from iterators for intercalation.

        Note that these can be float and will be randomly rounded during intercalation.
        """
        if self._iterator_frequencies is None:
            iterator_lens = np.array([v.n_obs_selected for v in self.carts.values()])  # eg. [10, 15, 20]
            # Compute ratios of sampling so that one batch is drawn from the smallest generator per intercalation cycle.
            # See also self.iterator.
            freq = iterator_lens / np.sum(iterator_lens)  # eg. [10/45, 15/45, 20/45]
            self._iterator_frequencies = freq
        return self._iterator_frequencies

    @property
    def var(self):
        return dict([(k, v.var) for k, v in self.carts.items()])

    @property
    def x(self):
        return dict([(k, v.x) for k, v in self.carts.items()])
