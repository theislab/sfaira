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

    def adaptor_torch_iter(self, loader, repeat, shuffle_buffer, batch_size=None, interleaved: bool = False, **kwargs):
        if interleaved:
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
            g = super().adaptor_torch_iter(loader=loader, repeat=repeat, shuffle_buffer=shuffle_buffer,
                                           batch_size=batch_size, **kwargs)
        return g

    @property
    def adata(self) -> Dict[str, anndata.AnnData]:
        """
        Assembles a dictionary of slices of this cart based on .obs_idx as anndata instances per organism.
        """
        return dict([(k, v.adata) for k, v in self.carts.items()])

    def iterator(self, repeat: int = 1, shuffle_buffer: int = 0):
        """
        Iterator over data matrix and meta data table, yields batches of data points.

        Note: uses same shuffle buffer size across organisms and within organism, these are separate buffers though!
        """

        def _iterator():
            if self.intercalated:
                iterator_frequencies = self.iterator_frequencies.tolist().copy()
                iterators = [v.iterator(repeat=1, shuffle_buffer=shuffle_buffer) for v in self.carts.values()]
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
            else:
                for gi in self.carts.values():
                    for x in gi.iterator():
                        yield x

        if shuffle_buffer > 2 and np.all([x.batch_size == 1 for x in self.carts.values()]):
            g_dataset = _ShuffleBuffer(_iterator, shuffle_buffer).iterator
        else:
            g_dataset = _iterator

        return _DatasetIteratorRepeater(g_dataset, n_repeats=repeat).iterator()

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
