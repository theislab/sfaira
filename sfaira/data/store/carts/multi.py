import anndata
from typing import Dict

import numpy as np
from sfaira.data.store.carts.base import CartBase
from sfaira.data.store.carts.single import CartSingle


class CartMulti(CartBase):

    """
    Cart for a DistributedStoreMultipleFeatureSpaceBase().
    """

    carts: Dict[str, CartSingle]
    intercalated: bool

    def __init__(self, carts: Dict[str, CartSingle], intercalated: bool = False):
        """

        :param carts: The generators to combine.
        :param intercalated: Whether to intercalate batches from both generators. Intercalates at frequency such that
            all generators terminate roughly at the same time.
            time.
        """
        self.carts = carts
        self.intercalated = intercalated
        self._ratios = None

    @property
    def adata(self) -> Dict[str, anndata.AnnData]:
        """
        Assembles a dictionary of slices of this cart based on .obs_idx as anndata instances per organism.
        """
        return dict([(k, v.adata) for k, v in self.carts.items()])

    @property
    def iterator(self) -> iter:

        if self.intercalated:
            ratios = self.ratios.copy()
            print(f"GENERATOR: intercalating generators at ratios {ratios}")

            def g():
                # Document which generators are still yielding batches:
                yielding = np.ones((ratios.shape[0],)) == 1.
                iterators = [v.iterator() for v in self.carts.values()]
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
                for gi in self.carts.values():
                    for x in gi.iterator():
                        yield x

        return g

    @property
    def n_batches(self) -> int:
        return np.sum([v.n_batches for v in self.carts.values()])

    @property
    def n_obs(self) -> int:
        """Total number of observations in cart."""
        return np.sum([v.n_obs for v in self.carts.values()])

    @property
    def n_obs_selected(self) -> int:
        """Total number of selected observations in cart."""
        return np.sum([v.n_obs_selected for v in self.carts.values()])

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
        self._ratios = None  # Reset ratios.

    @property
    def ratios(self):
        """
        Define relative drawing frequencies from iterators for intercalation.

        Note that these can be float and will be randomly rounded during intercalation.
        """
        if self._ratios is None:
            gen_lens = np.array([v.n_batches for v in self.carts.values()])  # eg. [10, 15, 20]
            # Compute ratios of sampling so that one batch is drawn from the smallest generator per intercalation cycle.
            # See also self.iterator.
            freq = gen_lens / np.min(gen_lens)  # eg. [1.0, 1.5, 2.0]
            self._ratios = freq
        return self._ratios

    @property
    def var(self):
        return dict([(k, v.var) for k, v in self.carts.items()])

    @property
    def x(self):
        return dict([(k, v.x) for k, v in self.carts.items()])
