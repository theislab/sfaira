import numpy as np
from typing import List, Tuple


def _get_batch_start_ends(idx: np.ndarray, batch_size: int):
    n_obs = len(idx)
    remainder = n_obs % batch_size if n_obs > 0 else 0
    n_batches = int(n_obs // batch_size + int(remainder > 0)) if n_obs > 0 else 0
    batch_starts_ends = [
        (int(x * batch_size), int(np.minimum((x * batch_size) + batch_size, n_obs)))
        for x in np.arange(0, n_batches)
    ]
    return batch_starts_ends


def _randomize_batch_start_ends(batch_starts_ends):
    batch_range = np.arange(0, len(batch_starts_ends))
    np.random.shuffle(batch_range)
    batch_starts_ends = [batch_starts_ends[i] for i in batch_range]
    return batch_starts_ends


class BatchDesignBase:

    def __init__(self, retrieval_batch_size: int, randomized_batch_access: bool, random_access: bool, **kwargs):
        self.retrieval_batch_size = retrieval_batch_size
        self._idx = None
        if randomized_batch_access and random_access:
            raise ValueError("Do not use randomized_batch_access and random_access.")
        self.randomized_batch_access = randomized_batch_access
        self.random_access = random_access

    @property
    def batch_bounds(self):
        """
        Protects property from changing.
        """
        return self._batch_bounds

    @property
    def idx(self):
        """
        Protects property from uncontrolled changing.
        Changes to _idx require changes to _batch_bounds.
        """
        return self._idx

    @idx.setter
    def idx(self, x):
        self._batch_bounds = _get_batch_start_ends(idx=x, batch_size=self.retrieval_batch_size)
        self._idx = np.sort(x)  # Sorted indices improve accession efficiency in some cases.

    @property
    def design(self) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
        """
        Yields index objects for one epoch of all data.

        These index objects are used by generators that have access to the data objects to build data batches.
        Randomization is performed anew with every call to this property.

        :returns: Tuple of:
            - Ordering of observations in epoch.
            - Batch start and end indices for batch based on ordering defined in first output.
        """
        raise NotImplementedError()


class BatchDesignBasic(BatchDesignBase):

    @property
    def design(self) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
        idx_proc = self.idx.copy()
        if self.random_access:
            np.random.shuffle(idx_proc)
        batch_bounds = self.batch_bounds.copy()
        if self.randomized_batch_access:
            batch_bounds = _randomize_batch_start_ends(batch_starts_ends=batch_bounds)
        return idx_proc, batch_bounds


class BatchDesignBalanced(BatchDesignBase):

    def __init__(self, grouping, group_weights: dict, randomized_batch_access: bool, random_access: bool,
                 **kwargs):
        """
        :param grouping: Group label for each entry in idx.
        :param group_weights: Group weight for each unique group in grouping. Does not have to normalise to a probability
            distribution but is normalised in this function. The outcome vector is always of length idx.
        """
        super(BatchDesignBalanced, self).__init__(randomized_batch_access=randomized_batch_access,
                                                  random_access=random_access, **kwargs)
        if randomized_batch_access:
            print("WARNING: randomized_batch_access==True is not a meaningful setting for BatchDesignBalanced.")
        if not random_access:
            print("WARNING: random_access==False is dangerous if you do not work with a large shuffle buffer "
                  "downstream of the sfaira generator.")
        # Create integer group assignment array.
        groups = np.sort(list(group_weights.keys()))
        grouping_int = np.zeros((grouping.shape[0],), dtype="int32") - 1
        for i, x in enumerate(groups):
            grouping_int[np.where(grouping == x)[0]] = i
        assert np.all(grouping_int >= 0)
        # Create sampling weights: Sampling weights are a probability distribution over groups.
        weight_group = np.array([group_weights[x] for x in groups])
        p_obs = np.asarray(weight_group[grouping_int], dtype="float64")
        p_obs = p_obs / np.sum(p_obs)
        if np.any(p_obs == 0.):
            raise ValueError(f"Down-sampling resulted in zero-probability weights on cells. "
                             f"Group weights: {weight_group}")
        self.p_obs = p_obs

    @property
    def design(self) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
        # Re-sample index vector.
        idx_proc = np.random.choice(a=self.idx, replace=True, size=len(self.idx), p=self.p_obs)
        if not self.random_access:  # Note: randomization is result from sampling above, need to revert if not desired.
            idx_proc = np.sort(idx_proc)
        batch_bounds = self.batch_bounds.copy()
        if self.randomized_batch_access:
            batch_bounds = _randomize_batch_start_ends(batch_starts_ends=batch_bounds)
        return idx_proc, batch_bounds


BATCH_SCHEDULE = {
    "base": BatchDesignBasic,
    "balanced": BatchDesignBalanced,
}
