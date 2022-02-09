from random import shuffle
from typing import List, Tuple

import numpy as np
from numba import njit


def _split_idx_along_arr_chunks(indicies: np.ndarray, chunks: Tuple[int]):
    assert np.all(np.diff(indicies) >= 0.)  # make sure array is sorted
    assert indicies.dtype == np.int64

    @njit
    def split_along_array_boundaries(idxs: np.ndarray, arr_boundaries: np.ndarray):
        idx_split = []  # variable to accumulate index groups

        partition_ctr = 0
        start_next_partition = arr_boundaries[partition_ctr]
        idxs_partition = [np.int64(x) for x in range(0)]  # temporary variable to accumulate idxs for a index group
        for idx in idxs:
            if idx < start_next_partition:
                # append idx to index group
                idxs_partition.append(idx)
            else:
                # find next bigger partition
                while idx >= start_next_partition:
                    partition_ctr += 1
                    start_next_partition = arr_boundaries[partition_ctr]
                # start a new index group list
                if len(idxs_partition) > 0:
                    idx_split.append(idxs_partition)
                    idxs_partition = [np.int64(x) for x in range(0)]
                # append idx to index group
                idxs_partition.append(idx)
        if len(idxs_partition) > 0:
            idx_split.append(idxs_partition)

        return idx_split

    return [np.array(xx) for xx in split_along_array_boundaries(indicies, np.cumsum(chunks))]


class BatchDesignBase:

    """
    Manages distribution of selected indices for a given data object over subsequent batches.

    This involves randomisation and possible meta-data dependent batching.
    This class is centred on the property `.design` which yields a list of observation indices (arrays), where each
    array is the set of indices that map to a batch and the sequence of batches in the list encodes the sequence of
    batches during an epoch over the data.
    """

    def __init__(self,
                 retrieval_batch_size: int,
                 randomized_batch_access: bool,
                 random_access: bool,
                 **kwargs):
        self.retrieval_batch_size = retrieval_batch_size
        self._batches = None
        self._idx = None
        self._batch_splits = None
        if randomized_batch_access and random_access:
            raise ValueError("Do not use randomized_batch_access and random_access.")
        self.randomized_batch_access = randomized_batch_access
        self.random_access = random_access

    @property
    def batchsplits(self) -> Tuple[int]:
        """
        Batch size of the yieleded batches.
        """
        return self._batch_splits

    @batchsplits.setter
    def batchsplits(self, batch_splits: Tuple[int]):
        self._batch_splits = batch_splits

    @property
    def n_batches(self) -> int:
        return len(self.design)

    @property
    def idx(self):
        """
        Protects property from uncontrolled changing.
        Changes to _idx require changes to _batch_bounds.
        """
        return self._idx

    @idx.setter
    def idx(self, x):
        self._idx = np.sort(x)  # idx has to be sorted for logic in subclasses

    @property
    def design(self) -> List[np.array]:
        """
        Yields index objects for one epoch of all data.

        These index objects are used by generators that have access to the data objects to build data batches.
        Randomization is performed anew with every call to this property.

        :returns: List[np.array]
            List of indices per batch
        """
        raise NotImplementedError()


class BatchDesignBasic(BatchDesignBase):

    """Standard batched access to data."""

    @property
    def design(self) -> List[np.ndarray]:

        if self.random_access:
            # shuffle idx for random access
            idx = np.random.permutation(self.idx)
            batches = np.array_split(idx, max(len(idx) // self.retrieval_batch_size, 1))
        elif self.randomized_batch_access:
            if self.batchsplits is not None:
                # if zarr chunk sizes are known -> split first along those and shuffle those batches
                batches = _split_idx_along_arr_chunks(self.idx, self.batchsplits)
                shuffle(batches)
                # accumulate smaller batches to the size of retrieval_batch_size
                batches = np.array_split(np.concatenate(batches), max(len(self.idx) // self.retrieval_batch_size, 1))
            else:
                batches = np.array_split(self.idx, max(len(self.idx) // self.retrieval_batch_size, 1))
                shuffle(batches)
        else:
            batches = np.array_split(self.idx, max(len(self.idx) // self.retrieval_batch_size, 1))

        return batches


class BatchDesignBalanced(BatchDesignBase):

    """
    Balanced batches across meta data partitions of data.
    """

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
            print("WARNING: randomized_batch_access==True is not a meaningful setting for BatchDesignBalanced. "
                  "Setting will be ignored!p")
        if not random_access:
            print("WARNING: random_access==False is dangerous if you do not work with a large shuffle buffer "
                  "downstream of the sfaira generator.")
        # Create integer group assignment array.
        groups = np.sort(list(group_weights.keys()))
        grouping_int = np.zeros((grouping.shape[0], ), dtype="int32") - 1
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
    def design(self) -> List[np.ndarray]:
        # select relevant probabilities and renormalize to prob vector
        p_obs = self.p_obs[self.idx]
        p_obs = p_obs / np.sum(p_obs)
        # Re-sample index vector.
        idx = np.random.choice(a=self.idx, replace=True, size=len(self.idx), p=p_obs)
        if not self.random_access:  # Note: randomization is result from sampling above, need to revert if not desired.
            idx = np.sort(idx)
        batches = np.array_split(idx, max(len(idx) // self.retrieval_batch_size, 1))

        return batches


class BatchDesignBlocks(BatchDesignBase):

    """
    Meta data-defined blocks of observations in each batch.

    Note that the batches do not necessarily contain equal numbers of observations! The yielded batch size depends
    solely on the number of observations in a meta data partition.
    """

    def __init__(self, grouping, random_access: bool, **kwargs):
        """

        :param grouping: List of group labels for each element in idx. Must have same length as data array.
        :param group_weights: Group weight for each unique group in grouping. Does not have to normalise to a probability
            distribution but is normalised in this function. The outcome vector is always of length idx.
        """
        super(BatchDesignBlocks, self).__init__(random_access=random_access, **kwargs)
        if not random_access:
            print("WARNING: random_access==False is dangerous if you do not work with a large shuffle buffer "
                  "downstream of the sfaira generator.")
        # Create integer group assignment array.
        self.grouping = grouping

    @property
    def grouping(self):
        return self._grouping

    @grouping.setter
    def grouping(self, x):
        self._grouping = x
        # Reset:
        self._groups = None
        self.idx_sorted = None

    @property
    def groups(self):
        if self._groups is None:
            self._groups = np.unique(self.grouping)
        return self._groups

    @property
    def idx(self):
        """
        Protects property from uncontrolled changing.
        Changes to _idx require changes to batch splitting.
        """
        return self._idx

    @idx.setter
    def idx(self, x):
        self._idx = x
        # Reset:
        self._groups = None
        idx_sorted = []
        for group in self.groups:
            # subselect by group and only keep indicies which are in the supplied index range
            idx = np.intersect1d(np.where(self.grouping == group)[0], x)
            idx_sorted.append(idx)
        self.idx_sorted = idx_sorted

    @property
    def design(self) -> List[np.ndarray]:
        idx_sorted = self.idx_sorted
        batches = []
        for idxs_group in idx_sorted:
            if len(idxs_group) > 0:
                if self.random_access:
                    # shuffle subgroups if random_access
                    idxs_group = np.random.permutation(idx_sorted)
                batches.append(idxs_group)
        if self.random_access:
            # shuffle blocks
            shuffle(batches)

        return batches


class BatchDesignFull(BatchDesignBase):

    """Emits full dataset as a single batch in each query."""

    @property
    def design(self) -> List[np.ndarray]:
        idx = np.arange(0, len(self.idx))
        if self.random_access:
            # shuffle idx for random access
            idx = np.random.permutation(idx)
        return [idx]


BATCH_SCHEDULE = {
    "base": BatchDesignBasic,
    "balanced": BatchDesignBalanced,
    "blocks": BatchDesignBlocks,
    "full": BatchDesignFull,
}
