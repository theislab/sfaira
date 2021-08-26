import numpy as np


def batch_designer_basic(idx, batch_size: int, randomized_batch_access: bool, random_access: bool):
    """
    Yields index objects for one epoch of all data.

    These index objects are used by generators that have access to the data objects to build data batches.

    :returns: Tuple of:
        - Ordering of observations in epoch.
        - Batch start and end indices for batch based on ordering defined in first output.
    """
    n_obs = len(idx)
    remainder = n_obs % batch_size
    n_batches = int(n_obs // batch_size + int(remainder > 0))
    batch_starts_ends = [
        (int(x * batch_size), int(np.minimum((x * batch_size) + batch_size, n_obs)))
        for x in np.arange(0, n_batches)
    ]
    batch_range = np.arange(0, len(batch_starts_ends))
    if randomized_batch_access:
        np.random.shuffle(batch_range)
    batch_starts_ends = [batch_starts_ends[i] for i in batch_range]
    idx_proc = idx.copy()
    if random_access:
        np.random.shuffle(idx_proc)
    return idx_proc, batch_starts_ends


def batch_designer_balanced(idx, batch_size: int, grouping, group_weights: dict, randomized_batch_access: bool, random_access: bool):
    batch_starts_ends = []
    idx_proc = []
    groups = np.unique(list(group_weights.keys()))
    weights = np.array([group_weights[x] for x in groups])
    p = weights / np.sum(weights)
    groups = np.random.choice(a=groups, replace=True, size=len(groups), p=p)
    for x in groups:
        idx_x = np.where(grouping == x)[0]
        n_obs = len(idx_x)
        batch_size_o = int(np.minimum(batch_size, n_obs))
        batch_starts_ends.append(np.array([(0, batch_size_o), ]))
        if balance_obs is None:
            p = np.ones_like(idx_x) / len(idx_x)
        else:
            if balance_obs not in self.obs.columns:
                raise ValueError(f"did not find column {balance_obs} in {self.organism}")
            val_meta_x = val_meta[idx_x]
            class_freq = dict([(y, np.mean(val_meta_x == y)) for y in np.unique(val_meta_x)])
            class_freq_x_by_obs = np.array([class_freq[y] for y in val_meta_x])
            damped_freq_coefficient = np.maximum(balance_damping, (1. - class_freq_x_by_obs))
            p = np.ones_like(idx_x) / len(idx_x) * damped_freq_coefficient
        idx_x_sample = np.random.choice(a=idx_x, replace=False, size=batch_size_o, p=p)
        idx_proc.append(idx_x_sample)
    idx_proc = np.asarray(idx_proc)
    return idx_proc, batch_starts_ends


BATCH_DESIGNS = {
    "base": batch_designer_basic,
    "balanced": batch_designer_balanced,
}
