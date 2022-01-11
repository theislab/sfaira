"""Script to benchmark performance of sfaira data generators."""

import os
import sys
import time
import warnings
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
import sfaira

# Set global variables.
print(f'sys.argv: {sys.argv}')

N_DRAWS = 320_000
LEN_IDX = 5_000_000
BATCH_SIZE = 1  # must be 0 or 1
OBS_KEYS = ['cell_type', 'cell_line', 'organism', 'organ']  # list of obs_keys to retrieve
RETRIEVAL_BATCH_SIZE = 65536  # number of samples to retrieve at once

path_store_h5ad = str(sys.argv[1])
path_store_dao = str(sys.argv[2])
path_out = str(sys.argv[3])

store_type = ['dao']
if path_store_h5ad.lower() != 'none':
    store_type.append("h5ad")

time_measurements_initiate = {'storage_format': [], 'instantiation_time': [], 'run': []}
memory_measurements_initiate = {'storage_format': [], 'memory_usage': [], 'run': []}
time_measurements = {
    'scenario': [], 'storage_format': [], 'data_access_type': [], 'varsubset': [], 'avg_time_per_sample': []
}


def _map_fn(x_sample, obs_sample):
    gene_expression = np.asarray(x_sample)
    obs = tuple(obs_sample[obs_key].to_numpy().reshape((-1, 1)) for obs_key in OBS_KEYS)
    x = (gene_expression,)
    y = (gene_expression,) + obs

    return x, y


def _time_gen(_store: sfaira.data.store.DistributedStoreSingleFeatureSpace,
              store_format: str,
              kwargs_generator: Dict[str, any],
              num_draws: int) -> List[float]:
    if store_format == "h5ad":
        del kwargs_generator["random_access"]
    if kwargs_generator["var_subset"]:
        gc = sfaira.versions.genomes.genomes.GenomeContainer(organism='Homo sapiens', release='104')
        gc.set(symbols=["VTA1", "MLXIPL", "BAZ1B", "RANBP9", "PPARGC1A", "DDX25", "CRYAB"])
        _store.genome_container = gc
    del kwargs_generator["var_subset"]
    _gen = (
        _store
        .checkout(**kwargs_generator)
        .iterator()
    )
    _measurements = []
    for _ in range(num_draws):
        _t0 = time.perf_counter()
        _ = next(_gen)
        _measurements.append(time.perf_counter() - _t0)

    return _measurements


def _create_generator_kwargs(index: np.ndarray,
                             var_subset: bool,
                             random_batch_access: bool,
                             random_access: bool):

    if random_access and random_batch_access:
        raise ValueError('You cannot select "random_access" and "random_batch_access" at the same time')

    return {
        "idx": index,
        "batch_size": BATCH_SIZE,
        "retrieval_batch_size": RETRIEVAL_BATCH_SIZE,
        "map_fn": _map_fn,
        "obs_keys": OBS_KEYS,
        "randomized_batch_access": random_batch_access,
        "random_access": random_access,
        "var_subset": var_subset,
        "return_dense": True
    }


# check if data sets contain the same datasets
if path_store_h5ad.lower() != 'none':
    store = sfaira.data.load_store(cache_path=path_store_dao, store_format="dao").stores['Homo sapiens']
    data_set_lengths_dao = {dataset: len(idx_arr) for dataset, idx_arr in store.indices.items()}
    store = sfaira.data.load_store(cache_path=path_store_h5ad, store_format="h5ad").stores['Homo sapiens']
    data_set_lengths_h5ad = {dataset: len(idx_arr) for dataset, idx_arr in store.indices.items()}

    for dataset in set(list(data_set_lengths_dao.keys()) + list(data_set_lengths_h5ad.keys())):
        if dataset not in data_set_lengths_dao:
            warnings.warn(f'{dataset} dataset missing in dao storage')
            continue
        if dataset not in data_set_lengths_h5ad:
            warnings.warn(f'{dataset} dataset missing in h5ad storage')
            continue
        n_cells_dao = data_set_lengths_dao[dataset]
        n_cells_h5ad = data_set_lengths_h5ad[dataset]
        if n_cells_dao != n_cells_h5ad:
            warnings.warn(
                f'{dataset} dataset has different lengths in dao (n={n_cells_dao} cells) storage '
                f'and h5ad storage (n={n_cells_h5ad} cells)'
            )


for store_type_i in store_type:
    print(f'Benchmarking {store_type_i} storage')
    path_store = path_store_h5ad if store_type_i == "h5ad" else path_store_dao

    print('Benchmarking storage instantiation')
    for i in range(3):
        t0 = time.perf_counter()
        store = sfaira.data.load_store(cache_path=path_store, store_format=store_type_i)
        # Include initialisation of generator in timing to time overhead generated here.
        _ = store.checkout(map_fn=_map_fn, obs_keys=OBS_KEYS).iterator()
        time_measurements_initiate['instantiation_time'].append(time.perf_counter() - t0)
        time_measurements_initiate['storage_format'].append(store_type_i)
        time_measurements_initiate['run'].append(i)
        memory_measurements_initiate['memory_usage'].append(np.sum(list(store.adata_memory_footprint.values())))
        memory_measurements_initiate['storage_format'].append(store_type_i)
        memory_measurements_initiate['run'].append(i)

    # Prepare benchmark
    store = sfaira.data.load_store(cache_path=path_store, store_format=store_type_i).stores['Homo sapiens']
    if BATCH_SIZE == 1:
        n_draws = int(N_DRAWS)
    else:
        n_draws = int(N_DRAWS * RETRIEVAL_BATCH_SIZE)

    for scenario in ['seq_idx', 'random_idx']:
        print(f'Benchmarking scenario: {scenario}')

        if scenario == 'seq_idx':
            idx = np.arange(0, min(LEN_IDX, store.n_obs), dtype=int)
        elif scenario == 'random_idx':
            idx = np.arange(0, min(LEN_IDX, store.n_obs), dtype=int)
            np.random.shuffle(idx)
        else:
            raise ValueError(f'scenario={scenario} is not defined')

        for data_access_type in ['sequential', 'random-batch-access', 'random-access']:
            for varsubset in [False, True]:

                time_measurements['scenario'].append(scenario)
                time_measurements['storage_format'].append(store_type_i)
                time_measurements['data_access_type'].append(data_access_type)
                time_measurements['varsubset'].append(varsubset)
                if data_access_type == 'sequential':
                    random_batch_access_ = False
                    random_access_ = False
                elif data_access_type == 'random-batch-access':
                    random_batch_access_ = True
                    random_access_ = False
                elif data_access_type == 'random-access':
                    random_batch_access_ = False
                    random_access_ = True
                else:
                    raise ValueError(f'data_access_type={data_access_type} is not supported')
                kwargs = _create_generator_kwargs(idx, varsubset, random_batch_access_, random_access_)
                measurements = _time_gen(store, store_type_i, kwargs, num_draws=min(n_draws, len(idx)))
                time_measurements['avg_time_per_sample'].append(np.mean(measurements))


# prepare results
instatiation_time_df = pd.DataFrame(time_measurements_initiate)
memory_usage_df = pd.DataFrame(memory_measurements_initiate)
res_df = pd.DataFrame(time_measurements).assign(avg_time_per_sample=lambda xx: xx.avg_time_per_sample * 10**6)

# save results to csv
res_df.to_csv(os.path.join(path_out, 'data_store_benchmark.csv'))
instatiation_time_df.to_csv(os.path.join(path_out, 'instantiation_time_benchmark.csv'))
memory_usage_df.to_csv(os.path.join(path_out, 'memory_usage_benchmark.csv'))

# create figures
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12, 12))

axs[0, 0].set_title('Storage instantiation time')
sb.barplot(
    x='storage_format', y='instantiation_time', data=instatiation_time_df, ax=axs[0, 0]
)
axs[0, 0].set_ylabel('time [s]')
axs[0, 0].set_yscale('log')

axs[0, 1].set_title('Storage memory footprint')
sb.barplot(x='storage_format', y='memory_usage', data=memory_usage_df, ax=axs[0, 1])
axs[0, 1].set_ylabel('memory usage [MB]')

axs[1, 0].set_title('Avg. time per sample [μs] | seq_idx & varsubset=False')
sb.barplot(
    x='storage_format',
    y='avg_time_per_sample',
    hue='data_access_type',
    data=res_df[res_df.varsubset.eq(False) & res_df.scenario.eq('seq_idx')],
    ax=axs[1, 0]
)
axs[1, 0].set_ylabel('avg. time [μs]')
axs[1, 0].set_yscale('log')

axs[1, 1].set_title('Avg. time per sample [μs] | seq_idx & varsubset=True')
sb.barplot(
    x='storage_format',
    y='avg_time_per_sample',
    hue='data_access_type',
    data=res_df[res_df.varsubset.eq(True) & res_df.scenario.eq('seq_idx')],
    ax=axs[1, 1]
)
axs[1, 1].set_ylabel('avg. time [μs]')
axs[1, 1].set_yscale('log')

axs[2, 0].set_title('Avg. time per sample [μs] | random_idx & varsubset=False')
sb.barplot(
    x='storage_format',
    y='avg_time_per_sample',
    hue='data_access_type',
    data=res_df[res_df.varsubset.eq(False) & res_df.scenario.eq('random_idx')],
    ax=axs[2, 0]
)
axs[2, 0].set_ylabel('avg. time [μs]')
axs[2, 0].set_yscale('log')

axs[2, 1].set_title('Avg. time per sample [μs] | random_idx & varsubset=True')
sb.barplot(
    x='storage_format',
    y='avg_time_per_sample',
    hue='data_access_type',
    data=res_df[res_df.varsubset.eq(True) & res_df.scenario.eq('random_idx')],
    ax=axs[2, 1]
)
axs[2, 1].set_ylabel('avg. time [μs]')
axs[2, 1].set_yscale('log')

# set y-scale to same range
y_lims = []
for i in range(1, 3):
    for j in range(0, 2):
        y_lims.append(axs[i, j].get_ylim()[1])
for i in range(1, 3):
    for j in range(0, 2):
        axs[i, j].set_ylim(1, max(y_lims))

plt.tight_layout()
plt.savefig(os.path.join(path_out, 'data_store_benchmark.pdf'))
