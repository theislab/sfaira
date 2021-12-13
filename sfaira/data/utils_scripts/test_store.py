import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sb
import sfaira
import sys
import time
from typing import List
from sfaira.data.store import DistributedStoreSingleFeatureSpace

# Set global variables.
print(f'sys.argv: {sys.argv}')

N_DRAWS = 100_000
BATCH_SIZE = 1  # must be 0 or 1
OBS_KEYS = ['cell_type', 'cell_line', 'organism', 'organ']  # list of obs_keys to retrieve
RETRIEVAL_BATCH_SIZE = 4096  # number of samples to retrieve at once
DAO_CHUNKSIZE = 128  # this has to match the chunksize of the dao storage

path_store_h5ad = str(sys.argv[1])
path_store_dao = str(sys.argv[2])
path_out = str(sys.argv[3])

store_type = []
kwargs = []
compression_kwargs = []
store_type.append("dao")
kwargs.append({"dense": True, "chunks": DAO_CHUNKSIZE})
compression_kwargs.append({"compressor": "default", "overwrite": True, "order": "C"})
if path_store_h5ad.lower() != 'none':
    store_type.append("h5ad")
    kwargs.append({"dense": False})
    compression_kwargs.append({})

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


def time_gen(_store: DistributedStoreSingleFeatureSpace, 
             store_format: str, 
             kwargs_generator,
             num_draws: int) -> List[float]:
    """
    Take samples from generator and measure time taken to generate each sample.
    """
    if store_format == "h5ad":
        del kwargs_generator["random_access"]
    if kwargs_generator["var_subset"]:
        gc = sfaira.versions.genomes.genomes.GenomeContainer(organism='Homo sapiens', release='104')
        gc.set(symbols=["VTA1", "MLXIPL", "BAZ1B", "RANBP9", "PPARGC1A", "DDX25", "CRYAB"])
        _store.genome_container = gc
    del kwargs_generator["var_subset"]
    _gen = (
        _store
        .generator(**kwargs_generator)
        .iterator()
    )
    _measurements = []
    for _ in range(num_draws):
        _t0 = time.perf_counter()
        _ = next(_gen)
        _measurements.append(time.perf_counter() - _t0)

    return _measurements


def create_generator_kwargs(index: np.ndarray,
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


def get_idx_dataset_start(_store: DistributedStoreSingleFeatureSpace, k_target):
    idx_ = {}
    counter = 0
    for k, v in _store.indices.items():
        if k in k_target:
            idx_[k] = counter
        counter += len(v)
    return [idx_[k] for k in k_target]


# Define data objects to be comparable:
store = sfaira.data.load_store(cache_path=path_store_dao, store_format="dao").stores['Homo sapiens']
k_datasets_dao = list(store.indices.keys())
k_datasets_dao = np.asarray(k_datasets_dao)[np.argsort([len(v) for v in store.indices.values()])].tolist()
if path_store_h5ad.lower() != 'none':
    store = sfaira.data.load_store(cache_path=path_store_h5ad, store_format="h5ad")
    store.subset(attr_key="organism", values="Homo sapiens")
    store = store.stores["Homo sapiens"]
    k_datasets_h5ad = list(store.indices.keys())
    # Only retain intersection of data sets while keeping order.
    k_datasets = [x for x in k_datasets_dao if x in k_datasets_h5ad]
else:
    k_datasets = k_datasets_dao
n_datasets = len(k_datasets)


print(f"Running benchmark on {n_datasets} data sets.")
for store_type_i, kwargs_i, compression_kwargs_i in zip(store_type, kwargs, compression_kwargs):
    path_store = path_store_h5ad if store_type_i == "h5ad" else path_store_dao
    print(f'Benchmarking {store_type_i} storage')

    print('Benchmarking storage instantiation')
    for i in range(3):
        t0 = time.perf_counter()
        store = sfaira.data.load_store(cache_path=path_store, store_format=store_type_i)
        # Include initialisation of generator in timing to time overhead generated here.
        _ = store.generator(map_fn=_map_fn, obs_keys=OBS_KEYS).iterator()
        time_measurements_initiate['instantiation_time'].append(time.perf_counter() - t0)
        time_measurements_initiate['storage_format'].append(store_type_i)
        time_measurements_initiate['run'].append(i)
        memory_measurements_initiate['memory_usage'].append(np.sum(list(store.adata_memory_footprint.values())))
        memory_measurements_initiate['storage_format'].append(store_type_i)
        memory_measurements_initiate['run'].append(i)

    # Prepare benchmark
    store = sfaira.data.load_store(cache_path=path_store, store_format=store_type_i).stores['Homo sapiens']
    idx_dataset_start = get_idx_dataset_start(_store=store, k_target=k_datasets)
    idx_dataset_end = [i + len(store.indices[x]) for i, x in zip(idx_dataset_start, k_datasets)]
    if BATCH_SIZE == 1:
        draws_per_dataset = int(N_DRAWS / n_datasets)
        n_draws = int(N_DRAWS)
    else:
        draws_per_dataset = int(N_DRAWS * RETRIEVAL_BATCH_SIZE / n_datasets)
        n_draws = int(N_DRAWS * RETRIEVAL_BATCH_SIZE)

    for scenario in ['seq_idx', 'random_idx']:
        print(f'Benchmarking scenario: {scenario}')

        if scenario == 'seq_idx':
            idx = np.concatenate(
                [np.arange(s, np.minimum(s + draws_per_dataset, e, dtype=int))
                 for s, e in zip(idx_dataset_start, idx_dataset_end)]
            )
        elif scenario == 'random_idx':
            idxs_per_dataset = [np.arange(s, np.minimum(s + draws_per_dataset, e), dtype=int)
                                for s, e in zip(idx_dataset_start, idx_dataset_end)]
            concatenated_idxs = np.concatenate(idxs_per_dataset)
            idx = np.random.choice(
                concatenated_idxs, size=min(n_draws, len(concatenated_idxs)), replace=False
            )
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
                kwargs = create_generator_kwargs(idx, varsubset, random_batch_access_, random_access_)
                measurements = time_gen(store, store_type_i, kwargs, num_draws=min(n_draws, len(idx)))
                time_measurements['avg_time_per_sample'].append(np.mean(measurements))

    print()


# prepare results
instatiation_time_df = pd.DataFrame(time_measurements_initiate)
memory_usage_df = pd.DataFrame(memory_measurements_initiate)
res_df = pd.DataFrame(time_measurements).assign(
    access_type=lambda xx: xx.scenario + '_' + xx.data_access_type,
    avg_time_per_sample=lambda xx: xx.avg_time_per_sample * 10**3
)

# save results to csv
res_df.to_csv(os.path.join(path_out, 'data_store_benchmark.csv'))

# create figures
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))

axs[0, 0].set_title('Storage instantiation time')
sb.barplot(
    x='storage_format', y='instantiation_time', data=instatiation_time_df, ax=axs[0, 0]
)
axs[0, 0].set_ylabel('time [s]')
axs[0, 0].set_yscale('log')

axs[0, 1].set_title('Storage memory footprint')
sb.barplot(x='storage_format', y='memory_usage', data=memory_usage_df, ax=axs[0, 1])
axs[0, 1].set_ylabel('memory usage [MB]')

axs[1, 0].set_title('Avg. time per sample [ms] | varsubset=False')
sb.barplot(
    x='storage_format',
    y='avg_time_per_sample',
    hue='access_type',
    data=res_df[res_df.varsubset == False],
    ax=axs[1, 0]
)
axs[1, 0].set_ylabel('avg. time [ms]')
axs[1, 0].set_yscale('log')

axs[1, 1].set_title('Avg. time per sample [ms] | varsubset=True')
sb.barplot(
    x='storage_format',
    y='avg_time_per_sample',
    hue='access_type',
    data=res_df[res_df.varsubset == True],
    ax=axs[1, 1]
)
axs[1, 1].set_ylabel('avg. time [ms]')
axs[1, 1].set_yscale('log')

plt.tight_layout()
plt.savefig(os.path.join(path_out, "data_store_benchmark.pdf"))
