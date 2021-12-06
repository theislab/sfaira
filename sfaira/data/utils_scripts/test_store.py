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

N_DRAWS = 100
BATCH_SIZE = 0
OBS_KEYS = ['cell_type', 'cell_line', 'organism', 'organ']
RETRIEVAL_BATCH_SIZE = 128
DAO_CHUNKSIZE = 128  # this is only relevant for DAO store

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

time_measurements_initiate = {store: [] for store in store_type}
memory_measurements_initiate = {store: [] for store in store_type}
time_measurements = {
    "load_sequential_from_one_dataset": {},
    "load_sequential_from_many_datasets": {},
    "load_random_from_one_dataset": {},
    "load_random_from_many_datasets": {},
    "load_sequential_from_one_dataset_todense": {},
    "load_sequential_from_many_datasets_todense": {},
    "load_random_from_one_dataset_todense": {},
    "load_random_from_many_datasets_todense": {},
    "load_sequential_from_one_dataset_todense_varsubet": {},
    "load_sequential_from_many_datasets_todense_varsubet": {},
    "load_random_from_one_dataset_todense_varsubet": {},
    "load_random_from_many_datasets_todense_varsubet": {}
}


def map_fn(x_sample, obs_sample):
    gene_expression = np.asarray(x_sample)
    obs = tuple(obs_sample[obs_key].to_numpy().reshape((-1, 1)) for obs_key in OBS_KEYS)
    x = (gene_expression, )
    y = (gene_expression, ) + obs

    return x, y


def time_gen(_store: DistributedStoreSingleFeatureSpace, store_format: str, kwargs_generator) -> List[float]:
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
    for _ in range(N_DRAWS):
        _t0 = time.perf_counter()
        _ = next(_gen)
        _measurements.append(time.perf_counter() - _t0)

    return _measurements


def get_idx_dataset_start(_store: DistributedStoreSingleFeatureSpace, k_target):
    idx_ = {}
    counter = 0
    for k, v in _store.indices.items():
        if k in k_target:
            idx_[k] = counter
        counter += len(v)
    return [idx_[k] for k in k_target]


# Define data objects to be comparable:
store = sfaira.data.load_store(cache_path=path_store_dao, store_format="dao")
store.subset(attr_key="organism", values="Homo sapiens")
store = store.stores["Homo sapiens"]
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
    for _ in range(3):
        t0 = time.perf_counter()
        store = sfaira.data.load_store(cache_path=path_store, store_format=store_type_i)
        # Include initialisation of generator in timing to time overhead generated here.
        _ = store.generator(map_fn=map_fn, obs_keys=OBS_KEYS).iterator()
        time_measurements_initiate[store_type_i].append(time.perf_counter() - t0)
        memory_measurements_initiate[store_type_i].append(np.sum(list(store.adata_memory_footprint.values())))

    # Prepare benchmark
    store = sfaira.data.load_store(cache_path=path_store, store_format=store_type_i)
    store.subset(attr_key="organism", values="Homo sapiens")
    store = store.stores["Homo sapiens"]
    idx_dataset_start = get_idx_dataset_start(_store=store, k_target=k_datasets)
    idx_dataset_end = [i + len(store.indices[x]) for i, x in zip(idx_dataset_start, k_datasets)]

    # Measure load_sequential_from_one_datasethttps://www.tensorflow.org/guide/data_performance_analysis?hl=en time.
    scenario = "load_sequential_from_one_dataset"
    print(f'Benchmarking scenario: {scenario}')
    for dense, varsubset in [(False, False), (True, False), (True, True)]:
        suffix = "_todense_varsubet" if dense and varsubset else "_todense" if dense and not varsubset else ""
        if BATCH_SIZE == 1:
            idx = np.arange(idx_dataset_start[0], idx_dataset_start[0] + N_DRAWS)
        else:
            idx = np.arange(idx_dataset_start[0], idx_dataset_start[0] + (N_DRAWS * RETRIEVAL_BATCH_SIZE))
        kwargs = {
            "idx": idx,
            "batch_size": BATCH_SIZE,
            "retrieval_batch_size": RETRIEVAL_BATCH_SIZE,
            "map_fn": map_fn,
            "obs_keys": OBS_KEYS,
            "return_dense": dense,
            "randomized_batch_access": False,
            "random_access": False,
            "var_subset": varsubset,
        }
        time_measurements[scenario + suffix][store_type_i] = time_gen(
            _store=store, store_format=store_type_i, kwargs_generator=kwargs
        )

    # Measure load_random_from_one_dataset time.
    scenario = "load_random_from_one_dataset"
    print(f'Benchmarking scenario: {scenario}')
    for dense, varsubset in [(False, False), (True, False), (True, True)]:
        suffix = "_todense_varsubet" if dense and varsubset else "_todense" if dense and not varsubset else ""
        if BATCH_SIZE == 1:
            idx = np.random.choice(
                np.arange(idx_dataset_start[0], np.maximum(idx_dataset_end[0], idx_dataset_start[0] + N_DRAWS)),
                size=N_DRAWS, replace=False
            )
        else:
            idx = np.random.choice(
                np.arange(
                    idx_dataset_start[0],
                    np.maximum(idx_dataset_end[0], idx_dataset_start[0] + (N_DRAWS * RETRIEVAL_BATCH_SIZE))
                ),
                size=(N_DRAWS * RETRIEVAL_BATCH_SIZE), replace=False
            )

        kwargs = {
            "idx": idx,
            "batch_size": BATCH_SIZE,
            "retrieval_batch_size": RETRIEVAL_BATCH_SIZE,
            "map_fn": map_fn,
            "obs_keys": OBS_KEYS,
            "return_dense": dense,
            "randomized_batch_access": False,
            "random_access": False,
            "var_subset": varsubset,
        }
        time_measurements[scenario + suffix][store_type_i] = time_gen(
            _store=store, store_format=store_type_i, kwargs_generator=kwargs
        )

    # Measure load_sequential_from_many_datasets time.
    scenario = "load_sequential_from_many_datasets"
    print(f'Benchmarking scenario: {scenario}')
    for dense, varsubset in [(False, False), (True, False), (True, True)]:
        suffix = "_todense_varsubet" if dense and varsubset else "_todense" if dense and not varsubset else ""
        if BATCH_SIZE == 1:
            idx = np.concatenate([np.arange(s, s + N_DRAWS) for s in idx_dataset_start])
        else:
            idx = np.concatenate([np.arange(s, s + (N_DRAWS * RETRIEVAL_BATCH_SIZE)) for s in idx_dataset_start])
        kwargs = {
            "idx": idx,
            "batch_size": BATCH_SIZE,
            "retrieval_batch_size": RETRIEVAL_BATCH_SIZE,
            "map_fn": map_fn,
            "obs_keys": OBS_KEYS,
            "return_dense": dense,
            "randomized_batch_access": False,
            "random_access": False,
            "var_subset": varsubset,
        }
        time_measurements[scenario + suffix][store_type_i] = time_gen(
            _store=store, store_format=store_type_i, kwargs_generator=kwargs
        )

    # Measure load_random_from_many_datasets time.
    scenario = "load_random_from_many_datasets"
    print(f'Benchmarking scenario: {scenario}')
    for dense, varsubset in [(False, False), (True, False), (True, True)]:
        suffix = "_todense_varsubet" if dense and varsubset else "_todense" if dense and not varsubset else ""
        if BATCH_SIZE == 1:
            idx = np.concatenate([
                np.random.choice(np.arange(s, np.maximum(e, s + N_DRAWS)), size=N_DRAWS, replace=False)
                for s, e in zip(idx_dataset_start, idx_dataset_end)
            ])
        else:
            idx = np.concatenate([
                np.random.choice(
                    np.arange(s, np.maximum(e, s + (N_DRAWS * RETRIEVAL_BATCH_SIZE))),
                    size=N_DRAWS * RETRIEVAL_BATCH_SIZE, replace=False
                )
                for s, e in zip(idx_dataset_start, idx_dataset_end)
            ])
        kwargs = {
            "idx": np.concatenate([
                np.random.choice(np.arange(s, np.maximum(e, s + N_DRAWS)), size=N_DRAWS, replace=False)
                for s, e in zip(idx_dataset_start, idx_dataset_end)
            ]),
            "batch_size": BATCH_SIZE,
            "retrieval_batch_size": RETRIEVAL_BATCH_SIZE,
            "map_fn": map_fn,
            "obs_keys": OBS_KEYS,
            "return_dense": dense,
            "randomized_batch_access": False,
            "random_access": False,
            "var_subset": varsubset,
        }
        time_measurements[scenario + suffix][store_type_i] = time_gen(
            _store=store, store_format=store_type_i, kwargs_generator=kwargs
        )


# create figures
ncols = 2
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(14, 12))
for i, x in enumerate([
    [
        "initialisation time",
    ],
    [
        "initialisation memory",
    ],
    [
        "load_sequential_from_one_dataset",
        "load_sequential_from_one_dataset_todense",
        "load_sequential_from_one_dataset_todense_varsubet",
    ],
    [
        "load_sequential_from_many_datasets",
        "load_sequential_from_many_datasets_todense",
        "load_sequential_from_many_datasets_todense_varsubet",
    ],
    [
        "load_random_from_one_dataset",
        "load_random_from_one_dataset_todense",
        "load_random_from_one_dataset_todense_varsubet",
    ],
    [
        "load_random_from_many_datasets",
        "load_random_from_many_datasets_todense",
        "load_random_from_many_datasets_todense_varsubet",
    ],
]):
    if i == 0 or i == 1:
        if i == 0:
            measurements_initiate = time_measurements_initiate
            ylabel = "log10 time sec"
            log = True
        else:
            measurements_initiate = memory_measurements_initiate
            ylabel = "memory MB"
            log = False
        df_sb = pd.concat([
            pd.DataFrame({
                ylabel: np.log(measurements_initiate[m]) / np.log(10) if log else measurements_initiate[m],
                "store": m,
                "draw": range(len(measurements_initiate[m])),
            })
            for m in measurements_initiate.keys()
        ], axis=0)
        sb.boxplot(
            data=df_sb,
            x="store", y=ylabel,
            ax=axs[i // ncols, i % ncols]
        )
        axs[i // ncols, i % ncols].set_title(x)
    elif len(x) > 0:
        df_sb = pd.concat([
            pd.concat([
                pd.DataFrame({
                    "log10 time sec": np.log(time_measurements[m][n]) / np.log(10),
                    "scenario": " ".join(m.split("_")[4:]),
                    "store": n,
                    "draw": range(len(time_measurements[m][n])),
                })
                for n in time_measurements[m].keys()
            ], axis=0)
            for m in x
        ], axis=0)
        # Could collapse draws to mean and put batch size on x.
        sb.lineplot(
            data=df_sb,
            x="draw", y="log10 time sec", hue="scenario", style="store",
            ax=axs[i // ncols, i % ncols]
        )
        axs[i // ncols, i % ncols].set_title(x[0])

plt.tight_layout()
plt.savefig(os.path.join(path_out, "data_store_benchmark.pdf"))
