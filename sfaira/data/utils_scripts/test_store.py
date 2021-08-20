import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sb
import sfaira
import sys
import time
from typing import List

# Set global variables.
print("sys.argv", sys.argv)

N_DRAWS = 10
BATCH_SIZES = [64]

path_store_h5ad = str(sys.argv[1])
path_store_dao = str(sys.argv[2])
path_out = str(sys.argv[3])

store_type = []
kwargs = []
compression_kwargs = []
if path_store_h5ad.lower() != "none":
    store_type.append("h5ad")
    kwargs.append({"dense": False})
    compression_kwargs.append({})

store_type.append("dao")
kwargs.append({"dense": True, "chunks": 128})
compression_kwargs.append({"compressor": "default", "overwrite": True, "order": "C"})

time_measurements_initiate = {}
memory_measurements_initiate = {}
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
    "load_random_from_many_datasets_todense_varsubet": {},
}


def time_gen(_store, store_format, kwargs) -> List[float]:
    """
    Take samples from generator and measure time taken to generate each sample.
    """
    if store_format == "h5ad":
        del kwargs["random_access"]
    if kwargs["var_subset"]:
        gc = sfaira.versions.genomes.genomes.GenomeContainer(assembly="Homo_sapiens.GRCh38.102")
        gc.subset(symbols=["VTA1", "MLXIPL", "BAZ1B", "RANBP9", "PPARGC1A", "DDX25", "CRYAB"])
        _store.genome_containers = gc
    del kwargs["var_subset"]
    _gen = _store.generator(**kwargs)()
    _measurements = []
    for _ in range(N_DRAWS):
        _t0 = time.time()
        _ = next(_gen)
        _measurements.append(time.time() - _t0)
    return _measurements


def get_idx_dataset_start(_store, k_target):
    idx = {}
    counter = 0
    for k, v in _store.indices.items():
        if k in k_target:
            idx[k] = counter
        counter += len(v)
    return [idx[k] for k in k_target]


# Define data objects to be comparable:
store = sfaira.data.load_store(cache_path=path_store_dao, store_format="dao")
store.subset(attr_key="organism", values="human")
store = store.stores["human"]
k_datasets_dao = list(store.indices.keys())
# Sort by size:
k_datasets_dao = np.asarray(k_datasets_dao)[np.argsort([len(v) for v in store.indices.values()])].tolist()
store = sfaira.data.load_store(cache_path=path_store_h5ad, store_format="h5ad")
store.subset(attr_key="organism", values="human")
store = store.stores["human"]
k_datasets_h5ad = list(store.indices.keys())
# Only retain intersection of data sets while keeping order.
k_datasets = [x for x in k_datasets_dao if x in k_datasets_h5ad]
n_datasets = len(k_datasets)
print(f"running benchmark on {n_datasets} data sets.")
for store_type_i, kwargs_i, compression_kwargs_i in zip(store_type, kwargs, compression_kwargs):
    path_store = path_store_h5ad if store_type_i == "h5ad" else path_store_dao

    # Measure initiate time.
    time_measurements_initiate[store_type_i] = []
    memory_measurements_initiate[store_type_i] = []
    for _ in range(3):
        t0 = time.time()
        store = sfaira.data.load_store(cache_path=path_store, store_format=store_type_i)
        # Include initialisation of generator in timing to time overhead generated here.
        _ = store.generator()
        time_measurements_initiate[store_type_i].append(time.time() - t0)
        memory_measurements_initiate[store_type_i].append(np.sum(list(store.adata_memory_footprint.values())))

    time_measurements["load_sequential_from_one_dataset"][store_type_i] = {}
    time_measurements["load_sequential_from_many_datasets"][store_type_i] = {}
    time_measurements["load_random_from_one_dataset"][store_type_i] = {}
    time_measurements["load_random_from_many_datasets"][store_type_i] = {}
    time_measurements["load_sequential_from_one_dataset_todense"][store_type_i] = {}
    time_measurements["load_sequential_from_many_datasets_todense"][store_type_i] = {}
    time_measurements["load_random_from_one_dataset_todense"][store_type_i] = {}
    time_measurements["load_random_from_many_datasets_todense"][store_type_i] = {}
    time_measurements["load_sequential_from_one_dataset_todense_varsubet"][store_type_i] = {}
    time_measurements["load_sequential_from_many_datasets_todense_varsubet"][store_type_i] = {}
    time_measurements["load_random_from_one_dataset_todense_varsubet"][store_type_i] = {}
    time_measurements["load_random_from_many_datasets_todense_varsubet"][store_type_i] = {}
    store = sfaira.data.load_store(cache_path=path_store, store_format=store_type_i)
    store.subset(attr_key="organism", values="human")
    store = store.stores["human"]
    idx_dataset_start = get_idx_dataset_start(_store=store, k_target=k_datasets)
    idx_dataset_end = [i + len(store.indices[x]) for i, x in zip(idx_dataset_start, k_datasets)]
    for bs in BATCH_SIZES:
        key_bs = "bs" + str(bs)
        print(key_bs)

        # Measure load_sequential_from_one_dataset time.
        scenario = "load_sequential_from_one_dataset"
        print(scenario)
        for dense_varsubset in [(False, False), (True, False), (True, True)]:
            dense, varsubset = dense_varsubset
            suffix = "_todense_varsubet" if dense and varsubset else "_todense" if dense and not varsubset else ""
            kwargs = {
                "idx": np.concatenate([
                    np.arange(idx_dataset_start[0] + bs * i, idx_dataset_start[0] + bs * (i + 1))
                    for i in range(N_DRAWS)]),
                "batch_size": bs,
                "return_dense": dense,
                "randomized_batch_access": False,
                "random_access": False,
                "var_subset": varsubset,
            }
            time_measurements[scenario + suffix][store_type_i][key_bs] = time_gen(
                _store=store, store_format=store_type_i, kwargs=kwargs)

        # Measure load_random_from_one_dataset time.
        scenario = "load_random_from_one_dataset"
        print(scenario)
        for dense_varsubset in [(False, False), (True, False), (True, True)]:
            dense, varsubset = dense_varsubset
            suffix = "_todense_varsubet" if dense and varsubset else "_todense" if dense and not varsubset else ""
            kwargs = {
                "idx": np.random.choice(
                    np.arange(idx_dataset_start[0], np.maximum(idx_dataset_end[0], idx_dataset_start[0] + bs * N_DRAWS)),
                    size=bs * N_DRAWS, replace=False),
                "batch_size": bs,
                "return_dense": dense,
                "randomized_batch_access": False,
                "random_access": False,
                "var_subset": varsubset,
            }
            time_measurements[scenario + suffix][store_type_i][key_bs] = time_gen(
                _store=store, store_format=store_type_i, kwargs=kwargs)

        # Measure load_sequential_from_many_datasets time.
        scenario = "load_sequential_from_many_datasets"
        print(scenario)
        for dense_varsubset in [(False, False), (True, False), (True, True)]:
            dense, varsubset = dense_varsubset
            suffix = "_todense_varsubet" if dense and varsubset else "_todense" if dense and not varsubset else ""
            kwargs = {
                "idx": np.concatenate([np.arange(s, s + bs) for s in idx_dataset_start]),
                "batch_size": bs,
                "return_dense": dense,
                "randomized_batch_access": False,
                "random_access": False,
                "var_subset": varsubset,
            }
            time_measurements[scenario + suffix][store_type_i][key_bs] = time_gen(
                _store=store, store_format=store_type_i, kwargs=kwargs)

        # Measure load_random_from_many_datasets time.
        scenario = "load_random_from_many_datasets"
        print(scenario)
        for dense_varsubset in [(False, False), (True, False), (True, True)]:
            dense, varsubset = dense_varsubset
            suffix = "_todense_varsubet" if dense and varsubset else "_todense" if dense and not varsubset else ""
            kwargs = {
                "idx": np.concatenate([
                    np.random.choice(np.arange(s, np.maximum(e, s + bs)), size=bs, replace=False)
                    for s, e in zip(idx_dataset_start, idx_dataset_end)]),
                "batch_size": bs,
                "return_dense": dense,
                "randomized_batch_access": False,
                "random_access": True,
                "var_subset": varsubset,
            }
            time_measurements[scenario + suffix][store_type_i][key_bs] = time_gen(
                _store=store, store_format=store_type_i, kwargs=kwargs)

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
                pd.concat([
                    pd.DataFrame({
                        "log10 time sec": np.log(time_measurements[m][n][o]) / np.log(10),
                        "scenario": " ".join(m.split("_")[4:]),
                        "store": n,
                        "batch size": o,
                        "scenario - batch size": " ".join(m.split("_")[4:]) + "_bs" + str(o),
                        "draw": range(len(time_measurements[m][n][o])),
                    })
                    for o in time_measurements[m][n].keys()
                ], axis=0)
                for n in time_measurements[m].keys()
            ], axis=0)
            for m in x
        ], axis=0)
        # Could collapse draws to mean and put batch size on x.
        sb.lineplot(
            data=df_sb,
            x="draw", y="log10 time sec", hue="scenario - batch size", style="store",
            ax=axs[i // ncols, i % ncols]
        )
        axs[i // ncols, i % ncols].set_title(x[0])
plt.tight_layout()
plt.savefig(os.path.join(path_out, "data_store_benchmark.pdf"))
