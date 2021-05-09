import os
import sfaira
import sys

# Set global variables.
print("sys.argv", sys.argv)

data_path = str(sys.argv[1])
path_meta = str(sys.argv[2])
path_cache = str(sys.argv[3])
path_store = str(sys.argv[4])
store_type = str(sys.argv[5]).lower()

# On disk format hyperparameters
if store_type == "h5ad":
    # Write sparse arrays in h5ad.
    compression_kwargs = {"dense": False, "compression_kwargs": {}}
elif store_type == "zarr":
    # Write dense arrays in zarr.
    compression_kwargs = {"dense": True, "chunks": 64,
                          "compression_kwargs": {"dtype": "float32", "compressor": "default", "overwrite": True,
                                                 "order": "C"}}
else:
    assert False, store_type

universe = sfaira.data.dataloaders.Universe(data_path=data_path, meta_path=path_meta, cache_path=path_cache)

for k, ds in universe.datasets.items():
    fn_store = os.path.join(path_store, ds.doi_cleaned_id + ".h5ad")
    if os.path.exists(fn_store):
        print(f"SCRIPT skipping {k}")
    else:
        print(f"SCRIPT loading {k}")
        ds.load(
            load_raw=False,
            allow_caching=True,
        )
        ds.streamline_features(
            remove_gene_version=True,
            match_to_reference={"human": "Homo_sapiens.GRCh38.102", "mouse": "Mus_musculus.GRCm38.102"},
            subset_genes_to_type="protein_coding"
        )
        ds.streamline_metadata(schema="sfaira", uns_to_obs=True, clean_obs=True, clean_var=True, clean_uns=True,
                               clean_obs_names=True)
        ds.write_distributed_store(dir_cache=path_store, store=store_type, **compression_kwargs)
        ds.clear()
