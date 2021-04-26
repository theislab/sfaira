import sfaira
import sys

# Set global variables.
print("sys.argv", sys.argv)

data_path = str(sys.argv[1])
path_meta = str(sys.argv[2])
path_cache = str(sys.argv[3])
path_store = str(sys.argv[4])

universe = sfaira.data.dataloaders.Universe(data_path=data_path, meta_path=path_meta, cache_path=path_cache)

for k, ds in universe.datasets.items():
    print(f"SCRIPT loading {k}")
    ds.load(
        match_to_reference=None,
        remove_gene_version=True,
        load_raw=False,
        allow_caching=True,
        set_metadata=False,
    )
    ds.streamline_features(remove_gene_version=True, match_to_reference=True, subset_genes_to_type="protein_coding")
    ds.streamline_metadata(schema="sfaira", uns_to_obs=False, clean_obs=False, clean_var=True, clean_uns=False,
                           clean_obs_names=True)
    ds.write_distributed_store(dir_cache=path_store, store="h5ad")
    ds.clear()
