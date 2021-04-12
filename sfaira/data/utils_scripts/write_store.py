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
    ds.streamline(
        format="sfaira",
        allow_uns_sfaira=True,
        clean_obs=True,
        clean_var=True,
        clean_uns=True,
    )
    ds.subset_genes(subset_type="protein_coding")
    ds.write_distributed_store(dir_cache=path_store, store="h5ad")
    ds.clear()
