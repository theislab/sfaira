import os
import sfaira
import sys

# Set global variables.
print("sys.argv", sys.argv)

data_path = str(sys.argv[1])
path_meta = str(sys.argv[2])
path_cache = str(sys.argv[3])
processes = int(str(sys.argv[4]))

ds = sfaira.data.dataloaders.Universe(
    data_path=data_path, meta_path=path_meta, cache_path=path_cache
)
# Write meta data, cache and test load from cache:
for x in ds.dataset_groups:
    for k, v in x.datasets.items():
        print(f"SCRIPT: loading {x} {k}")
        try:
            # Initial load and cache writing:
            # Only run this if data set was not already cached to speed up resumed jobs.
            if not os.path.exists(v.cache_fn):
                v.load(load_raw=False, allow_caching=True)
            # Write meta data, cache.
            v.write_meta(fn_meta=None, dir_out=path_meta)
            # Test load from cache.
            v.load(load_raw=False, allow_caching=False)
            v.clear()
        except ValueError as e:
            # Do not abort upon ValueErrors, such as from cell type map bugs.
            print(f"SCRIPT WARNING: TO-FIX: ValueError in {k}: {e}")
