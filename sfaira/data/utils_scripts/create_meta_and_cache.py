import sfaira
import sys
import tensorflow as tf

print(tf.__version__)


def write_meta(args0, args1):
    # Write meta data, cache.
    args0.write_meta(fn_meta=None, dir_out=args1)
    # Test load from cache.
    args0.load(
        remove_gene_version=True,
        load_raw=False,
        allow_caching=False,
    )
    return None


# Set global variables.
print("sys.argv", sys.argv)

data_path = str(sys.argv[1])
path_meta = str(sys.argv[2])
path_cache = str(sys.argv[3])
processes = int(str(sys.argv[4]))

ds = sfaira.data.dataloaders.DatasetSuperGroupSfaira(
    data_path=data_path, meta_path=path_meta, cache_path=path_cache
)
ds = ds.flatten()
# Write meta data, cache and test load from cache:
ds.load(
    annotated_only=False,
    match_to_reference=None,
    remove_gene_version=True,
    load_raw=True,
    allow_caching=True,
    processes=processes,
    func=write_meta,
    kwargs_func={"args1": path_meta},
)
