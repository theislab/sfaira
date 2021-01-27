import sfaira
import sys
import tensorflow as tf

print(tf.__version__)


def write_meta(args0, args1):
    # Write meta data, cache.
    args0.write_meta(fn_meta=None, dir_out=args1, fn_data=None)
    # Test load from cache.
    args0.load(
        celltype_version=None,
        match_to_reference=True,
        remove_gene_version=True,
        load_raw=False,
        allow_caching=False,
    )
    return None


# Set global variables.
print("sys.argv", sys.argv)

path = str(sys.argv[1])
path_meta = str(sys.argv[2])
path_cache = str(sys.argv[3])
processes = int(str(sys.argv[4]))

ds = sfaira.data.dataloaders.DatasetSuperGroupSfaira(
    path=path, meta_path=path_meta, cache_path=path_cache
)
dsg = ds.flatten()  # need to flatten in this case to parallelise across Groups and not just within.
# Write meta data, cache and test load from cache:
dsg.load(
    annotated_only=False,
    match_to_reference=None,
    remove_gene_version=True,
    load_raw=True,
    allow_caching=True,
    processes=processes,
    func=write_meta,
    kwargs_func={"args1": path_meta},
)
