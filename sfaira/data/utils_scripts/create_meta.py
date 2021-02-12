import sfaira
import sys
import tensorflow as tf

print(tf.__version__)


def write_meta(args0, args1):
    args0.write_meta(fn_meta=None, dir_out=args1, fn_data=None)
    return None


# Set global variables.
print("sys.argv", sys.argv)

data_path = str(sys.argv[1])
path_meta = str(sys.argv[2])
processes = int(str(sys.argv[3]))

ds = sfaira.data.dataloaders.DatasetSuperGroupSfaira(
    data_path=data_path, meta_path=path_meta, cache_path=path_meta
)
dsg = ds.flatten()  # need to flatten in this case to parallelise across Groups and not just within.
dsg.load(
    annotated_only=False,
    match_to_reference=None,
    remove_gene_version=True,
    load_raw=True,
    allow_caching=False,
    processes=processes,
    func=write_meta,
    kwargs_func={"args1": path_meta},
)
