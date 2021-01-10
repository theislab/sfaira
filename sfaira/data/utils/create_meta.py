import sfaira
import sys
import tensorflow as tf

print(tf.__version__)


def write_meta(args0, args1):
    args0.write_meta(dir_out=args1)
    return None


# Set global variables.
print("sys.argv", sys.argv)

path = str(sys.argv[1])
path_meta = str(sys.argv[2])
processes = int(str(sys.argv[3]))

print("start initialising")
ds = sfaira.data.dataloaders.DatasetSuperGroupSfaira(
    path=path, meta_path=path_meta, cache_path=path_meta
)
print("start flattening")
dsg = ds.flatten()  # need to flatten in this case to parallelise across Groups and not just within.
print("start loading")
dsg.load(
    celltype_version=None,
    annotated_only=False,
    match_to_reference=None,
    remove_gene_version=True,
    load_raw=True,
    allow_caching=False,
    processes=processes,
    func=write_meta,
    kwargs_func={"args1": path_meta},
)
