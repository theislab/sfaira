import sfaira
import sys
import tensorflow as tf


def write_meta(x, path_meta):
    print("entering function")
    x.write_meta(dir_out=path_meta)
    return None

print(tf.__version__)

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
    kwargs_func={"path_meta": path_meta},
)
