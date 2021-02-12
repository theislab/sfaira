import os
import sfaira
import sys
import tensorflow as tf

print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

path = str(sys.argv[1])
fn = str(sys.argv[2])
genome = str(sys.argv[3])

path_meta = os.path.join(path, "meta")
ds = sfaira.data.dataloaders.DatasetSuperGroupSfaira(
    data_path=path, meta_path=path_meta, cache_path=path_meta
)
ds.subset(key="organism", values=["mouse"])
ds.load_all_tobacked(
    fn_backed=fn,
    genome=genome,
    shuffled=False,
    as_dense=False,
    annotated_only=False
)
