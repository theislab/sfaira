import sfaira
import sys
import tensorflow as tf

print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

path = str(sys.argv[1])
path_meta = str(sys.argv[2])

ds = sfaira.data.dataloaders.DatasetSuperGroupSfaira(
    path=path, meta_path=path_meta, cache_path=path_meta
)
for dsg in list(ds.dataset_groups):
    for k in dsg.ids:
        print(k)
        dsg.datasets[k].write_meta(dir_out=path_meta)
