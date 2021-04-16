import os
import sfaira
import sys
import tensorflow as tf

# Any data loader here to extract path:
from sfaira.data import DistributedStore

print(tf.__version__)


# Set global variables.
print("sys.argv", sys.argv)

store_path = str(sys.argv[1])
config_path = str(sys.argv[2])
out_path = str(sys.argv[3])


for f in os.listdir(config_path):
    fn = os.path.join(config_path, f)
    if os.path.isfile(fn):  # only files
        # Narrow down to supported file types:
        if f.split(".")[-1] == "csv" and f.startswith("config_"):
            print(f"Writing {f}")
            organism = f.split("_")[1]
            organ = f.split("_")[2]
            store = DistributedStore(cache_path=store_path)
            store.load_config(fn=fn)
            store.write_config(os.path.join(config_path, f"targets_{organism}_{organ}.csv"))
