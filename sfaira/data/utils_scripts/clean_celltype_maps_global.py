import sys
import tensorflow as tf

# Any data loader here to extract path:
from sfaira.data.dataloaders.loaders import DatasetSuperGroupLoaders

print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

data_path = str(sys.argv[1])
path_meta = str(sys.argv[2])
path_cache = str(sys.argv[3])
processes = int(str(sys.argv[4]))

dsgl = DatasetSuperGroupLoaders(
    data_path=data_path,
    meta_path=path_meta,
    cache_path=path_cache
)

for x in dsgl.dataset_groups:
    print(x.ids)
    x.clean_ontology_class_map()
