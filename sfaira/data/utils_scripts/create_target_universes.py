import numpy as np
import os
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
        if f.split(".")[-1] == "pickle" and f.startswith("config_"):
            print(f"Writing {f}")
            organism = f.split("_")[1]
            organ = f.split("_")[2]
            store = DistributedStore(cache_path=store_path)
            store.load_config(fn=fn)
            celltypes_found = set([])
            for adata in store.adatas:
                celltypes_found = celltypes_found.union(set(adata.obs["cell_ontology_class"].values.tolist()))
            celltypes_found = np.sort(list(celltypes_found - set([
                store._adata_ids_sfaira.unknown_celltype_identifier,
                store._adata_ids_sfaira.not_a_cell_celltype_identifier
            ]))).tolist()
            if len(celltypes_found) == 0:
                print(f"WARNING: No cells found for {organism} {organ}, skipping.")
            else:
                celltypes_found = store.celltypes_universe.onto_cl.get_effective_leaves(x=celltypes_found)
                store.celltypes_universe.write_target_universe(
                    fn=os.path.join(config_path, f"targets_{organism}_{organ}.csv"),
                    x=celltypes_found,
                )
