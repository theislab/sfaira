import os
import sys

# Any data loader here to extract path:
from sfaira.consts import OC
from sfaira.data import load_store
from sfaira.versions.metadata import CelltypeUniverse, OntologyCl
from sfaira.data.dataloaders.obs_utils import is_custom


# Set global variables.
print("sys.argv", sys.argv)

store_path = str(sys.argv[1])
config_path = str(sys.argv[2])
out_path = str(sys.argv[3])

col_name_annot = "cell_ontology_class"
cl_branch = "v2021-02-01"

for f in os.listdir(config_path):
    fn = os.path.join(config_path, f)
    if os.path.isfile(fn):  # only files
        # Narrow down to supported file types:
        if f.split(".")[-1] == "pickle" and f.startswith("config_"):
            print(f"Writing target universe for {f}")
            organism = f.split("_")[1]
            organ = f.split("_")[2].split(".")[0]
            store = load_store(cache_path=store_path)
            store.load_config(fn=fn)
            celltypes_found = set([])
            for k, idx in store.indices.items():
                if col_name_annot not in store.adatas[k].obs.columns:
                    print(f"WARNING: annotation column {col_name_annot} not found in  {k}, skipping.")
                else:
                    if len(idx) > 0:
                        celltypes_found = celltypes_found.union(
                            set(store.adatas[k].obs[col_name_annot].values[idx].tolist())
                        )
            celltypes_found = sorted(list(celltypes_found - {store._adata_ids_sfaira.unknown_metadata_identifier,
                                                             store._adata_ids_sfaira.not_a_cell_celltype_identifier}))
            celltypes_found = [
                x for x in celltypes_found
                if not is_custom(x, store._adata_ids_sfaira)
            ]
            if len(celltypes_found) == 0:
                print(f"WARNING: No cells found for {organism} {organ}, skipping.")
            else:
                celltypes_universe = CelltypeUniverse(
                    cl=OntologyCl(branch=cl_branch),
                    uberon=OC.organ,
                )
                celltypes_found = celltypes_universe.onto_cl.get_effective_leaves(x=celltypes_found)
                celltypes_universe.write_target_universe(
                    fn=os.path.join(out_path, f"targets_{organism}_{organ}.csv"),
                    x=celltypes_found,
                )
