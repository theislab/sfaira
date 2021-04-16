import os
import sfaira
import sys

# Set global variables.
print("sys.argv", sys.argv)

data_path = str(sys.argv[1])
path_meta = str(sys.argv[2])
path_cache = str(sys.argv[3])
path_out = str(sys.argv[4])
schema = str(sys.argv[5])
dois = str(sys.argv[6])

path_cache = path_cache if path_cache != "None" else None

for x in dois.split(","):
    ds = sfaira.data.dataloaders.Universe(
        data_path=data_path, meta_path=path_meta, cache_path=path_cache
    )
    ds.subset(key="doi", values=[x])
    ds.load(
        match_to_reference=None,
        remove_gene_version=True,
        load_raw=False,
        allow_caching=True,
        set_metadata=False,
    )
    if schema == "cellxgene":
        ds.streamline_features(remove_gene_version=True, match_to_reference=True, subset_genes_to_type=None)
    ds.streamline_metadata(format=schema.lower(), allow_uns_sfaira=True, clean_obs=False, clean_var=True, clean_uns=False)
    assert len(ds.dataset_groups) == 1, len(ds.dataset_groups)
    dsg = ds.dataset_groups[0]
    for k, v in dsg.datasets.items():
        fn = v.doi_cleaned_id + ".h5ad"
        v.adata.write_h5ad(os.path.join(path_out, fn))
