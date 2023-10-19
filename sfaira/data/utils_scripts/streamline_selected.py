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

path_cache = path_cache if path_cache.lower() != "none" else None
path_meta = path_meta if path_meta.lower() != "none" else None

for doi in dois.split(","):
    ds = sfaira.data.dataloaders.Universe(
        data_path=data_path, meta_path=path_meta, cache_path=path_cache
    )
    ds.subset(key="doi", values=[doi])
    ds.load(
        load_raw=False,
        allow_caching=True,
    )
    if schema == "cellxgene":
        ds.streamline_var(match_to_release=None, schema="cellxgene:" + "3_0_0")
        ds.streamline_obs_uns(
            schema=schema.lower(),
            clean_obs=False,
            clean_var=True,
            clean_uns=True,
            clean_obs_names=False,
            keep_orginal_obs=False,
            keep_symbol_obs=True,
            keep_id_obs=True,
        )
        ds.collapse_counts()
    assert len(ds.dataset_groups) == 1, len(ds.dataset_groups)
    dsg = ds.dataset_groups[0]
    for k, v in dsg.datasets.items():
        fn = v.id + ".h5ad"
        dir_name = v.directory_formatted_doi
        if not os.path.exists(os.path.join(path_out, dir_name)):
            os.makedirs(os.path.join(path_out, dir_name))
        v.adata.write_h5ad(os.path.join(path_out, dir_name, fn))
