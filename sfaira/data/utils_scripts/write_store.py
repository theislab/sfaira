import argparse
import os

import sfaira


parser = argparse.ArgumentParser()
parser.add_argument('--data_path', type=str)
parser.add_argument('--path_meta', type=str)
parser.add_argument('--path_cache', type=str)
parser.add_argument('--path_store', type=str)
parser.add_argument('--store_type', type=str)
parser.add_argument('--chunks', type=int, default=128)
parser.add_argument('--shuffle_data', type=lambda x: str(x).lower() in ['true', '1', 'yes'], default=True)
args = parser.parse_args()

# On disk format hyperparameters
if args.store_type == "h5ad":
    # Write sparse arrays in h5ad.
    kwargs = {"dense": False}
    compression_kwargs = {}
elif args.store_type == "dao":
    # Write dense arrays in zarr.
    kwargs = {"dense": True, "chunks": args.chunks, "shuffle_data": args.shuffle_data}
    compression_kwargs = {"compressor": "default", "overwrite": True, "order": "C"}
else:
    assert False, args.store_type

universe = sfaira.data.dataloaders.Universe(data_path=args.data_path,
                                            meta_path=args.path_meta,
                                            cache_path=args.path_cache)

for k, ds in universe.datasets.items():
    if args.store_type == "h5ad":
        fn_store = os.path.join(args.path_store, ds.doi_cleaned_id + ".h5ad")
    elif args.store_type == "dao":
        fn_store = os.path.join(args.path_store, ds.doi_cleaned_id)
    else:
        assert False
    if os.path.exists(fn_store):
        print(f"SCRIPT skipping {k}")
    else:
        print(f"SCRIPT loading {k}")
        ds.load(
            load_raw=False,
            allow_caching=True,
        )
        ds.streamline_features(
            remove_gene_version=True,
            match_to_release={"Homo sapiens": "104", "Mus musculus": "104"},
            subset_genes_to_type="protein_coding"
        )
        ds.streamline_metadata(schema="sfaira", clean_obs=True, clean_var=True, clean_uns=True, clean_obs_names=True)
        ds.write_distributed_store(dir_cache=args.path_store, store_format=args.store_type,
                                   compression_kwargs=compression_kwargs, **kwargs)
        ds.clear()
