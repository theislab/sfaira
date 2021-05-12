import anndata
import numpy as np
import os

from sfaira.data import DistributedStore, Universe


def simulate_anndata(genes, n_obs, targets=None, assays=None) -> anndata.AnnData:
    """
    Simulate basic data example.

    :return: AnnData instance.
    """
    data = anndata.AnnData(
        np.random.randint(low=0, high=100, size=(n_obs, len(genes))).astype(np.float32)
    )
    if assays is not None:
        data.obs["assay_sc"] = [
            assays[np.random.randint(0, len(assays))]
            for i in range(n_obs)
        ]
    if targets is not None:
        data.obs["cell_ontology_class"] = [
            targets[np.random.randint(0, len(targets))]
            for i in range(n_obs)
        ]
    data.var["ensembl"] = genes
    return data


def cached_store_writing(dir_data, dir_meta, assembly, organism: str = "mouse", organ: str = "lung",
                         store_format: str = "h5ad") -> str:
    """
    Writes a store if it does not already exist.

    :return: Path to store.
    """
    store_path = os.path.join(dir_data, "store")
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key="organism", values=[organism])
    ds.subset(key="organ", values=[organ])
    # Only load files that are not already in cache.
    anticipated_files = np.unique([
        v.doi[0] if isinstance(v.doi, list) else v.doi for k, v in ds.datasets.items()
        if (not os.path.exists(os.path.join(store_path, v.doi_cleaned_id + "." + store_format)) and store_format == "h5ad") or
           (not os.path.exists(os.path.join(store_path, v.doi_cleaned_id)) and store_format == "zarr")
    ]).tolist()
    ds.subset(key="doi", values=anticipated_files)
    ds.load(allow_caching=True)
    ds.streamline_features(remove_gene_version=True, match_to_reference={organism: assembly},
                           subset_genes_to_type="protein_coding")
    ds.streamline_metadata(schema="sfaira", uns_to_obs=True, clean_obs=True, clean_var=True, clean_uns=True,
                           clean_obs_names=True)
    ds.write_distributed_store(dir_cache=store_path, store_format=store_format, dense=store_format == "zarr")
    return store_path
