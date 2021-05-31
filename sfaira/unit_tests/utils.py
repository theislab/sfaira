import anndata
import numpy as np
import os
from typing import Tuple, Union

from sfaira.consts import AdataIdsSfaira, OCS
from sfaira.data import Universe
from sfaira.versions.metadata import OntologyOboCustom


def simulate_anndata(genes, n_obs, targets=None, assays=None, obo: Union[None, OntologyOboCustom] = None) -> \
        anndata.AnnData:
    """
    Simulate basic data example.

    :return: AnnData instance.
    """
    adata_ids_sfaira = AdataIdsSfaira()
    data = anndata.AnnData(
        np.random.randint(low=0, high=100, size=(n_obs, len(genes))).astype(np.float32)
    )
    if assays is not None:
        data.obs[adata_ids_sfaira.assay_sc] = [
            assays[np.random.randint(0, len(assays))]
            for _ in range(n_obs)
        ]
    if targets is not None:
        data.obs[adata_ids_sfaira.cellontology_class] = [
            targets[np.random.randint(0, len(targets))]
            for _ in range(n_obs)
        ]
        if obo is None:
            data.obs[adata_ids_sfaira.cellontology_id] = [
                OCS.cellontology_class.convert_to_id(x)
                if x not in [adata_ids_sfaira.unknown_celltype_identifier,
                             adata_ids_sfaira.not_a_cell_celltype_identifier]
                else x
                for x in data.obs[adata_ids_sfaira.cellontology_class].values
            ]
        else:
            data.obs[adata_ids_sfaira.cellontology_id] = [
                obo.convert_to_id(x)
                if x not in [adata_ids_sfaira.unknown_celltype_identifier,
                             adata_ids_sfaira.not_a_cell_celltype_identifier]
                else x
                for x in data.obs[adata_ids_sfaira.cellontology_class].values
            ]
    data.var[adata_ids_sfaira.gene_id_ensembl] = genes
    return data


def cached_store_writing(dir_data, dir_meta, assembly, organism: str = "mouse", organ: str = "lung",
                         store_format: str = "h5ad", return_ds: bool = False) -> Union[str, Tuple[str, Universe]]:
    """
    Writes a store if it does not already exist.

    :return: Path to store.
    """
    adata_ids_sfaira = AdataIdsSfaira()
    store_path = os.path.join(dir_data, "store")
    if not os.path.exists(store_path):
        os.mkdir(store_path)
    ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
    ds.subset(key=adata_ids_sfaira.organism, values=[organism])
    ds.subset(key=adata_ids_sfaira.organ, values=[organ])
    # Only load files that are not already in cache.
    anticipated_files = np.unique([
        v.doi[0] if isinstance(v.doi, list) else v.doi for k, v in ds.datasets.items()
        if (not os.path.exists(os.path.join(store_path, v.doi_cleaned_id + "." + store_format)) and
            store_format == "h5ad") or
           (not os.path.exists(os.path.join(store_path, v.doi_cleaned_id)) and store_format == "dao")
    ]).tolist()
    ds.subset(key=adata_ids_sfaira.doi, values=anticipated_files)
    ds.load(allow_caching=True)
    ds.streamline_features(remove_gene_version=True, match_to_reference={organism: assembly},
                           subset_genes_to_type="protein_coding")
    ds.streamline_metadata(schema="sfaira", uns_to_obs=True, clean_obs=True, clean_var=True, clean_uns=True,
                           clean_obs_names=True)
    if store_format == "zarr":
        compression_kwargs = {"compressor": "default", "overwrite": True, "order": "C"}
    else:
        compression_kwargs = {}
    ds.write_distributed_store(dir_cache=store_path, store_format=store_format, dense=store_format == "dao",
                               chunks=128, compression_kwargs=compression_kwargs)
    if return_ds:
        ds = Universe(data_path=dir_data, meta_path=dir_meta, cache_path=dir_data)
        ds.subset(key=adata_ids_sfaira.organism, values=[organism])
        ds.subset(key=adata_ids_sfaira.organ, values=[organ])
        ds.load(allow_caching=True)
        ds.streamline_features(remove_gene_version=True, match_to_reference={organism: assembly},
                               subset_genes_to_type="protein_coding")
        ds.streamline_metadata(schema="sfaira", uns_to_obs=True, clean_obs=True, clean_var=True, clean_uns=True,
                               clean_obs_names=True)
        return store_path, ds
    else:
        return store_path
