import anndata
import numpy as np
import scipy.sparse
from typing import Union

from sfaira.consts.adata_fields import AdataIdsCellxgene

DEFAULT_CELLXGENE_VERSION = "1_1_0"


def cellxgene_export_adaptor(adata: anndata.AnnData, adata_ids: AdataIdsCellxgene, version: Union[None, str]) \
        -> anndata.AnnData:
    """
    Projects a streamlined data set to the export-ready cellxgene format.
    """
    if version is None:
        version = DEFAULT_CELLXGENE_VERSION
    if version == "1_1_0":
        return cellxgene_export_adaptor_1_1_0(adata=adata, adata_ids=adata_ids)
    else:
        raise ValueError(f"Did not recognise cellxgene schema version {version}")


def cellxgene_export_adaptor_1_1_0(adata: anndata.AnnData, adata_ids: AdataIdsCellxgene) -> anndata.AnnData:
    """
    Cellxgene-schema 1.1.0
    """
    # Check input object characteristics:
    cellxgene_cols = [getattr(adata_ids, x) for x in adata_ids.ontology_constrained] + \
                     [getattr(adata_ids, x) for x in adata_ids.obs_keys
                      if x not in adata_ids.ontology_constrained] + \
                     [getattr(adata_ids, x) + adata_ids.onto_id_suffix
                      for x in adata_ids.ontology_constrained]
    for x in cellxgene_cols:
        if x not in adata.obs.columns:
            raise ValueError(f"Cannot streamlined data set {adata.uns['id']} to cellxgene format because meta data {x} "
                             f"is missing and the corresponding .obs column could not be written.\n"
                             f"Columns found were {adata.obs.columns}.")
    # 1) Modify .uns
    adata.uns["layer_descriptions"] = {"X": "raw"}
    adata.uns["version"] = {
        "corpora_encoding_version": "0.1.0",
        "corpora_schema_version": "1.1.0",
    }
    adata.uns["contributors"] = {
        "name": "sfaira",
        "email": "https://github.com/theislab/sfaira/issues",
        "institution": "sfaira",
    }
    # TODO port this into organism ontology handling.
    # Infer organism from adata object.
    organism = np.unique(adata.obs[adata_ids.organism].values)[0]
    if organism == "mouse":
        adata.uns["organism"] = "Mus musculus"
        adata.uns["organism_ontology_term_id"] = "NCBITaxon:10090"
    elif organism == "human":
        adata.uns["organism"] = "Homo sapiens"
        adata.uns["organism_ontology_term_id"] = "NCBITaxon:9606"
    else:
        raise ValueError(f"organism {organism} currently not supported by cellxgene schema")
    # 2) Modify .obs
    # Correct unknown cell type entries:
    adata.obs[adata_ids.cell_type] = [
        x if x not in [adata_ids.unknown_metadata_identifier, adata_ids.not_a_cell_celltype_identifier]
        else "native cell" for x in adata.obs[adata_ids.cell_type]]
    adata.obs[adata_ids.cell_type + adata_ids.onto_id_suffix] = [
        x if x not in [adata_ids.unknown_metadata_identifier, adata_ids.not_a_cell_celltype_identifier]
        else "CL:0000003" for x in adata.obs[adata_ids.cell_type + adata_ids.onto_id_suffix]]
    # Reorder data frame to put ontology columns first:
    adata.obs = adata.obs[cellxgene_cols + [x for x in adata.obs.columns if x not in cellxgene_cols]]
    # 3) Modify .X
    # Check if .X is counts: The conversion are based on the assumption that .X is csr.
    assert isinstance(adata.X, scipy.sparse.csr_matrix), type(adata.X)
    count_values = np.unique(np.asarray(adata.X.todense()))
    if not np.all(count_values % 1. == 0.):
        print(f"WARNING: not all count entries were counts, "
              f"the maximum deviation from integer is "
              f"{np.max([x % 1. if x % 1. < 0.5 else 1. - x % 1. for x in count_values])}. "
              f"The count matrix is rounded.")
        adata.X.data = np.rint(adata.X.data)
    return adata


def cellxgene_export_adaptor_2_0_0(adata: anndata.AnnData, adata_ids: AdataIdsCellxgene) -> anndata.AnnData:
    """
    Cellxgene-schema 2.0.0
    """
    adata = cellxgene_export_adaptor_1_1_0(adata=adata, adata_ids=adata_ids)
    adata.var["feature_biotype"] = "gene"
