import anndata
import numpy as np
import pandas as pd
import scipy.sparse
from typing import Union

from sfaira.consts.adata_fields import AdataIdsCellxgene

DEFAULT_CELLXGENE_VERSION = "2_0_0"


def cellxgene_export_adaptor(adata: anndata.AnnData, adata_ids: AdataIdsCellxgene, version: Union[None, str],
                             **kwargs) -> anndata.AnnData:
    """
    Projects a streamlined data set to the export-ready cellxgene format.
    """
    if version is None:
        version = DEFAULT_CELLXGENE_VERSION
    if version == "2_0_0":
        return cellxgene_export_adaptor_2_0_0(adata=adata, adata_ids=adata_ids, **kwargs)
    else:
        raise ValueError(f"Did not recognise cellxgene schema version {version}")


def cellxgene_export_adaptor_1_1_0(adata: anndata.AnnData, adata_ids: AdataIdsCellxgene, **kwargs) -> anndata.AnnData:
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


def cellxgene_export_adaptor_2_0_0(adata: anndata.AnnData, adata_ids: AdataIdsCellxgene, obs_keys_batch,
                                   mask_portal_fields: bool = True, **kwargs) -> anndata.AnnData:
    """
    Cellxgene-schema 2.0.0.

    Documented here: https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/2.0.0/schema.md
    """
    obs_keys_autofill = [getattr(adata_ids, x) for x in adata_ids.ontology_constrained]
    uns_keys_autofill = []
    var_keys_autofill = [adata_ids.feature_symbol, adata_ids.feature_reference]
    raw_var_keys_remove = [adata_ids.feature_is_filtered]
    x_is_raw = True
    # 1) Modify .uns
    if x_is_raw:
        adata.uns["layer_descriptions"] = {"X": "raw"}
        adata.uns["X_normalization"] = "none"
        adata.uns["X_approximate_distribution"] = "count"
    else:
        raise NotImplementedError()
    adata.uns["schema_version"] = "2.0.0"
    adata.uns[adata_ids.author] = {
        "name": "sfaira",
        "email": "https://github.com/theislab/sfaira/issues",
        "institution": "sfaira",
    }
    if obs_keys_batch is not None:
        adata.uns["batch_condition"] = obs_keys_batch.split("*")
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
    cellxgene_cols = [getattr(adata_ids, x) for x in adata_ids.ontology_constrained] + \
                     [getattr(adata_ids, x) for x in adata_ids.obs_keys
                      if x not in adata_ids.ontology_constrained] + \
                     [getattr(adata_ids, x) + adata_ids.onto_id_suffix
                      for x in adata_ids.ontology_constrained]
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
    # 4) Modify .var:
    adata.var[adata_ids.feature_biotype] = "gene"
    adata.var[adata_ids.feature_reference] = adata.uns[adata_ids.organism]
    adata.var[adata_ids.feature_is_filtered] = False
    # Modify ensembl ID writing:
    #adata.var[adata_ids.feature_id] = ["G:".join(x.split("G")) for x in adata.var[adata_ids.feature_id]]
    #adata.var.index = ["G:".join(x.split("G")) for x in adata.var.index]
    # 5) Take out elements that are auto-filled by cellxgene upload interface:
    if mask_portal_fields:
        for k in obs_keys_autofill:
            del adata.obs[k]
        for k in uns_keys_autofill:
            del adata.uns[k]
        for k in var_keys_autofill:
            del adata.var[k]
    # 6) Modify .raw:
    # Note that this depends on var cleaning in step 5.
    if not x_is_raw:
        var_raw = pd.DataFrame(dict([(k, v) for k, v in adata.var.items() if k not in raw_var_keys_remove]))
        adata.raw = anndata.AnnData(adata.X, var=var_raw)
    print(adata.var)
    return adata
