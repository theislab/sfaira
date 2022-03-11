import anndata
import numpy as np
import pandas as pd
import scanpy as sc
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
    if organism == "musmusculus":
        adata.uns["organism"] = "Mus musculus"
        adata.uns["organism_ontology_term_id"] = "NCBITaxon:10090"
    elif organism == "homosapiens":
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


def cellxgene_export_adaptor_2_0_0(adata: anndata.AnnData, adata_ids: AdataIdsCellxgene, layer_key_counts: str,
                                   layer_key_proc: str, obs_keys_batch, mask_portal_fields: bool = True,
                                   title: Union[str, None] = None, **kwargs) -> anndata.AnnData:
    """
    Cellxgene-schema 2.0.0.

    Documented here: https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/2.0.3/schema.md
    """
    obs_keys_autofill = [getattr(adata_ids, x) for x in adata_ids.ontology_constrained]
    uns_keys_autofill = []
    var_keys_autofill = [adata_ids.feature_symbol]
    raw_var_keys_autofill = [adata_ids.feature_symbol]
    # 1) Modify .uns
    if layer_key_proc == "X":
        add_proc = False
        adata.uns["X_normalization"] = "custom"
    elif layer_key_counts == "X":
        assert layer_key_proc is None
        # This implies that a auto log-normalised layer will be added below.
        add_proc = True
        adata.uns["X_normalization"] = "log_normalized"
        adata.uns["X_approximate_distribution"] = "normal"
    else:
        assert False
    adata.uns["schema_version"] = "2.0.0"
    adata.uns[adata_ids.author] = {
        "name": "sfaira",
        "email": "https://github.com/theislab/sfaira/issues",
        "institution": "sfaira",
    }
    if obs_keys_batch is not None:
        adata.uns["batch_condition"] = obs_keys_batch.split("*")
        if adata_ids.tech_sample in adata.obs.columns:
            # Delete curated columns as cellxgene does not require this in obs, it only requires the reference to the
            # original column in .uns.
            del adata.obs[adata_ids.tech_sample]
    if title is not None:
        if not isinstance(title, str):
            raise ValueError(f"found type {type(title)} for title, require string or None")
        adata.uns[adata_ids.title] = title
    # Remove overloading:
    adata.uns = dict(adata.uns)
    # 2) Modify .obs
    # a) Correct unknown cell type entries:
    adata.obs[adata_ids.cell_type] = [
        x if x not in [adata_ids.unknown_metadata_identifier, adata_ids.not_a_cell_celltype_identifier]
        else "native cell" for x in adata.obs[adata_ids.cell_type]]
    adata.obs[adata_ids.cell_type + adata_ids.onto_id_suffix] = [
        x if x not in [adata_ids.unknown_metadata_identifier, adata_ids.not_a_cell_celltype_identifier]
        else "CL:0000003" for x in adata.obs[adata_ids.cell_type + adata_ids.onto_id_suffix]]
    # b) Correct unknown disease entries:
    adata.obs[adata_ids.disease] = [
        x if x not in [adata_ids.unknown_metadata_identifier] else "healthy"
        for x in adata.obs[adata_ids.disease]]
    adata.obs[adata_ids.disease + adata_ids.onto_id_suffix] = [
        x if x not in [adata_ids.unknown_metadata_identifier] else "PATO:0000461"
        for x in adata.obs[adata_ids.disease + adata_ids.onto_id_suffix]]
    # Reorder data frame to put ontology columns first:
    cellxgene_cols = [getattr(adata_ids, x) for x in adata_ids.ontology_constrained] + \
                     [getattr(adata_ids, x) for x in adata_ids.obs_keys
                      if x not in adata_ids.ontology_constrained] + \
                     [getattr(adata_ids, x) + adata_ids.onto_id_suffix
                      for x in adata_ids.ontology_constrained]
    adata.obs = adata.obs[cellxgene_cols + [x for x in adata.obs.columns if x not in cellxgene_cols]]
    # 3) Modify .X.
    # Check if .X is counts: The conversion are based on the assumption that .X is csr.
    assert isinstance(adata.X, scipy.sparse.csr_matrix), type(adata.X)
    # Add standard processed count layer (log-normalized) if only counts are supplied:
    # Note that processed counts are in .X if they are supplied along-side raw counts at this point already (this
    # happens in feature streamlining). Raw counts are only not in .X if they processed counts are not supplied.
    # In this exception case, we generate log-normalized data suited for plotting here.
    if add_proc:
        adata.raw = adata
        # Log-normalise values in .X
        sc.pp.normalize_per_cell(adata)
        sc.pp.log1p(adata)
        # This key may be written by log1p above:
        if "log1p" in adata.uns.keys():
            del adata.uns["log1p"]
    # TODO skipped for now:
    # Round numbers in count object:
    # if layer_key_proc is not None:
    #    count_values = np.unique(np.asarray(adata.layers[adata_ids.layer_counts].todense()))
    #    if not np.all(count_values % 1. == 0.):
    #        print(f"WARNING: not all count entries were counts, "
    #              f"the maximum deviation from integer is "
    #              f"{np.max([x % 1. if x % 1. < 0.5 else 1. - x % 1. for x in count_values])}. "
    #              f"The count matrix is rounded.")
    #        adata.layers[adata_ids.layer_counts].data = np.rint(adata.layers[adata_ids.layer_counts].data)
    # 4) Modify .var:
    adata.var[adata_ids.feature_biotype] = pd.Categorical(["gene" for _ in range(adata.n_vars)])
    adata.var.index = adata.var[adata_ids.feature_id].values
    # TODO remove this after clean up such that this does not happend after streamlining:
    if adata_ids.feature_id in adata.var.columns:
        del adata.var[adata_ids.feature_id]
    if adata.raw is not None:
        adata.raw.var[adata_ids.feature_biotype] = pd.Categorical(["gene" for _ in range(adata.raw.n_vars)])
        if adata_ids.feature_id in adata.raw.var.columns:
            del adata.raw.var[adata_ids.feature_id]
    if adata.raw is not None and adata.n_vars < adata.raw.n_vars:
        # extend X by zero columns and add filtered attribute:
        raise NotImplementedError()
    else:
        adata.var[adata_ids.feature_is_filtered] = False
        # Assert that genes are ordered the same way in raw and processed:
        assert np.all(adata.var.index == adata.raw.var.index)
    # Modify ensembl ID writing:
    # adata.var[adata_ids.feature_id] = ["G:".join(x.split("G")) for x in adata.var[adata_ids.feature_id]]
    # adata.var.index = ["G:".join(x.split("G")) for x in adata.var.index]
    # 5) Take out elements that are auto-filled by cellxgene upload interface:
    if mask_portal_fields:
        for k in obs_keys_autofill:
            del adata.obs[k]
        for k in uns_keys_autofill:
            del adata.uns[k]
        for k in var_keys_autofill:
            del adata.var[k]
        for k in raw_var_keys_autofill:
            if adata.raw is not None and k in adata.raw.var.columns:
                del adata.raw.var[k]
    # 6) Check if default embedding is present, add in otherwise:
    # First check if any pre-computed embedding is given.
    # If that is not the case, compute a default UMAP.
    # Define hierarchy of embeddings accepted as defaults, first one matched will be chosen:
    default_embedding_names = ["X_umap", "X_tsne", "X_draw_graph_fa"]
    if adata.uns[adata_ids.default_embedding] is None or \
            adata.uns[adata_ids.default_embedding] == adata_ids.unknown_metadata_identifier:
        found_default = False
        counter = 0
        while not found_default and counter < len(default_embedding_names):
            if default_embedding_names[counter] in adata.obsm.keys():
                adata.uns[adata_ids.default_embedding] = default_embedding_names[counter]
                found_default = True
            counter += 1
        if not found_default and adata.n_obs > 10:
            adata_embedding = adata.copy()
            sc.pp.pca(adata_embedding)
            sc.pp.neighbors(adata_embedding)
            sc.tl.umap(adata_embedding)
            adata.obsm["X_umap"] = adata_embedding.obsm["X_umap"]
            adata.uns[adata_ids.default_embedding] = "X_umap"
    return adata
