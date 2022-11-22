import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from typing import Union

from sfaira.consts.adata_fields import AdataIdsCellxgene, AdataIdsCellxgene_v3_0_0
from sfaira.consts.schema import DEFAULT_SCHEMA
from sfaira.data.dataloaders.ontology_access import get_ontology
from sfaira.versions.metadata.base import OntologyEfo


def cellxgene_export_adaptor(adata: anndata.AnnData, adata_ids: AdataIdsCellxgene, version: Union[None, str],
                             **kwargs) -> anndata.AnnData:
    """
    Projects a streamlined data set to the export-ready cellxgene format.
    """
    if version is None:
        version = DEFAULT_SCHEMA
    if version == "2_0_0":
        return cellxgene_export_adaptor_2_0_0(adata=adata, adata_ids=adata_ids, **kwargs)
    elif version == "3_0_0":
        return cellxgene_export_adaptor_3_0_0(adata=adata, adata_ids=adata_ids, **kwargs)
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
    if adata.raw is not None:
        adata.raw.var[adata_ids.feature_biotype] = pd.Categorical(["gene" for _ in range(adata.raw.n_vars)])
    if adata.raw is not None and adata.n_vars < adata.raw.n_vars:
        # extend X by zero columns and add filtered attribute:
        raise NotImplementedError()
    elif adata.raw is not None:
        adata.var[adata_ids.feature_is_filtered] = False
        # Assert that genes are ordered the same way in raw and processed:
        assert np.all(adata.var.index == adata.raw.var.index)
    # Modify ensembl ID writing:
    # adata.var[adata_ids.feature_id] = ["G:".join(x.split("G")) for x in adata.var[adata_ids.feature_id]]
    # adata.var.index = ["G:".join(x.split("G")) for x in adata.var.index]
    # 5) Take out elements that are autofilled by cellxgene upload interface:
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


def match_supsension_and_efo(adata, adata_ids, efo_ontology: OntologyEfo, valid_combinations: dict):
    # Create matching table that takes hierarchical relation in EFO into account: valid_combinations contains
    # high level nodes and their matching to suspension type.
    efo_suspension_table = adata.obs[
        [adata_ids.assay_sc, adata_ids.assay_sc + adata_ids.onto_id_suffix, adata_ids.suspension_type]
    ].drop_duplicates()
    efo_map = {}
    for efo in efo_suspension_table[adata_ids.assay_sc + adata_ids.onto_id_suffix].values:
        efo_map[efo] = "na"
        for k, v in valid_combinations.items():
            if efo_ontology.is_a(is_=efo, a_=k):
                if efo in efo_map.keys():
                    raise Warning(f"found multiple suspension matches for EFO {efo}.")
                efo_map[efo] = v[0]
    # Picks first match, ie prioritises cell if both cell and nucleus are possible.
    adata.obs[adata_ids.suspension_type] = [
        efo_map[x] for x in adata.obs[adata_ids.assay_sc + adata_ids.onto_id_suffix].values
    ]
    return adata


VALID_EFO_SUS = {
    "3_0_0": {
        "EFO:0030080": ["cell", "nucleus"],
        "EFO:0007045": ["nucleus"],
        "EFO:0010010": ["cell", "nucleus"],
        "EFO:0009294": ["cell"],
        "EFO:0008720": ["nucleus"],
        "EFO:0008722": ["cell", "nucleus"],
        "EFO:0030002": ["cell"],
        "EFO:0008853": ["cell"],
        "EFO:0030026": ["nucleus"],
        "EFO:0010550": ["cell", "nucleus"],
        "EFO:0008919": ["cell"],
        "EFO:0010184": ["cell", "nucleus"],
        "EFO:0009918": ["na"],
        "EFO:0008939": ["nucleus"],
        "EFO:0030027": ["nucleus"],
        "EFO:0700000": ["na"],
        "EFO:0030005": ["na"],
    }
}


def cellxgene_export_adaptor_3_0_0(adata: anndata.AnnData, adata_ids: AdataIdsCellxgene_v3_0_0,
                                   layer_key_counts: str, layer_key_proc: str,
                                   obs_keys_batch, mask_portal_fields: bool = True, title: Union[str, None] = None,
                                   **kwargs) -> anndata.AnnData:
    """
    Cellxgene-schema 3.0.0.

    Documented here: https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md
    """
    ontology_constrained_with_id = [x for x in adata_ids.ontology_constrained if x not in ['suspension_type']]
    # These fields are filled by cellxgene based on other fields:
    obs_keys_autofill = [getattr(adata_ids, x) for x in ontology_constrained_with_id]
    uns_keys_autofill = []
    var_keys_autofill = [adata_ids.feature_symbol]
    raw_var_keys_autofill = [adata_ids.feature_symbol]

    # 1) Modify .uns
    if layer_key_proc == "X":
        add_proc = False
    elif layer_key_counts == "X":
        assert layer_key_proc is None
        # This implies that a auto log-normalised layer will be added below.
        add_proc = True
        adata.uns["X_approximate_distribution"] = "normal"
    else:
        assert False
    adata.uns["schema_version"] = "3.0.0"
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
    # Make donor categorical.
    adata.obs[adata_ids.individual] = pd.Categorical(adata.obs[adata_ids.individual].values.tolist())
    # Suspension is a custom ontology and does not require an ID. The column itself must be categorical.
    # Add this annotation here if it was not set before.
    if np.all(adata.obs[adata_ids.suspension_type].values == adata_ids.unknown_metadata_identifier):
        adata = match_supsension_and_efo(adata=adata, adata_ids=adata_ids, efo_ontology=get_ontology(k="assay_sc"),
                                         valid_combinations=VALID_EFO_SUS["3_0_0"])
    adata.obs[adata_ids.suspension_type] = pd.Categorical(adata.obs[adata_ids.suspension_type].values.tolist())
    del adata.obs[adata_ids.suspension_type + adata_ids.onto_id_suffix]
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
                      for x in ontology_constrained_with_id]
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
    # 4) Modify .var:
    if adata.raw is not None and adata.n_vars < adata.raw.n_vars:
        # extend X by zero columns and add filtered attribute:
        raise NotImplementedError()
    elif adata.raw is not None:
        adata.var[adata_ids.feature_is_filtered] = False
        # Assert that genes are ordered the same way in raw and processed:
        assert np.all(adata.var.index == adata.raw.var.index)
    # Modify ensembl ID writing:
    # adata.var[adata_ids.feature_id] = ["G:".join(x.split("G")) for x in adata.var[adata_ids.feature_id]]
    # adata.var.index = ["G:".join(x.split("G")) for x in adata.var.index]
    # 5) Take out elements that are autofilled by cellxgene upload interface:
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
