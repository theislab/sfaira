from typing import Dict, List, Union

import anndata
import numpy as np
import pandas as pd

from sfaira.consts import AdataIds, OCS
from sfaira.data.dataloaders.export_adaptors import cellxgene_export_adaptor
from sfaira.versions.metadata import Ontology, OntologyHierarchical


def streamline_obs_uns(adata: anndata.AnnData,
                       adata_input_ids: AdataIds,
                       adata_target_ids: AdataIds,
                       annotation_container,
                       organism: str,
                       schema: str = "sfaira",
                       clean_obs: Union[bool, list, np.ndarray, tuple] = True,
                       clean_uns: bool = True,
                       clean_obs_names: bool = True,
                       dataset_id: str = "",
                       keep_orginal_obs: Union[bool, list, np.ndarray, tuple] = False,
                       keep_symbol_obs: bool = True,
                       keep_id_obs: bool = True,
                       ontology_class_maps: Dict[str, pd.DataFrame] = {},
                       **kwargs) -> anndata.AnnData:
    """
    Streamline the adata instance to a defined output schema.

    Output format are saved in ADATA_FIELDS* classes.

    Note on ontology-controlled meta data:
    These are defined for a given format in `ADATA_FIELDS*.ontology_constrained`.
    They may appear in three different formats:
        - original (free text) annotation
        - ontology symbol
        - ontology ID
    During streamlining, these ontology-controlled meta data are projected to all of these three different formats.
    The initially annotated column may be any of these and is defined as "{attr}_obs_col".
    The resulting three column per meta data item are named:
        - ontology symbol: "{ADATA_FIELDS*.attr}"
        - ontology ID: {ADATA_FIELDS*.attr}_{ADATA_FIELDS*.onto_id_suffix}"
        - original (free text) annotation: "{ADATA_FIELDS*.attr}_{ADATA_FIELDS*.onto_original_suffix}"

    :param adata:
    :param adata_input_ids:
    :param adata_target_ids:
    :param annotation_container:
    :param organism:
    :param schema: Export format.
        - "sfaira"
        - "cellxgene"
    :param clean_obs: Whether to delete non-streamlined fields in .obs, .obsm and .obsp.
         Alternatively, list of .obs fields to remove.
    :param clean_uns: Whether to delete non-streamlined fields in .uns.
    :param clean_obs_names: Whether to replace obs_names with a string comprised of dataset id and an increasing
        integer.
    :param dataset_id: Identifier of dataset used in logging.
    :param keep_orginal_obs: For ontology-constrained .obs columns, whether to keep a column with original
        annotation. Alternatively, list of original fields to keep, others are removed.
    :param keep_symbol_obs: For ontology-constrained .obs columns, whether to keep a column with ontology symbol
        annotation.
    :param keep_id_obs: For ontology-constrained .obs columns, whether to keep a column with ontology ID annotation.
    :return:
    """
    if isinstance(keep_orginal_obs, tuple):
        keep_orginal_obs = list(keep_orginal_obs)
    elif isinstance(keep_orginal_obs, np.ndarray):
        keep_orginal_obs = keep_orginal_obs.tolist()

    experiment_batch_labels = [getattr(adata_input_ids, x) for x in adata_input_ids.batch_keys]

    # Prepare new .uns dict:
    uns_new = {}
    for k in adata_target_ids.uns_keys:
        if hasattr(annotation_container, k) and getattr(annotation_container, k) is not None:
            val = getattr(annotation_container, k)
        elif hasattr(annotation_container, f"{k}_obs_key") and getattr(annotation_container, f"{k}_obs_key") is not None:
            val = adata.obs[getattr(annotation_container, f"{k}_obs_key")].unique().astype(str).tolist()
        elif getattr(adata_input_ids, k) in adata.obs.columns:
            val = adata.obs[getattr(adata_input_ids, k)].unique().astype(str).tolist()
        else:
            val = None
        while hasattr(val, '__len__') and not isinstance(val, str) and len(val) == 1:  # Unpack nested lists/tuples.
            val = val[0]
        uns_new[getattr(adata_target_ids, k)] = val
    if clean_uns:
        adata.uns = uns_new
    else:
        adata.uns.update(uns_new)

    # Prepare new .obs dataframe
    # Queried meta data may be:
    # 1) in .obs
    #   a) for an ontology-constrained meta data item
    #       I) as free annotation with a term map to an ontology
    #       II) as column with ontology symbols
    #       III) as column with ontology IDs
    #   b) for a non-ontology-constrained meta data item:
    #       I) as free annotation
    # 2) in .uns
    #   b) as elements that are ontology symbols
    #   c) as elements that are ontology IDs
    # .obs annotation takes priority over .uns annotation if both are present.
    # The output columns are:
    # - for an ontology-constrained meta data item "attr":
    #   *  symbols:         "attr"
    #   *  IDs:             "attr" + adata_input_ids.onto_id_suffix
    #   *  original labels: "attr" + adata_input_ids.onto_original_suffix
    # - for a non-ontology-constrained meta data item "attr":
    #   *  original labels: "attr" + adata_input_ids.onto_original_suffix
    obs_new = pd.DataFrame(index=adata.obs.index)
    obs_to_delete = []
    for k in [x for x in adata_target_ids.obs_keys]:
        if k in experiment_batch_labels and getattr(annotation_container, f"{k}_obs_key") is not None:
            # Handle batch-annotation columns which can be provided as a combination of columns separated by an
            # asterisk.
            # The queried meta data are always:
            # 1b-I) a combination of existing columns in .obs
            old_cols = getattr(annotation_container, f"{k}_obs_key")
            batch_cols = []
            for batch_col in old_cols.split("*"):
                if batch_col in adata.obs.columns:
                    batch_cols.append(batch_col)
                else:
                    # This should not occur in single data set loaders (see warning below) but can occur in
                    # streamlined data loaders if not all instances of the streamlined data sets have all columns
                    # in .obs set.
                    print(f"WARNING: attribute {batch_col} of data set {dataset_id} was not found in columns: "
                          f"{adata.obs.columns}.")
            # Build a combination label out of all columns used to describe this group.
            # Add data set label into this label so that these groups are unique across data sets.
            val = [
                dataset_id + "_".join([str(xxx) for xxx in xx])
                for xx in zip(*[adata.obs[batch_col].values.tolist() for batch_col in batch_cols])
            ]
        else:
            # Locate annotation.
            if hasattr(annotation_container, f"{k}_obs_key") and getattr(annotation_container, f"{k}_obs_key") is not None and \
                    getattr(annotation_container, f"{k}_obs_key") in adata.obs.columns:
                # Last and-clause to check if this column is included in data sets. This may be violated if data
                # is obtained from a database which is not fully streamlined.
                # Look for 1a-* and 1b-I
                k_custom = getattr(annotation_container, f"{k}_obs_key")
                val = adata.obs[k_custom].values.tolist()
                if not (isinstance(keep_orginal_obs, list) and k_custom in keep_orginal_obs):
                    obs_to_delete.append(k_custom)
            else:
                # Look for 2a, 2b
                val = getattr(annotation_container, k)
                if val is None:
                    val = adata_input_ids.unknown_metadata_identifier
                # Unpack nested lists/tuples:
                while hasattr(val, '__len__') and not isinstance(val, str) and len(val) == 1:
                    val = val[0]
                val = [val] * adata.n_obs
        # Identify annotation: disambiguate 1a-I, 1a-II, 1a-III, 1b-I.
        if k in adata_input_ids.ontology_constrained:
            # 1a-*.
            if isinstance(get_ontology(k=k, organism=organism), OntologyHierarchical) and np.all([
                get_ontology(k=k, organism=organism).is_a_node_name(x) or
                x == adata_input_ids.unknown_metadata_identifier
                for x in np.unique(val)
            ]):  # 1a-II)
                new_col = getattr(adata_target_ids, k)
                validation_ontology = get_ontology(k=k, organism=organism)
            elif isinstance(get_ontology(k=k, organism=organism), OntologyHierarchical) and np.all([
                get_ontology(k=k, organism=organism).is_a_node_id(x) or x == adata_input_ids.unknown_metadata_identifier
                for x in np.unique(val)
            ]):  # 1a-III)
                new_col = getattr(adata_target_ids, k) + adata_input_ids.onto_id_suffix
                validation_ontology = None
            else:  # 1a-I)
                new_col = getattr(adata_target_ids, k) + adata_input_ids.onto_original_suffix
                validation_ontology = None
        else:
            # 1b-I.
            new_col = getattr(adata_target_ids, k)
            validation_ontology = get_ontology(k=k, organism=organism)
        # Check values for validity:
        value_protection(attr=new_col, allowed=validation_ontology, attempted=[
            x for x in np.unique(val)
            if x not in [
                adata_input_ids.unknown_metadata_identifier,
            ]
        ], dataset_id=dataset_id)
        obs_new[new_col] = val
        # For ontology-constrained meta data, the remaining columns are added after .obs cleaning below.
    if isinstance(clean_obs, bool) and clean_obs:
        if adata.obsm is not None:
            del adata.obsm
        if adata.obsp is not None:
            del adata.obsp
        adata.obs = obs_new
    elif isinstance(clean_obs, bool) and not clean_obs:
        index_old = adata.obs.index.copy()
        # Add old columns in if they are not duplicated in target obs column space, even if this column is not
        # defined. This would result in the instance accessing this column assuming it was streamlined.
        adata.obs = pd.concat([
            obs_new,
            pd.DataFrame(dict([(k, v) for k, v in adata.obs.items()
                               if (k not in adata_target_ids.controlled_meta_keys and
                                   k not in obs_to_delete)]))
        ], axis=1)
        adata.obs.index = index_old
    elif isinstance(clean_obs, list) or isinstance(clean_obs, tuple) or isinstance(clean_obs, np.ndarray):
        index_old = adata.obs.index.copy()
        # Add old columns in if they are not duplicated in target obs column space, even if this column is not
        # defined. This would result in the instance accessing this column assuming it was streamlined.
        adata.obs = pd.concat([
            obs_new,
            pd.DataFrame(dict([(k, v) for k, v in adata.obs.items()
                               if (k not in adata_target_ids.controlled_meta_keys and
                                   k not in clean_obs and
                                   k not in obs_to_delete)]))
        ], axis=1)
        adata.obs.index = index_old
    else:
        raise ValueError(clean_obs)
    for k in [x for x in adata_target_ids.obs_keys if x in adata_target_ids.ontology_constrained]:
        # Add remaining output columns for ontology-constrained meta data.
        adata = impute_ontology_cols_obs(adata=adata, attr=k, adata_input_ids=adata_input_ids,
                                         adata_target_ids=adata_target_ids, ontology_class_maps=ontology_class_maps,
                                         organism=organism)
        # Delete attribute-specific columns that are not desired.
        col_name = getattr(adata_target_ids, k) + adata_target_ids.onto_id_suffix
        if not keep_id_obs and col_name in adata.obs.columns:
            del adata.obs[col_name]
        col_name = getattr(adata_target_ids, k) + adata_target_ids.onto_original_suffix
        if (isinstance(keep_orginal_obs, list) or not keep_orginal_obs) and col_name in adata.obs.columns:
            del adata.obs[col_name]
        col_name = getattr(adata_target_ids, k)
        if not keep_symbol_obs and col_name in adata.obs.columns:
            del adata.obs[col_name]
    if clean_obs_names:
        adata.obs.index = [f"{dataset_id}_{i}" for i in range(1, adata.n_obs + 1)]

    # Make sure that correct unknown_metadata_identifier is used in .uns, .obs and .var metadata
    unknown_old = adata_input_ids.unknown_metadata_identifier
    unknown_new = adata_target_ids.unknown_metadata_identifier
    adata.obs = adata.obs.replace({None: unknown_new})
    adata.obs = adata.obs.replace({unknown_old: unknown_new})
    adata.var = adata.var.replace({None: unknown_new})
    adata.var = adata.var.replace({unknown_old: unknown_new})
    for k in adata.uns_keys():
        if adata.uns[k] is None or adata.uns[k] == unknown_old:
            adata.uns[k] = unknown_new

    # Add additional hard-coded description changes for cellxgene schema:
    if schema.startswith("cellxgene"):
        schema_version = schema.split(":")[-1] if ":" in schema else None
        adata = cellxgene_export_adaptor(adata=adata,
                                         adata_ids=adata_target_ids,
                                         layer_key_counts=annotation_container.layer_counts,
                                         layer_key_proc=annotation_container.layer_processed,
                                         obs_keys_batch=annotation_container.tech_sample_obs_key,
                                         version=schema_version,
                                         **kwargs)
    return adata


def impute_ontology_cols_obs(adata: anndata.AnnData,
                             adata_input_ids: AdataIds,
                             adata_target_ids: AdataIds,
                             attr: str,
                             ontology_class_maps: Dict[str, pd.DataFrame],
                             organism: str) -> anndata.AnnData:
    """
    Add missing ontology defined columns (symbol, ID, original) for a given ontology.

    1) If original column is non-empty and symbol and ID are empty:
        orginal column is projected to ontology and both symbol and ID are inferred.
        Note that in this case, a label map is required.
    2) If ID column is non-empty or symbol is non-empty, an error is thrown.
        a) If ID column is non-empty and symbol is empty, symbol is inferred.
        b) If ID column is empty and symbol is non-empty, ID is inferred.
        c) If ID column is non-empty and non-symbol is empty, symbol is inferred and over-written.
            Note that this setting allows usage of data sets which were streamlined with a different ontology
            version.
        In all cases original is kept if it is set and is set to symbol otherwise.
    3) If original, ID and symbol columns are empty, no action is taken (meta data item was not set).
    """
    ontology = get_ontology(k=attr, organism=organism)
    col_symbol = getattr(adata_target_ids, attr)
    col_id = getattr(adata_target_ids, attr) + adata_input_ids.onto_id_suffix
    col_original = getattr(adata_target_ids, attr) + adata_input_ids.onto_original_suffix
    if ontology is None:
        # Fill with invalid ontology identifiers if no ontology was found.
        adata.obs[col_id] = \
            [adata_input_ids.invalid_metadata_identifier for _ in range(adata.n_obs)]
        adata.obs[col_original] = \
            [adata_input_ids.invalid_metadata_identifier for _ in range(adata.n_obs)]
        adata.obs[col_symbol] = \
            [adata_input_ids.invalid_metadata_identifier for _ in range(adata.n_obs)]
    else:
        # Note that for symbol and ID, the columns may be filled but not streamlined according to the ontology,
        # in that case the corresponding meta data is defined as absent.
        # Check which level of meta data annotation is present.
        # Symbols:
        symbol_col_present = col_symbol in adata.obs.columns
        symbol_col_streamlined = np.all([
            ontology.is_a_node_name(x) or x == adata_input_ids.unknown_metadata_identifier
            for x in np.unique(adata.obs[col_symbol].values)]) if symbol_col_present else False
        symbol_present = symbol_col_present and symbol_col_streamlined
        # IDs:
        id_col_present = col_id in adata.obs.columns
        id_col_streamlined = np.all([
            ontology.is_a_node_id(x) or x == adata_input_ids.unknown_metadata_identifier
            for x in np.unique(adata.obs[col_id].values)]) if id_col_present else False
        id_present = id_col_present and id_col_streamlined
        # Original annotation (free text):
        original_present = col_original in adata.obs.columns
        if original_present and not symbol_present and not id_present:  # 1)
            adata = project_free_to_ontology(adata=adata, attr=attr, adata_ids=adata_target_ids,
                                             ontology_class_maps=ontology_class_maps, organism=organism)
        if symbol_present or id_present:  # 2)
            if symbol_present and not id_present:  # 2a)
                adata = project_ontology_ids_obs(adata=adata, attr=attr, from_id=False, adata_ids=adata_target_ids,
                                                 organism=organism)
            if not symbol_present and id_present:  # 2b)
                adata = project_ontology_ids_obs(adata=adata, attr=attr, from_id=True, adata_ids=adata_target_ids,
                                                 organism=organism)
            if symbol_present and id_present:  # 2c)
                adata = project_ontology_ids_obs(adata=adata, attr=attr, from_id=True, adata_ids=adata_target_ids,
                                                 organism=organism)
            if not original_present:
                val = adata.obs[col_symbol]
                adata.obs[col_original] = val
    return adata


def project_free_to_ontology(adata: anndata.AnnData,
                             adata_ids: AdataIds,
                             attr: str,
                             ontology_class_maps: Dict[str, pd.DataFrame],
                             organism: str,
                             dataset_id: str = ""):
    """
    Project free text cell type names to ontology based on mapping table.

    ToDo: add ontology ID setting here.
    :param dataset_id: Identifier of dataset used in logging.
    """
    ontology_map = ontology_class_maps[attr]
    col_original = getattr(adata_ids, attr) + adata_ids.onto_original_suffix
    labels_original = adata.obs[col_original].values
    if ontology_map is not None:  # only if this was defined
        labels_mapped = [
            ontology_map[x] if x in ontology_map.keys()
            else x for x in labels_original
        ]
        # Convert unknown celltype placeholders (needs to be hardcoded here as placeholders are also hardcoded in
        # conversion tsv files
        placeholder_conversion = {
            "unknown": adata_ids.unknown_metadata_identifier,
            "NOT_A_CELL": adata_ids.not_a_cell_celltype_identifier,
        }
        labels_mapped = [
            placeholder_conversion[x] if x in placeholder_conversion.keys()
            else x for x in labels_mapped
        ]
        map_exceptions = [adata_ids.unknown_metadata_identifier]
        if attr == "cell_type":
            map_exceptions.append(adata_ids.not_a_cell_celltype_identifier)
        # Validate mapped IDs based on ontology:
        # This aborts with a readable error if there was a target in the mapping file that doesnt match the ontology
        # This protection blocks progression in the unit test if not deactivated.
        value_protection(
            attr=attr,
            allowed=get_ontology(k=attr, organism=organism),
            attempted=[x for x in list(set(labels_mapped)) if x not in map_exceptions],
        )
        # Add cell type IDs into object:
        # The IDs are not read from a source file but inferred based on the class name.
        # TODO this could be changed in the future, this allows this function to be used both on cell type name
        #  mapping files with and without the ID in the third column.
        # This mapping blocks progression in the unit test if not deactivated.
        adata.obs[getattr(adata_ids, attr)] = labels_mapped
        adata = project_ontology_ids_obs(adata=adata, attr=attr, map_exceptions=map_exceptions, from_id=False,
                                         adata_ids=adata_ids, organism=organism)
    else:
        # Assumes that the original labels are the correct ontology symbols, because of a lack of ontology,
        # ontology IDs cannot be inferred.
        # TODO is this necessary in the future?
        adata.obs[getattr(adata_ids, attr)] = labels_original
        adata.obs[getattr(adata_ids, attr) + adata_ids.onto_id_suffix] = \
            [adata_ids.unknown_metadata_identifier] * adata.n_obs
    adata.obs[getattr(adata_ids, attr) + adata_ids.onto_original_suffix] = labels_original
    return adata


def project_ontology_ids_obs(adata: anndata.AnnData,
                             attr: str,
                             adata_ids: AdataIds,
                             organism: str,
                             map_exceptions: Union[None, List[str]] = None,
                             map_exceptions_value=None,
                             from_id: bool = False) -> anndata.AnnData:
    """
    Project ontology names to IDs for a given ontology in .obs entries.

    :param attr: name of obs_column containing names to convert or python list containing these values
    :param map_exceptions: list of values that should not be mapped.
        Defaults to unknown meta data identifier defined in ID object if None.
    :param map_exceptions_value: placeholder target value for values excluded from mapping.
        Defaults to unknown meta data identifier defined in ID object if None.
    :param from_id: Whether to output ontology symbol or ID.
    :return:
    """
    ontology = get_ontology(k=attr, organism=organism)
    assert ontology is not None, f"cannot project value for {attr} because ontology is None"
    assert isinstance(attr, (str, list)), f"argument key_in needs to be of type str or list. Supplied" \
                                          f"type: {type(attr)}"
    map_exceptions = map_exceptions if map_exceptions is not None else [adata_ids.unknown_metadata_identifier]
    map_exceptions = [x.lower() for x in map_exceptions]
    if map_exceptions_value is None:
        # TODO this may be simplified in the future once all unknown meta data labels are the same.
        if attr == "cell_type":
            map_exceptions_value = adata_ids.unknown_metadata_identifier
        else:
            map_exceptions_value = adata_ids.unknown_metadata_identifier
    col_name = getattr(adata_ids, attr)
    if from_id:
        col_name += adata_ids.onto_id_suffix
    input_values = adata.obs[col_name].values
    map_vals = dict([
        (x, ontology.convert_to_name(x)) if from_id else
        (x, ontology.convert_to_id(x))
        for x in np.unique([
            xx for xx in input_values
            if (xx.lower() not in map_exceptions and xx is not None)
        ])
    ])
    output_values = [
        map_vals[x] if x in map_vals.keys() else map_exceptions_value
        for x in input_values
    ]
    key_out = getattr(adata_ids, attr) if from_id else getattr(adata_ids, attr) + adata_ids.onto_id_suffix
    adata.obs[key_out] = output_values
    return adata


def value_protection(allowed: Union[Ontology, None],
                     attempted,
                     attr: str,
                     dataset_id: str = ""):
    """
    Check whether value is from set of allowed values.

    Does not check if allowed is None.
    Cleans entry to term name if ontology ID is provided.

    :param attr: Attribute to set.
    :param allowed: Constraint for values of `attr`.
        Either ontology instance used to constrain entries, or list of allowed values.
    :param attempted: Value(s) to attempt to set in `attr`.
    :param dataset_id: Dataset ID used for logging.
    :return:
    """
    if not isinstance(attempted, list):
        if isinstance(attempted, np.ndarray):
            attempted_ls = attempted.tolist()
        elif isinstance(attempted, tuple):
            attempted_ls = list(attempted)
        else:
            attempted_ls = [attempted]
    else:
        attempted_ls = attempted
    attempted_clean = []
    for x in attempted_ls:
        if allowed is None:
            attempted_clean.append(x)
        elif isinstance(allowed, Ontology):
            if attr == "disease" and (x.lower() == "normal" or x.lower() == "healthy"):
                # TODO required because of missing streamlining between sfaira and 10x, remove in future.
                attempted_clean.append("healthy")
            elif x in allowed.node_names:
                attempted_clean.append(x)
            else:
                if isinstance(allowed, OntologyHierarchical) and x in allowed.node_ids:
                    attempted_clean.append(allowed.convert_to_name(x))
                else:
                    raise ValueError(f"'{x}' is not a valid entry for {attr} in data set {dataset_id}.")
        else:
            raise ValueError(f"argument allowed of type {type(allowed)} is not a valid entry for {attr}.")
    # Flatten attempts if only one was made:
    if len(attempted_clean) == 1:
        attempted_clean = attempted_clean[0]
    return attempted_clean


def get_ontology(k, organism: str) -> Union[OntologyHierarchical, None]:
    # Use global instance of ontology container:
    ocs = OCS

    x = getattr(ocs, k) if hasattr(ocs, k) else None
    if x is not None and isinstance(x, dict):
        assert isinstance(organism, str), organism
        # Check if organism-specific option is available, otherwise choose generic option:
        if organism in x.keys():
            k = organism
        else:
            k = ocs.key_other
            assert k in x.keys(), x.keys()  # Sanity check on dictionary keys.
        x = x[k]
        assert x is None or isinstance(x, Ontology), x  # Sanity check on dictionary element.
    return x
