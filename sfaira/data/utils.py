import anndata
import numpy as np
import pandas as pd
import scipy.sparse
from typing import Dict, List, Union
import yaml


from sfaira.versions.metadata.maps import prepare_ontology_map


def map_freetext_to_ontology(
        queries: Union[str, List[str]],
        onto: str,
        always_return_dict: bool = False,
        anatomical_constraint: Union[str, None] = None,
        choices_for_perfect_match: bool = False,
        include_synonyms: bool = True,
        keep_strategy: bool = False,
        n_suggest: int = 4,
        omit_target_list: list = [],
        threshold_for_partial_matching: float = 90.,
        **kwargs
) -> Union[List[str], Dict[str, List[str]], str]:
    """
    Map free text node name to ontology node names via sfaira cell type matching.

    For details, see also sfaira.versions.metadata.CelltypeUniverse.prepare_celltype_map_fuzzy()

    :param queries: Free text node label which is to be matched to ontology nodes.
        Can also be a list of strings to query.
    :param onto: Name of ontology to map to.
    :param always_return_dict: Also return a dictionary over queries if only one query was given.
    :param anatomical_constraint: Whether to require suggestions to be within a target anatomy defined within UBERON.
    :param choices_for_perfect_match: Whether to give additional matches if a perfect match was found. Note that there
        are cases in which an apparent perfect match corresponds to a general term which could be specified knowing the
        anatomic location of the sample. If this is False and a perfect match is found, only this perfect match is
        returned as a string, rather than as a list.
    :param include_synonyms: Whether to include synonyms of nodes in string search.
    :param keep_strategy: Whether to keep search results structured by search strategy.
        For details, see also sfaira.versions.metadata.CelltypeUniverse.prepare_celltype_map_fuzzy()
    :param n_suggest: Number of cell types to suggest per search strategy.
    :param omit_target_list: Ontology nodes to not match to.
    :param threshold_for_partial_matching: Maximum fuzzy match score below which lenient matching (ratio) is
        extended through partial_ratio.
    :param kwargs: Additional parameters to CelltypeUniverse.
    :return: List over queries, each entry is:
        A list of high priority matches or perfect match (see choices_for_perfect_match) or, if keep_strategy,
        dictionary of lists of search strategies named by strategy name. If a search strategy yields perfect matches, it
        does not return a list of strings but just a single string.
    """
    if isinstance(queries, str):
        queries = [queries]
    matches_to_return = {}
    matches = prepare_ontology_map(
        source=queries,
        onto=onto,
        match_only=False,
        include_synonyms=include_synonyms,
        anatomical_constraint=anatomical_constraint,
        choices_for_perfect_match=choices_for_perfect_match,
        omit_list=[],
        omit_target_list=omit_target_list,
        n_suggest=n_suggest,
        threshold_for_partial_matching=threshold_for_partial_matching,
    )
    matches = matches[0]
    # Prepare the output:
    for x, matches_i in zip(queries, matches):
        matches_i = matches_i
        # Flatten list of lists:
        # Flatten dictionary of lists and account for string rather than list entries.
        if len(matches_i.values()) == 1 and isinstance(list(matches_i.values())[0], str):
            matches_flat = list(matches_i.values())[0]
        else:
            matches_flat = []
            for xx in matches_i.values():
                if isinstance(xx, list):
                    matches_flat.extend(xx)
                else:
                    assert isinstance(xx, str)
                    matches_flat.append(xx)
        if not choices_for_perfect_match and x in matches_flat:
            matches_to_return.update({x: x})
        elif keep_strategy:
            matches_to_return.update({x: matches_i})
        else:
            matches_to_return.update({x: matches_flat})
    # Only return a list over queries if more than one query was given.
    if len(queries) == 1 and not always_return_dict:
        return matches_to_return[queries[0]]
    else:
        return matches_to_return


def read_yaml(fn) -> Dict[str, Dict[str, Union[str, int, bool]]]:
    """
    Read data yaml file.

    Matches format names to Dataset attribute names.

    :param fn: YAML file name.
    :return: Dictionary of dictionaries of names of Dataset attributes and their values.

        - "attr": Data set attributes.
        - "meta": Meta information of yaml and representation.
    """
    with open(fn) as f:
        yaml_dict = yaml.safe_load(f)
    attr_dict = {}
    meta_dict = {}
    for k, v in yaml_dict.items():
        if k not in ["dataset_structure", "meta"]:
            attr_dict.update(v)
        else:
            meta_dict.update(v)
    return {"attr": attr_dict, "meta": meta_dict}


def collapse_matrix(adata: Union[anndata.AnnData, anndata.Raw], var_column: str) -> anndata.AnnData:
    """
    Collapses (sum) features with the same var_name in a provided var column.

    Does not retain .varm if duplicated var_names are found.
    keeps .var column of first occurrence of duplicated variables.

    :param adata: Input anndata instance with potential duplicated var names.
    :param var_column: column name in .var that contains the duplicated features of interest
    :return: Processed anndata instance without duplicated var names.
    """
    old_index = adata.var.index.tolist() if var_column == "index" else adata.var[var_column].tolist()
    new_index = list(np.unique(old_index))
    if len(new_index) < len(old_index):
        idx_map = np.array([np.where(x == np.array(old_index))[0] for x in new_index])
        # Build initial matrix from first match.
        data = adata.X[:, np.array([x[0] for x in idx_map])].copy()
        # Add additional matched (duplicates) on top:
        for i, idx in enumerate(idx_map):
            if len(idx) > 1:
                data[:, i] = data[:, i] + adata.X[:, idx[1:]].sum(axis=1)

        # Remove varm and populate var with first occurrence only:
        obs_names = adata.obs_names
        if isinstance(adata, anndata.AnnData):
            adata = anndata.AnnData(
                X=data,
                obs=adata.obs,
                obsm=adata.obsm,
                var=adata.var.iloc[[old_index.index(x) for x in new_index]],
                uns=adata.uns
            )
        elif isinstance(adata, anndata.Raw):
            adata = anndata.AnnData(
                X=data,
                var=adata.var.iloc[[old_index.index(x) for x in new_index]],
            )
        else:
            assert False
        adata.obs_names = obs_names
    return adata


def subset_adata_genes(x, var, feature_id_var_key, feature_symbol_var_key, subset_ids_ensg, subset_ids_symbol):
    # Convert data matrix to csc matrix
    if isinstance(x, np.ndarray):
        # Change NaN to zero. This occurs for example in concatenation of anndata instances.
        if np.any(np.isnan(x)):
            x[np.isnan(x)] = 0
        x = scipy.sparse.csc_matrix(x)
    elif isinstance(x, scipy.sparse.spmatrix):
        x = x.tocsc()
    else:
        raise ValueError(f"Data type {type(x)} not recognized.")

    # Remove unmapped genes
    data_ids_ensg = var.index.values if feature_id_var_key == "index" else var[feature_id_var_key].values
    idx_feature_kept = np.where([x.upper() in subset_ids_ensg for x in data_ids_ensg])[0]
    data_ids_kept = data_ids_ensg[idx_feature_kept]
    x = x[:, idx_feature_kept]
    # Build map of subset_ids to features in x:
    idx_feature_map = np.array([subset_ids_ensg.index(x) for x in data_ids_kept])
    # Create reordered feature matrix based on reference and convert to csr
    x_new = scipy.sparse.csc_matrix((x.shape[0], len(subset_ids_ensg)), dtype=x.dtype)
    # copying this over to the new matrix in chunks of size `steps` prevents a strange scipy error:
    # ... scipy/sparse/compressed.py", line 922, in _zero_many i, j, offsets)
    # ValueError: could not convert integer scalar
    step = 500
    if step < len(idx_feature_map):
        i = 0
        for i in range(0, len(idx_feature_map), step):
            x_new[:, idx_feature_map[i:i + step]] = x[:, i:i + step]
        x_new[:, idx_feature_map[i + step:]] = x[:, i + step:]
    else:
        x_new[:, idx_feature_map] = x
    x_new = x_new.tocsr()

    # Create new var dataframe
    if feature_symbol_var_key == "index":
        var_index = subset_ids_symbol
        var_data = {feature_id_var_key: subset_ids_ensg}
    elif feature_id_var_key == "index":
        var_index = subset_ids_ensg
        var_data = {feature_symbol_var_key: subset_ids_symbol}
    else:
        var_index = None
        var_data = {feature_symbol_var_key: subset_ids_symbol,
                    feature_id_var_key: subset_ids_ensg}
    var_new = pd.DataFrame(data=var_data, index=var_index)
    return x_new, var_new
