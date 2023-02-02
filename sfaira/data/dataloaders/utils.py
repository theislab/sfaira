import gzip
import os
import shutil
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
        omit_prefix_list=[],
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


def buffered_decompress(fn_compressed):
    if fn_compressed.endswith(".zip"):
        dir_decompressed = fn_compressed.split(".zip")[0]
        decompressor = "shututil"
    elif fn_compressed.endswith(".tar.gz"):
        dir_decompressed = fn_compressed.split(".tar.gz")[0]
        decompressor = "shututil"
    elif fn_compressed.endswith(".tar"):
        dir_decompressed = fn_compressed.split(".tar")[0]
        decompressor = "shututil"
    elif fn_compressed.endswith(".gz"):
        dir_decompressed = fn_compressed.split(".gz")[0]
        decompressor = "gzip"
    else:
        dir_decompressed = os.path.dirname(fn_compressed)
        decompressor = "shututil"
    if not os.path.exists(dir_decompressed):
        if decompressor == "shututil":
            shutil.unpack_archive(filename=fn_compressed, extract_dir=dir_decompressed)
        elif decompressor == "gzip":
            with gzip.open(fn_compressed, 'rb') as f_in:
                with open(dir_decompressed, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            assert False, decompressor
    return dir_decompressed
