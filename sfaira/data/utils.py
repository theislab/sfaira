from typing import Dict, List, Union

from sfaira.versions.metadata import CelltypeUniverse


def map_celltype_to_ontology(
        source: str,
        organism: str,
        include_synonyms: bool = True,
        anatomical_constraint: Union[str, None] = None,
        n_suggest: int = 4,
        choices_for_perfect_match: bool = True,
        keep_strategy: bool = False,
        **kwargs
) -> Union[List[str], Dict[str, List[str]], str]:
    """
    Map free text node name to ontology node names via sfaira cell type matching.

    For details, see also sfaira.versions.metadata.CelltypeUniverse.prepare_celltype_map_fuzzy()

    :param source: Free text node label which is to be matched to ontology nodes. Must not be a a list or tuple.
    :param organism: Organism, defines ontology extension used.
    :param include_synonyms: Whether to include synonyms of nodes in string search.
    :param anatomical_constraint: Whether to require suggestions to be within a target anatomy defined within UBERON.
    :param n_suggest: Number of cell types to suggest per search strategy.
    :param choices_for_perfect_match: Whether to give additional matches if a perfect match was found. Note that there
        are cases in which an apparent perfect match corresponds to a general term which could be specified knowing the
        anatomic location of the sample. If this is False and a perfect match is found, only this perfect match is
        returned as a string, rather than as a list.
    :param keep_strategy: Whether to keep search results structured by search strategy.
        For details, see also sfaira.versions.metadata.CelltypeUniverse.prepare_celltype_map_fuzzy()
    :param **kwargs: Additional parameters to CelltypeUniverse.
    :return: List of high priority matches or perfect match (see choices_for_perfect_match) or, if keep_strategy,
        dictionary of lists of search strategies named by strategy name.
    """
    assert isinstance(source, str)
    cu = CelltypeUniverse(organism=organism, **kwargs)
    matches = cu.prepare_celltype_map_fuzzy(
        source=[source],
        match_only=False,
        include_synonyms=include_synonyms,
        anatomical_constraint=anatomical_constraint,
        omit_list=[],
        n_suggest=n_suggest,
    )[0][0]
    # Flatten list of lists:
    matches_flat = [x for xx in matches.values() for x in xx]
    if not choices_for_perfect_match and source in matches_flat:
        return source
    elif keep_strategy:
        return matches
    else:
        return matches_flat
