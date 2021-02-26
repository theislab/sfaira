import pytest
from typing import Union

from sfaira.data.utils import map_celltype_to_ontology


@pytest.mark.parametrize("trial_cell_type_labels",
                         ["type B pancreatic cell", "beta", ["type B pancreatic cell", "beta"]])
@pytest.mark.parametrize("choices_for_perfect_match", [True, False])
@pytest.mark.parametrize("anatomical_constraint", [None, "pancreas"])
def test_map_celltype_to_ontology(
        trial_cell_type_labels: str,
        choices_for_perfect_match: bool,
        anatomical_constraint: Union[str, None]
):
    trial_cell_type_labels = [trial_cell_type_labels] if isinstance(trial_cell_type_labels, str) \
        else trial_cell_type_labels
    perfectly_matched_query = ["type B pancreatic cell" == x for x in trial_cell_type_labels]
    matches = map_celltype_to_ontology(
        queries=trial_cell_type_labels,
        organism="human",
        include_synonyms=True,
        anatomical_constraint=anatomical_constraint,
        choices_for_perfect_match=choices_for_perfect_match,
        always_return_dict=False,
    )
    for x, y in zip(trial_cell_type_labels, perfectly_matched_query):
        if isinstance(matches, dict):  # dictionary over queries with list of matches as value each
            if y and not choices_for_perfect_match:
                assert isinstance(matches[x], str), matches
                assert matches[x] == "type B pancreatic cell"
            else:
                assert isinstance(matches[x], list), matches
                assert "type B pancreatic cell" in matches[x]
        else:  # matches for single query
            if y and not choices_for_perfect_match:
                assert isinstance(matches, str), matches
                assert matches == "type B pancreatic cell"
            else:
                assert isinstance(matches, list), matches
                assert "type B pancreatic cell" in matches
