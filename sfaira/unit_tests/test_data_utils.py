import pytest
from typing import Union

from sfaira.data.utils import map_celltype_to_ontology


@pytest.mark.parametrize("perfectly_matched_query", [True, False])
@pytest.mark.parametrize("choices_for_perfect_match", [True, False])
@pytest.mark.parametrize("anatomical_constraint", [None, "pancreas"])
def test_map_celltype_to_ontology(
        perfectly_matched_query: bool,
        choices_for_perfect_match: bool,
        anatomical_constraint: Union[str, None]
):
    trial_cell_type = "type B pancreatic cell" if perfectly_matched_query else "beta"
    x = map_celltype_to_ontology(
        queries=[trial_cell_type],
        organism="human",
        include_synonyms=True,
        anatomical_constraint=anatomical_constraint,
        choices_for_perfect_match=choices_for_perfect_match
    )
    if perfectly_matched_query and (not choices_for_perfect_match or not anatomical_constraint):
        assert isinstance(x[trial_cell_type], str), x
        assert x[trial_cell_type] == "type B pancreatic cell"
    else:
        assert isinstance(x[trial_cell_type], list), x
        assert "type B pancreatic cell" in x[trial_cell_type]
