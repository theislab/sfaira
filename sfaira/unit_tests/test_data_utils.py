import pytest

from sfaira.data.utils import map_celltype_to_ontology


@pytest.mark.parametrize("perfect_match", [True, False])
def test_map_celltype_to_ontology(perfect_match: bool):
    trial_cell_type = "T cell" if perfect_match else "T cells"
    x = map_celltype_to_ontology(source=trial_cell_type, organism="human", choices_for_perfect_match=False)
    print(x)
    if perfect_match:
        assert isinstance(x, str)
    else:
        assert isinstance(x, list)
