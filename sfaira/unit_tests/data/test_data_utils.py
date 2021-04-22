import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse
from typing import Union

from sfaira.data.utils import map_celltype_to_ontology, collapse_matrix


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


@pytest.mark.parametrize("data", ["scipy.sparse.lil_matrix", "scipy.sparse.csr_matrix", "numpy"])
@pytest.mark.parametrize("duplications", [False, True])
def test_collapse_matrix(
        data: str,
        duplications: bool,
):
    """
    Tests collapse_matrix.

    Tests if:

        - matrix type is maintained
        - row sums are maintained
        - feature dimension is correct

    :param data: Data format.
    :param duplications: Whether feature names are duplicated.
    :return:
    """
    x = np.asarray(np.random.randint(0, 100, size=(10, 10)), dtype="float32")
    if data == "scipy.sparse.lil_matrix":
        x = scipy.sparse.lil_matrix(x)
    elif data == "scipy.sparse.csr_matrix":
        x = scipy.sparse.csr_matrix(x)
    elif data == "numpy":
        pass
    else:
        assert False
    if duplications:
        # Create triplicate and duplicate gene names:
        index = ["g" + str(i) for i in range(2)] + ["g" + str(i) for i in range(3)] + \
                ["g" + str(i) for i in range(x.shape[1] - 3 - 2)]
    else:
        index = ["g" + str(i) for i in range(x.shape[1])]
    adata = anndata.AnnData(x, var=pd.DataFrame({"var_column": index}))
    adata2 = collapse_matrix(adata=adata, var_column="var_column")
    assert adata.X.shape[0] == adata2.X.shape[0], "observation dimension mismatch"
    assert adata.X.dtype == adata2.X.dtype, "type mismatch"
    assert adata2.X.shape[1] == len(np.unique(adata.var["var_column"])), "feature dimension mismatch"
    assert np.all(np.asarray(adata.X.sum()).flatten() == np.asarray(adata2.X.sum().flatten())), \
        "total count mismatch"
    assert np.all(np.asarray(adata.X.sum(axis=1)).flatten() == np.asarray(adata2.X.sum(axis=1).flatten())), \
        "observation-wise count mismatch"
