import numpy as np
import pytest
from typing import Tuple, Union

from sfaira.versions.genomes import GenomeContainer, translate_id_to_symbols, translate_symbols_to_id

ASSEMBLY = "Mus_musculus.GRCm38.102"

"""
GenomeContainer.
"""


@pytest.mark.parametrize("assembly", [ASSEMBLY])
def test_gc_init(assembly: Union[str]):
    """
    Tests different modes of initialisation for fatal errors.
    """
    gc = GenomeContainer(assembly=assembly)
    assert gc.organism == "mus_musculus"


@pytest.mark.parametrize("subset", [
    ({"biotype": "protein_coding"}, 21936),
    ({"biotype": "lincRNA"}, 5629),
    ({"biotype": "protein_coding,lincRNA"}, 21936 + 5629),
    ({"symbols": "Gnai3,Pbsn,Cdc45"}, 3),
    ({"ensg": "ENSMUSG00000000003,ENSMUSG00000000028"}, 2)
])
def test_gc_subsetting(subset: Tuple[dict, int]):
    """
    Tests if genome container is subsetted correctly.
    """
    gc = GenomeContainer(assembly="Mus_musculus.GRCm38.102")
    gc.subset(**subset[0])
    assert gc.n_var == subset[1]
    assert len(gc.ensembl) == subset[1]
    assert len(gc.symbols) == subset[1]
    assert len(gc.biotype) == subset[1]
    if list(subset[0].keys())[0] == "protein_coding":
        assert np.all(gc.biotype == "protein_coding")


"""
Utils.
"""


@pytest.mark.parametrize("genes", [
    ("Adora3", "ENSMUSG00000000562"),  # single string
    (["Adora3", "Timp1"], ["ENSMUSG00000000562", "ENSMUSG00000001131"]),  # list of strings
    (["ADORA3", "timp1"], ["EnsmusG00000000562", "ENSMUSG00000001131"]),  # list of strings with weird capitalization
])
def test_translate_id_to_symbols(genes):
    """
    Tests translate_id_to_symbols and translate_symbols_to_id for translation errors.
    """
    x, y = genes
    y_hat = translate_symbols_to_id(x=x, assembly="Mus_musculus.GRCm38.102")
    # Correct target spelling of y:
    y = [z.upper() for z in y] if isinstance(y, list) else y.upper()
    assert np.all(y_hat == y)
    y, x = genes
    y_hat = translate_id_to_symbols(x=x, assembly="Mus_musculus.GRCm38.102")
    # Correct target spelling of y:
    y = [z[0].upper() + z[1:].lower() for z in y] if isinstance(y, list) else y[0].upper() + y[1:].lower()
    assert np.all(y_hat == y)
