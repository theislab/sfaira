import numpy as np
import pytest
from typing import Tuple, Union

from sfaira.versions.genomes import GenomeContainer

"""
GenomeContainer
"""


@pytest.mark.parametrize("organism", ["mouse"])
@pytest.mark.parametrize("assembly", [None, "Mus_musculus.GRCm38.102"])
def test_gc_init(organism: Union[str, None], assembly: Union[str, None]):
    """
    Tests different modes of initialisation for fatal errors.
    """
    gc = GenomeContainer(organism=organism, assembly=assembly)
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
    gc = GenomeContainer(organism=None, assembly="Mus_musculus.GRCm38.102")
    gc.subset(**subset[0])
    assert gc.n_var == subset[1]
    assert len(gc.ensembl) == subset[1]
    assert len(gc.symbols) == subset[1]
    assert len(gc.biotype) == subset[1]
    if list(subset[0].keys())[0] == "protein_coding":
        assert np.all(gc.biotype == "protein_coding")
