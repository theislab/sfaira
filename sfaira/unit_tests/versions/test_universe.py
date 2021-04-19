import os

from sfaira.versions.metadata import CelltypeUniverse, OntologyCl, OntologyUberon

"""
CelltypeUniverse
"""


def test_universe_io():
    tmp_fn = "./universe_tempp.csv"
    targets = ["stromal cell", "lymphocyte", "T-helper 1 cell", "T-helper 17 cell"]
    cl = OntologyCl(branch="v2021-02-01")
    uberon = OntologyUberon()
    cu = CelltypeUniverse(cl=cl, uberon=uberon)
    cu.write_target_universe(fn=tmp_fn, x=targets)
    cu.load_target_universe(fn=tmp_fn)
    os.remove(tmp_fn)
    leaves = cu.leaves
    assert set(leaves) == set(targets), (leaves, targets)
