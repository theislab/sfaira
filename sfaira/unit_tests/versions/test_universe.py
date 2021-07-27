import os

from sfaira.versions.metadata import CelltypeUniverse, OntologyCl, OntologyUberon

dir_temp = os.path.join(os.path.dirname(__file__), "temp")

"""
CelltypeUniverse
"""


def test_universe_io():
    if not os.path.exists(dir_temp):
        os.mkdir(dir_temp)
    tmp_fn = os.path.join(dir_temp, "universe_temp.csv")
    targets = ["stromal cell", "lymphocyte", "T-helper 1 cell", "T-helper 17 cell"]
    leaves_target = ["stromal cell", "T-helper 1 cell", "T-helper 17 cell"]
    cl = OntologyCl(branch="v2021-02-01")
    uberon = OntologyUberon()
    cu = CelltypeUniverse(cl=cl, uberon=uberon)
    cu.write_target_universe(fn=tmp_fn, x=targets)
    cu.load_target_universe(fn=tmp_fn)
    os.remove(tmp_fn)
    leaves = cu.onto_cl.convert_to_name(cu.onto_cl.leaves)
    assert set(leaves) == set(leaves_target), (leaves, leaves_target)
