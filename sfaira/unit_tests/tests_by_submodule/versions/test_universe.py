import os

from sfaira.consts.schema import DEFAULT_SCHEMA, ONTOLOGY_VERSIONS
from sfaira.versions.metadata import CelltypeUniverse, OntologyCl, OntologyUberon
from sfaira.unit_tests import DIR_TEMP

"""
CelltypeUniverse
"""


def test_universe_io():
    if not os.path.exists(DIR_TEMP):
        os.mkdir(DIR_TEMP)
    tmp_fn = os.path.join(DIR_TEMP, "universe_temp.csv")
    targets = ["stromal cell", "lymphocyte", "T-helper 1 cell", "T-helper 17 cell"]
    leaves_target = ["stromal cell", "T-helper 1 cell", "T-helper 17 cell"]
    cl = OntologyCl(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_CL"])
    uberon = OntologyUberon(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_UBERON"])
    cu = CelltypeUniverse(cl=cl, uberon=uberon)
    cu.write_target_universe(fn=tmp_fn, x=targets)
    cu.load_target_universe(fn=tmp_fn)
    os.remove(tmp_fn)
    leaves = cu.onto_cl.convert_to_name(cu.onto_cl.leaves)
    assert set(leaves) == set(leaves_target), (leaves, leaves_target)
