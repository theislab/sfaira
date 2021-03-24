from sfaira.versions.metadata import OntologyUberon, OntologyCelltypes, OntologyMondo, OntologyMmusdv, OntologyHsapdv

"""
CL
"""


def test_cl_loading():
    _ = OntologyCelltypes(branch="v2021-02-01")


def test_cl_subsetting():
    oc = OntologyCelltypes(branch="v2021-02-01")
    assert oc.is_a(query="T cell", reference="lymphocyte")
    assert oc.is_a(query="lymphocyte", reference="lymphocyte")
    assert not oc.is_a(query="lymphocyte", reference="T cell")


"""
Hancestro
"""


#def test_hancestro_loading():
#    _ = OntologyHancestro()


"""
Hsapdv
"""


def test_hsapdv_loading():
    _ = OntologyHsapdv()


"""
MONDO
"""


def test_mondo_loading():
    _ = OntologyMondo()


"""
Mmusdv
"""


def test_mmusdv_loading():
    _ = OntologyMmusdv()


"""
UBERON
"""


def test_uberon_loading():
    _ = OntologyUberon()


def test_uberon_subsetting():
    ou = OntologyUberon()
    assert ou.is_a(query="lobe of lung", reference="lung")
    assert ou.is_a(query="lobe of lung", reference="lobe of lung")
    assert not ou.is_a(query="lung", reference="lobe of lung")

    assert ou.is_a(query="lobar bronchus", reference="lung")
    assert ou.is_a(query="lobar bronchus", reference="lobar bronchus")
    assert not ou.is_a(query="lung", reference="lobar bronchus")

