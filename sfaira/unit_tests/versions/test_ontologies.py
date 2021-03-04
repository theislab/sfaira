from sfaira.versions.metadata import OntologyUberon, OntologyCelltypes, OntologyMmusdv, OntologyHsapdv, \
    OntologyHancestro


def test_cl():
    _ = OntologyCelltypes(branch="v2021-02-01")


def test_uberon():
    _ = OntologyUberon()


def test_uberon_subsetting():
    ou = OntologyUberon()
    assert ou.is_a(query="lung", reference="lobe of lung")
    assert ou.is_a(query="lung", reference="lobe of lung")
    assert not ou.is_a(query="lobe of lung", reference="lung")

    assert ou.is_a(query="lung", reference="lobar bronchus")
    assert ou.is_a(query="lung", reference="lobar bronchus")
    assert not ou.is_a(query="lobar bronchus", reference="lung")
