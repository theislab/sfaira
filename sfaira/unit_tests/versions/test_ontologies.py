from sfaira.versions.metadata import OntologyUberon, OntologyCelltypes, OntologyMmusdv, OntologyHsapdv, \
    OntologyHancestro


def test_cl():
    _ = OntologyCelltypes(branch="v2021-02-01")


def test_uberon():
    _ = OntologyUberon()
