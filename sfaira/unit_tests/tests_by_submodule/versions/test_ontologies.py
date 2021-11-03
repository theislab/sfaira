import numpy as np
from sfaira.consts.ontologies import DEFAULT_CL, DEFAULT_HSAPDV, DEFAULT_MONDO, DEFAULT_MMUSDV, DEFAULT_PATO, \
    DEFAULT_NCBITAXON, DEFAULT_UBERON
from sfaira.versions.metadata import OntologyUberon, OntologyCl, OntologyHancestro, OntologyHsapdv, OntologyMondo, \
    OntologyMmusdv, OntologyEfo, OntologyTaxon, OntologySex, OntologyUberonLifecyclestage

"""
OntologyCelltypes
"""


def test_cl_loading():
    """
    Tests if ontology can be initialised.
    """
    _ = OntologyCl(branch=DEFAULT_CL, recache=True)
    _ = OntologyCl(branch=DEFAULT_CL, recache=False)


def test_cl_is_a():
    """
    Tests if is-a relationships work correctly.
    """
    oc = OntologyCl(branch=DEFAULT_CL)
    assert oc.is_a(query="T cell", reference="lymphocyte")
    assert oc.is_a(query="lymphocyte", reference="lymphocyte")
    assert not oc.is_a(query="lymphocyte", reference="T cell")


def test_effective_leaves():
    """
    Tests if node sets can be mapped to effective leaf sets via `OntologyCelltypes.get_effective_leaves()`
    """
    oc = OntologyCl(branch=DEFAULT_CL)
    x = oc.get_effective_leaves(x=[
        "CD4-positive helper T cell", "lymphocyte", "stromal cell", "T cell", "T-helper 1 cell",
        "T-helper 17 cell"
    ])
    x = oc.convert_to_name(x)
    assert set(x) == {"stromal cell", "T-helper 1 cell", "T-helper 17 cell"}, x


def test_map_leaves():
    """
    Tests if nodes can be mapped to leave nodes in ontology.
    """
    oc = OntologyCl(branch=DEFAULT_CL)
    leaf_map_1 = oc.convert_to_name(oc.map_to_leaves(node="CD4-positive helper T cell", include_self=True))
    leaf_map_2 = oc.map_to_leaves(node="CD4-positive helper T cell", include_self=True, return_type="idx")
    assert len(leaf_map_1) == 7
    assert np.all(leaf_map_2 == np.sort([oc.convert_to_name(oc.leaves).index(x) for x in list(leaf_map_1)]))


def test_set_leaves():
    """
    Tests if ontology behaves correctly if leaf nodes were reset.
    """
    oc = OntologyCl(branch=DEFAULT_CL, use_developmental_relationships=False)
    targets = ["stromal cell", "T-helper 1 cell", "T-helper 17 cell"]
    oc.leaves = targets
    leaves = oc.convert_to_name(oc.leaves)
    assert set(leaves) == set(targets), leaves
    assert len(oc.node_ids) == 22
    assert np.all([x in oc.convert_to_name(oc.node_ids) for x in targets]), oc.convert_to_name(oc.node_ids)
    leaf_map_1 = oc.convert_to_name(oc.map_to_leaves(node="lymphocyte"))
    leaf_map_2 = oc.map_to_leaves(node="lymphocyte", return_type="idx")
    leaf_map_3 = oc.convert_to_name(oc.map_to_leaves(node="T-helper 1 cell"))
    leaf_map_4 = oc.map_to_leaves(node="T-helper 1 cell", return_type="idx")
    assert set(leaf_map_1) == {"T-helper 1 cell", "T-helper 17 cell"}
    assert np.all(leaf_map_2 == np.sort([oc.convert_to_name(oc.leaves).index(x) for x in list(leaf_map_1)]))
    assert set(leaf_map_3) == {"T-helper 1 cell"}
    assert np.all(leaf_map_4 == np.sort([oc.convert_to_name(oc.leaves).index(x) for x in list(leaf_map_3)]))


def test_reset_root():
    """
    Tests if root can be reset correctly.
    """
    oc = OntologyCl(branch=DEFAULT_CL, use_developmental_relationships=False)
    oc.reset_root(root="T cell")
    assert "T-helper 1 cell" in oc.node_names
    assert "T cell" in oc.node_names
    assert "lymphocyte" not in oc.node_names


"""
OntologyEfo
"""


def test_sclc_loading():
    """
    Tests if ontology can be initialised.
    """
    _ = OntologyEfo(recache=True)
    _ = OntologyEfo(recache=False)


def test_sclc_nodes():
    """
    Tests for presence and absence of a few commonly mistaken nodes.
    """
    sclc = OntologyEfo()
    assert "10x technology" in sclc.node_names
    assert "10x 3' v3" in sclc.node_names
    assert "Smart-like" in sclc.node_names
    assert "Smart-seq2" in sclc.node_names
    assert "sci-Plex" in sclc.node_names
    assert "single cell library construction" in sclc.node_names


def test_sclc_is_a():
    """
    Tests if is-a relationships work correctly.
    """
    sclc = OntologyEfo()
    assert sclc.is_a(query="10x 3' v3", reference="10x technology")
    assert sclc.is_a(query="10x 3' v3", reference="10x 3' transcription profiling")
    assert not sclc.is_a(query="10x technology", reference="10x 3' transcription profiling")
    assert sclc.is_a(query="10x 3' v3", reference="single cell library construction")
    assert sclc.is_a(query="sci-Plex", reference="single cell library construction")
    assert not sclc.is_a(query="sci-Plex", reference="10x technology")


"""
Hancestro
"""


def test_hancestro_loading():
    _ = OntologyHancestro(recache=True)
    _ = OntologyHancestro(recache=False)


"""
Hsapdv
"""


def test_hsapdv_loading():
    _ = OntologyHsapdv(branch=DEFAULT_HSAPDV, recache=True)
    _ = OntologyHsapdv(branch=DEFAULT_HSAPDV, recache=False)


"""
MONDO
"""


def test_mondo_loading():
    _ = OntologyMondo(branch=DEFAULT_MONDO, recache=True)
    _ = OntologyMondo(branch=DEFAULT_MONDO, recache=False)


"""
Mmusdv
"""


def test_mmusdv_loading():
    _ = OntologyMmusdv(branch=DEFAULT_MMUSDV, recache=True)
    _ = OntologyMmusdv(branch=DEFAULT_MMUSDV, recache=False)


"""
NCBI Taxon
"""


def test_taxon_loading():
    _ = OntologyTaxon(branch=DEFAULT_NCBITAXON, recache=True)
    _ = OntologyTaxon(branch=DEFAULT_NCBITAXON, recache=False)


"""
Sex
"""


def test_sex_loading():
    _ = OntologySex(branch=DEFAULT_PATO, recache=True)
    _ = OntologySex(branch=DEFAULT_PATO, recache=False)


"""
UBERON
"""


def test_uberon_loading():
    _ = OntologyUberon(branch=DEFAULT_UBERON, recache=True)
    _ = OntologyUberon(branch=DEFAULT_UBERON, recache=False)


def test_uberon_lcs_loading():
    _ = OntologyUberonLifecyclestage(branch=DEFAULT_UBERON, recache=False)


def test_uberon_subsetting():
    ou = OntologyUberon(branch=DEFAULT_UBERON)
    assert ou.is_a(query="lobe of lung", reference="lung")
    assert ou.is_a(query="lobe of lung", reference="lobe of lung")
    assert not ou.is_a(query="lung", reference="lobe of lung")

    assert ou.is_a(query="lobar bronchus", reference="lung")
    assert ou.is_a(query="lobar bronchus", reference="lobar bronchus")
    assert not ou.is_a(query="lung", reference="lobar bronchus")

    assert ou.is_a(query="adipose tissue of abdominal region", reference="adipose tissue")
    assert ou.is_a(query="adipose tissue", reference="adipose tissue")
    assert not ou.is_a(query="adipose tissue", reference="adipose tissue of abdominal region")
