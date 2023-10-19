import numpy as np
from sfaira.consts.schema import DEFAULT_SCHEMA, ONTOLOGY_VERSIONS
from sfaira.versions.metadata import OntologyUberon, OntologyCl, OntologyHancestro, OntologyHsapdv, OntologyList, \
    OntologyMondo, OntologyMmusdv, OntologyEfo, OntologyTaxon, OntologyPato, OntologyUberonLifecyclestage


"""
Base
"""


def test_ontologylist():
    """
    Tests if ontology can be initialised.
    """
    cl = OntologyCl(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_CL"], recache=False)
    # Check list mode:
    cl_list0 = OntologyList(terms=cl.node_ids)
    # IDs and names should be the same:
    assert np.all(np.asarray(cl_list0.node_ids) == np.asarray(cl_list0.node_names))
    assert np.all([len(cl_list0.get_ancestors(node=x)) == 0 for x in cl_list0.node_names])
    assert np.all([x in cl.node_ids for x in cl_list0.node_ids])
    assert np.all([x in cl.node_ids for x in cl_list0.node_names])
    # Check dict mode:
    cl_list1 = OntologyList(terms=cl.nodes_dict)
    # IDs and names should not be the same:
    assert np.all(np.asarray(cl_list1.node_ids) != np.asarray(cl_list1.node_names))
    assert np.all([len(cl_list1.get_ancestors(node=x)) == 0 for x in cl_list1.node_names])
    assert np.all([x in cl.node_ids for x in cl_list1.node_ids])
    assert np.all([x in cl.node_names for x in cl_list1.node_names])


"""
OntologyCelltypes
"""


def test_cl_loading():
    """
    Tests if ontology can be initialised.
    """
    _ = OntologyCl(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_CL"], recache=True)
    _ = OntologyCl(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_CL"], recache=False)


def test_cl_is_a():
    """
    Tests if is-a relationships work correctly.
    """
    oc = OntologyCl(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_CL"])
    assert oc.is_a(is_="T cell", a_="lymphocyte")
    assert oc.is_a(is_="lymphocyte", a_="lymphocyte")
    assert not oc.is_a(is_="lymphocyte", a_="T cell")


def test_effective_leaves():
    """
    Tests if node sets can be mapped to effective leaf sets via `OntologyCelltypes.get_effective_leaves()`
    """
    oc = OntologyCl(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_CL"])
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
    oc = OntologyCl(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_CL"])
    leaf_map_1 = oc.convert_to_name(oc.map_to_leaves(node="CD4-positive helper T cell", include_self=True))
    leaf_map_2 = oc.map_to_leaves(node="CD4-positive helper T cell", include_self=True, return_type="idx")
    assert len(leaf_map_1) == 7
    assert np.all(leaf_map_2 == np.sort([oc.convert_to_name(oc.leaves).index(x) for x in list(leaf_map_1)]))


def test_set_leaves():
    """
    Tests if ontology behaves correctly if leaf nodes were reset.
    """
    oc = OntologyCl(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_CL"])
    targets = ["stromal cell", "T-helper 1 cell", "T-helper 17 cell"]
    oc.leaves = targets
    leaves = oc.convert_to_name(oc.leaves)
    assert set(leaves) == set(targets), leaves
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
    oc = OntologyCl(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_CL"])
    assert len(oc.node_ids) == len(oc.node_names)
    n0 = len(oc.node_ids)
    oc.reset_root(root="T cell")
    assert len(oc.node_ids) == len(oc.node_names)
    n1 = len(oc.node_ids)
    assert n1 < n0
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
    _ = OntologyEfo(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_EFO"], recache=True)
    _ = OntologyEfo(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_EFO"], recache=False)


def test_sclc_nodes():
    """
    Tests for presence and absence of a few commonly mistaken nodes.
    """
    sclc = OntologyEfo(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_EFO"])
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
    sclc = OntologyEfo(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_EFO"])
    assert sclc.is_a(is_="10x 3' v3", a_="10x technology")
    assert sclc.is_a(is_="10x 3' v3", a_="10x 3' transcription profiling")
    assert not sclc.is_a(is_="10x technology", a_="10x 3' transcription profiling")
    assert sclc.is_a(is_="10x 3' v3", a_="single cell library construction")
    assert sclc.is_a(is_="sci-Plex", a_="single cell library construction")
    assert not sclc.is_a(is_="sci-Plex", a_="10x technology")


"""
Hancestro
"""


def test_hancestro_loading():
    _ = OntologyHancestro(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_HANCESTRO"], recache=True)
    _ = OntologyHancestro(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_HANCESTRO"], recache=False)


"""
Hsapdv
"""


def test_hsapdv_loading():
    _ = OntologyHsapdv(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_HSAPDV"], recache=True)
    _ = OntologyHsapdv(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_HSAPDV"], recache=False)


"""
MONDO
"""


def test_mondo_loading():
    _ = OntologyMondo(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_MONDO"], recache=True)
    _ = OntologyMondo(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_MONDO"], recache=False)


"""
Mmusdv
"""


def test_mmusdv_loading():
    _ = OntologyMmusdv(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_MMUSDV"], recache=True)
    _ = OntologyMmusdv(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_MMUSDV"], recache=False)


"""
NCBI Taxon
"""


def test_taxon_loading():
    _ = OntologyTaxon(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_NCBITAXON"], recache=True)
    _ = OntologyTaxon(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_NCBITAXON"], recache=False)


"""
Sex
"""


def test_sex_loading():
    _ = OntologyPato(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_PATO"], recache=True)
    _ = OntologyPato(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_PATO"], recache=False)


"""
UBERON
"""


def test_uberon_loading():
    _ = OntologyUberon(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_UBERON"], recache=True)
    _ = OntologyUberon(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_UBERON"], recache=False)


def test_uberon_lcs_loading():
    _ = OntologyUberonLifecyclestage(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_UBERON"], recache=False)


def test_uberon_subsetting():
    ou = OntologyUberon(branch=ONTOLOGY_VERSIONS[DEFAULT_SCHEMA]["VERSION_UBERON"])
    assert ou.is_a(is_="lobe of lung", a_="lung")
    assert ou.is_a(is_="lobe of lung", a_="lobe of lung")
    assert not ou.is_a(is_="lung", a_="lobe of lung")

    assert ou.is_a(is_="lobar bronchus", a_="lung")
    assert ou.is_a(is_="lobar bronchus", a_="lobar bronchus")
    assert not ou.is_a(is_="lung", a_="lobar bronchus")

    assert ou.is_a(is_="adipose tissue of abdominal region", a_="adipose tissue")
    assert ou.is_a(is_="adipose tissue", a_="adipose tissue")
    assert not ou.is_a(is_="adipose tissue", a_="adipose tissue of abdominal region")
