import abc
import networkx
import numpy as np
import obonet
import os
import requests
from typing import Dict, List, Tuple, Union
import warnings

from sfaira.consts.adata_fields import AdataIdsSfaira
from sfaira.versions.metadata.extensions import ONTOLOGIY_EXTENSION_HUMAN, ONTOLOGIY_EXTENSION_MOUSE

FILE_PATH = __file__

"""
Ontology managament classes.

We consider any structured collection of meta data identifiers an ontology and define classes to interact with such
data here.

- All classes inherit from Ontology()
- Onotlogies can be read as follows:
    - from string lists which are typically hardcoded in sfaira (OntologyList),
    - from .obo files which are emitted by obofoundry for example (OntologyObo)),
    - ToDo from .owl files which are emitted from EBI for example (OntologyOwl)),
    - from the EBI web API via direct queries (OntologyEbi)).

ToDo explain usage of ontology extension.
"""


class Ontology:
    leaves: List[str]

    @abc.abstractmethod
    def node_names(self):
        pass

    @abc.abstractmethod
    def map_node_suggestion(self, x: str, include_synonyms: bool = True, n_suggest: int = 10):
        """
        Map free text node name to ontology node names via fuzzy string matching.

        :param x: Free text node label which is to be matched to ontology nodes.
        :param include_synonyms: Whether to search for meaches in synonyms field of node instances, too.
        :return List of proposed matches in ontology.
        """
        pass

    def is_node(self, x: str):
        return x in self.node_names

    def validate_node(self, x: str):
        if not self.is_node(x=x):
            suggestions = self.map_node_suggestion(x=x, include_synonyms=False)
            raise ValueError(f"Node label {x} not found. Did you mean any of {suggestions}?")


class OntologyList(Ontology):
    """
    Basic unordered ontology container
    """

    def __init__(
            self,
            terms: Union[List[Union[str, bool, int]]],
            **kwargs
    ):
        self.nodes = terms

    @property
    def node_names(self) -> List[str]:
        return self.nodes

    def map_node_suggestion(self, x: str, include_synonyms: bool = True, n_suggest: int = 10):
        """
        Map free text node name to ontology node names via fuzzy string matching.

        :param x: Free text node label which is to be matched to ontology nodes.
        :param include_synonyms: Whether to search for meaches in synonyms field of node instances, too.
        :param n_suggest: number of suggestions returned
        :return List of proposed matches in ontology.
        """
        from fuzzywuzzy import fuzz
        scores = np.array([
            np.max([
                fuzz.ratio(x.lower(), y.lower())
            ])
            for y in self.node_names
        ])
        # Suggest top n_suggest hits by string match:
        return [self.node_names[i] for i in np.argsort(scores)[-n_suggest:]][::-1]

    def synonym_node_properties(self) -> List[str]:
        return []

    def is_a(self, query: str, reference: str) -> bool:
        """
        Checks if query node is reference node.

        Note that there is no notion of ancestors for list ontologies.

        :param query: Query node name. Node ID or name.
        :param reference: Reference node name. Node ID or name.
        :return: If query node is reference node or an ancestor thereof.
        """
        return query == reference


class OntologyEbi(Ontology):
    """
    Recursively assembles ontology by querying EBI web interface.

    Not recommended for large ontologies.
    Yields unstructured list of terms.
    """

    def __init__(
            self,
            ontology: str,
            root_term: str,
            additional_terms: Union[Dict[str, Dict[str, str]], None] = None,
            **kwargs
    ):
        """

        :param ontology:
        :param root_term:
        :param additional_terms: Dictionary with additional terms, values should be

            - "name" necessary
            - "description" optional
            - "synonyms" optional
            - "has_children" optional
        :param kwargs:
        """
        def get_url(iri):
            return f"https://www.ebi.ac.uk/ols/api/ontologies/{ontology}/terms/" \
                   f"http%253A%252F%252Fwww.ebi.ac.uk%252F{ontology}%252F{iri}/children"

        def recursive_search(iri):
            print(requests.get(get_url(iri=iri)).json().keys())
            terms = requests.get(get_url(iri=iri)).json()["_embedded"]["terms"]
            nodes_new = {}
            for x in terms:
                nodes_new[x["iri"].split("/")[-1]] = {
                    "name": x["label"],
                    "description": x["description"],
                    "synonyms": x["synonyms"],
                    "has_children": x["has_children"],
                }
                if x["has_children"]:
                    nodes_new.update(recursive_search(iri=x["iri"].split("/")[-1]))
            return nodes_new

        self.nodes = recursive_search(iri=root_term)
        self.nodes.update(additional_terms)

    @property
    def node_names(self) -> List[str]:
        return [v["name"] for k, v in self.nodes.items()]

    def map_node_suggestion(self, x: str, include_synonyms: bool = True, n_suggest: int = 10):
        """
        Map free text node name to ontology node names via fuzzy string matching.

        :param x: Free text node label which is to be matched to ontology nodes.
        :param include_synonyms: Whether to search for meaches in synonyms field of node instances, too.
        :return List of proposed matches in ontology.
        """
        from fuzzywuzzy import fuzz
        scores = np.array([
            np.max(
                [
                    fuzz.partial_ratio(x.lower(), v["name"].lower())
                ] + [
                    fuzz.partial_ratio(x.lower(), yyy.lower())
                    for yy in self.synonym_node_properties if yy in v.keys() for yyy in v[yy]
                ]
            ) if include_synonyms else
            np.max([
                fuzz.partial_ratio(x.lower(), v["name"].lower())
            ])
            for k, v in self.nodes.items()
        ])
        # Suggest top n_suggest hits by string match:
        return [self.node_names[i] for i in np.argsort(scores)[-n_suggest:]][::-1]

    def synonym_node_properties(self) -> List[str]:
        return ["synonyms"]

# class OntologyOwl(Ontology):
#
#    onto: owlready2.Ontology
#
#    def __init__(
#            self,
#            owl: str,
#            **kwargs
#    ):
#        self.onto = owlready2.get_ontology(owl)
#        self.onto.load()
#        # ToDo build support here
#
#    @property
#    def node_names(self):
#        pass


class OntologyObo(Ontology):

    graph: networkx.MultiDiGraph
    leaves: List[str]

    def __init__(
            self,
            obo: str,
            **kwargs
    ):
        self.graph = obonet.read_obo(obo)

    def _check_graph(self):
        if not networkx.is_directed_acyclic_graph(self.graph):
            warnings.warn("DAG was broken")

    @property
    def nodes(self) -> List[Tuple[str, dict]]:
        return list(self.graph.nodes.items())

    @property
    def nodes_dict(self) -> dict:
        return self.graph.nodes.items()

    @property
    def node_names(self) -> List[str]:
        return [x["name"] for x in self.graph.nodes.values()]

    @property
    def node_ids(self) -> List[str]:
        return list(self.graph.nodes())

    def id_from_name(self, x: str) -> str:
        self.validate_node(x=x)
        return [k for k, v in self.graph.nodes.items() if v["name"] == x][0]

    def set_leaves(self, nodes: list = None):
        # ToDo check that these are not include parents of each other!
        if nodes is not None:
            for x in nodes:
                assert x in self.graph.nodes, f"{x} not found"
            self.leaves = nodes
        else:
            self.leaves = self.get_all_roots()

    def get_all_roots(self) -> List[str]:
        return [x for x in self.graph.nodes() if self.graph.in_degree(x) == 0]

    def get_ancestors(self, node: str) -> List[str]:
        if node not in self.node_ids:
            node = self.id_from_name(node)
        return list(networkx.ancestors(self.graph, node))

    def is_a(self, query: str, reference: str) -> bool:
        """
        Checks if query node is reference node or an ancestor thereof.

        :param query: Query node name. Node ID or name.
        :param reference: Reference node name. Node ID or name.
        :return: If query node is reference node or an ancestor thereof.
        """
        if query not in self.node_ids:
            query = self.id_from_name(query)
        if reference not in self.node_ids:
            reference = self.id_from_name(reference)
        return query in self.get_ancestors(node=reference) or query == reference

    def map_to_leaves(self, node: str, return_type: str = "elements", include_self: bool = True):
        """
        Map a given list of nodes to leave nodes.

        :param node:
        :param return_type:

            "elements": names of mapped leave nodes
            "idx": indicies in leave note list of of mapped leave nodes
        :param include_self: whether to include node itself
        :return:
        """
        assert self.leaves is not None
        ancestors = self.get_ancestors(node)
        if include_self:
            ancestors = ancestors + [node]
        if return_type == "elements":
            return [x for x in self.leaves if x in ancestors]
        if return_type == "idx":
            return np.array([i for i, (x, y) in enumerate(self.leaves) if x in ancestors])

    @abc.abstractmethod
    def synonym_node_properties(self) -> List[str]:
        pass

    def map_node_suggestion(self, x: str, include_synonyms: bool = True, n_suggest: int = 10):
        """
        Map free text node name to ontology node names via fuzzy string matching.

        :param x: Free text node label which is to be matched to ontology nodes.
        :param include_synonyms: Whether to search for meaches in synonyms field of node instances, too.
        :return List of proposed matches in ontology.
        """
        from fuzzywuzzy import fuzz
        scores = np.array([
            np.max(
                [
                    fuzz.ratio(x.lower().strip("'").strip("\""), y[1]["name"].lower())
                ] + [
                    fuzz.ratio(x.lower().strip("'").strip("\"").strip("]").strip("["), yyy.lower())
                    for yy in self.synonym_node_properties if yy in y[1].keys() for yyy in y[1][yy]
                ]
            ) if "synonym" in y[1].keys() and include_synonyms else
            np.max([
                fuzz.ratio(x.lower().strip("'").strip("\""), y[1]["name"].lower())
            ])
            for y in self.nodes
        ])
        # Suggest top n_suggest hits by string match:
        return [self.nodes[i][1]["name"] for i in np.argsort(scores)[-n_suggest:]][::-1]


class OntologyExtendedObo(OntologyObo):
    """
    Basic .obo ontology extended by additional nodes and edges without breaking DAG.
    """

    def __init__(self, obo, **kwargs):
        super().__init__(obo=obo, **kwargs)
        # ToDo distinguish here:
        self.add_extension(dict_ontology=ONTOLOGIY_EXTENSION_HUMAN)

    def add_extension(self, dict_ontology: Dict[str, List[str]]):
        """
        Extend ontology by additional edges and nodes defined in a dictionary.

        Checks that DAG is not broken after graph assembly.

        :param dict_ontology: Dictionary of nodes and edges to add to ontology. Parsing:

            - keys: parent nodes (which must be in ontology)
            - values: children nodes (which can be in ontology), must be given as list of stringd.
                If these are in the ontology, an edge is added, otherwise, an edge and the node are added.
        :return:
        """
        for k, v in dict_ontology.items():
            assert isinstance(v, list), "dictionary values should be list of strings"
            # Check that parent node is present:
            if k not in self.nodes:
                raise ValueError(f"key {k} was not in reference ontology")
            # Check if edge is added only, or edge and node.
            for child_node in v:
                if child_node not in self.nodes:  # Add node.
                    self.graph.add_node(child_node)
                # Add edge.
                self.graph.add_edge(k, child_node)
        # Check that DAG was not broken:
        self._check_graph()

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]


class OntologyUberon(OntologyExtendedObo):

    def __init__(
            self,
            **kwargs
    ):
        super().__init__(obo="http://purl.obolibrary.org/obo/uberon.obo")

        # Clean up nodes:
        nodes_to_delete = []
        for k, v in self.graph.nodes.items():
            # Only retain nodes which are  "anatomical collection" 'UBERON:0034925':
            # ToDo this seems to narrow, need to check if we need to constrain the nodes we use.
            if "name" not in v.keys():
                nodes_to_delete.append(k)
        for k in nodes_to_delete:
            self.graph.remove_node(k)

        # Clean up edges:
        # The graph object can hold different types of edges,
        # and multiple types are loaded from the obo, not all of which are relevant for us:
        # All edge types (based on previous download, assert below that this is not extended):
        edge_types = [
            'aboral_to',
            'adjacent_to',
            'anastomoses_with',
            'anterior_to',
            'anteriorly_connected_to',
            'attaches_to',
            'attaches_to_part_of',
            'bounding_layer_of',
            'branching_part_of',
            'channel_for',
            'channels_from',
            'channels_into',
            'composed_primarily_of',
            'conduit_for',
            'connected_to',
            'connects',
            'contains',
            'continuous_with',
            'contributes_to_morphology_of',
            'deep_to',
            'developmentally_induced_by',
            'developmentally_replaces',
            'develops_from',  # developmental DAG -> include because it reflects the developmental hierarchy
            'develops_from_part_of',  # developmental DAG -> include because it reflects the developmental hierarchy
            'develops_in',
            'directly_develops_from',  # developmental DAG -> include because it reflects the developmental hierarchy
            'distal_to',
            'distally_connected_to',
            'distalmost_part_of',
            'dorsal_to',
            'drains',
            'ends',
            'ends_with',
            'existence_ends_during',
            'existence_ends_during_or_before',
            'existence_ends_with',
            'existence_starts_and_ends_during',
            'existence_starts_during',
            'existence_starts_during_or_after',
            'existence_starts_with',
            'extends_fibers_into',
            'filtered_through',
            'has_boundary',
            'has_component',
            'has_developmental_contribution_from',
            'has_fused_element',
            'has_member',
            'has_muscle_antagonist',
            'has_muscle_insertion',
            'has_muscle_origin',
            'has_part',
            'has_potential_to_develop_into',
            'has_potential_to_developmentally_contribute_to',
            'has_skeleton',
            'immediate_transformation_of',
            'immediately_anterior_to',
            'immediately_deep_to',
            'immediately_posterior_to',
            'immediately_preceded_by',
            'immediately_superficial_to',
            'in_anterior_side_of',
            'in_central_side_of',
            'in_deep_part_of',
            'in_distal_side_of',
            'in_dorsal_side_of',
            'in_innermost_side_of',
            'in_lateral_side_of',
            'in_left_side_of',
            'in_outermost_side_of',
            'in_posterior_side_of',
            'in_proximal_side_of',
            'in_right_side_of',
            'in_superficial_part_of',
            'in_ventral_side_of',
            'indirectly_supplies',
            'innervated_by',
            'innervates',
            'intersects_midsagittal_plane_of',
            'is_a',
            'layer_part_of',
            'located_in',  # anatomic DAG -> include because it reflects the anatomic coarseness / hierarchy
            'location_of',
            'lumen_of',
            'luminal_space_of',
            'overlaps',
            'part_of',  # anatomic DAG -> include because it reflects the anatomic coarseness / hierarchy
            'postaxialmost_part_of',
            'posterior_to',
            'posteriorly_connected_to',
            'preaxialmost_part_of',
            'preceded_by',
            'precedes',
            'produced_by',
            'produces',
            'protects',
            'proximal_to',
            'proximally_connected_to',
            'proximalmost_part_of',
            'seeAlso',
            'serially_homologous_to',
            'sexually_homologous_to',
            'skeleton_of',
            'starts',
            'starts_with',
            'subdivision_of',
            'superficial_to',
            'supplies',
            'surrounded_by',
            'surrounds',
            'transformation_of',
            'tributary_of',
            'trunk_part_of',
            'ventral_to'
        ]
        edges_to_delete = []
        for i, x in enumerate(self.graph.edges):
            assert x[2] in edge_types, x
            if x[2] not in [
                "develops_from",
                "located_in",
                "part_of",
            ]:
                edges_to_delete.append((x[0], x[1]))
        for x in edges_to_delete:
            self.graph.remove_edge(u=x[0], v=x[1])
        self._check_graph()

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym", "latin term", "has relational adjective"]


class OntologyCelltypes(OntologyExtendedObo):

    def __init__(
            self,
            branch: str,
            **kwargs
    ):
        if os.name == "nt":  # if running on windows, do not download obo file, but rather pass url directly to obonet
            obofile = f"https://raw.github.com/obophenotype/cell-ontology/{branch}/cl.obo"
        else:
            # Identify cache:
            folder = FILE_PATH.split(os.sep)[:-4]
            folder.insert(1, os.sep)
            ontology_cache_dir = os.path.join(*folder, "cache", "ontologies", "cl")
            fn = f"{branch}_cl.obo"
            obofile = os.path.join(ontology_cache_dir, fn)
            # Download if necessary:
            if not os.path.isfile(obofile):
                def download_cl():
                    url = f"https://raw.github.com/obophenotype/cell-ontology/{branch}/cl.obo"
                    print(f"Downloading: {fn}")
                    if not os.path.exists(ontology_cache_dir):
                        os.makedirs(ontology_cache_dir)
                    r = requests.get(url, allow_redirects=True)
                    open(obofile, 'wb').write(r.content)
                download_cl()

        super().__init__(obo=obofile)

        # Clean up nodes:
        nodes_to_delete = []
        for k, v in self.graph.nodes.items():
            # Some terms are not associated with the namespace cell but are cell types,
            # we identify these based on their ID nomenclature here.
            if ("namespace" in v.keys() and v["namespace"] not in ["cell", "cl"]) or \
                    ("namespace" not in v.keys() and str(k)[:2] != "CL"):
                nodes_to_delete.append(k)
            elif "name" not in v.keys():
                nodes_to_delete.append(k)
        for k in nodes_to_delete:
            self.graph.remove_node(k)

        # Clean up edges:
        # The graph object can hold different types of edges,
        # and multiple types are loaded from the obo, not all of which are relevant for us:
        # All edge types (based on previous download, assert below that this is not extended):
        edge_types = [
            'is_a',  # nomenclature DAG -> include because of annotation coarseness differences
            'derives_from',
            'develops_from',  # developmental DAG -> include because of developmental differences
            'has_part',  # ?
            'develops_into',  # inverse developmental DAG -> do not include
            'part_of',
            'RO:0002120',  # ?
            'RO:0002103',  # ?
            'lacks_plasma_membrane_part',  # ?
        ]
        edges_to_delete = []
        for i, x in enumerate(self.graph.edges):
            assert x[2] in edge_types, x
            if x[2] not in ["is_a", "develops_from"]:
                edges_to_delete.append((x[0], x[1]))
        for x in edges_to_delete:
            self.graph.remove_edge(u=x[0], v=x[1])
        self._check_graph()

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]


# use OWL for OntologyHancestro


class OntologyHsapdv(OntologyExtendedObo):

    def __init__(
            self,
            **kwargs
    ):
        super().__init__(obo="http://purl.obolibrary.org/obo/hsapdv.obo")

        # Clean up nodes:
        nodes_to_delete = []
        for k, v in self.graph.nodes.items():
            if "name" not in v.keys():
                nodes_to_delete.append(k)
        for k in nodes_to_delete:
            self.graph.remove_node(k)

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]


class OntologyMmusdv(OntologyExtendedObo):

    def __init__(
            self,
            **kwargs
    ):
        super().__init__(obo="http://purl.obolibrary.org/obo/mmusdv.obo")

        # Clean up nodes:
        nodes_to_delete = []
        for k, v in self.graph.nodes.items():
            if "name" not in v.keys():
                nodes_to_delete.append(k)
        for k in nodes_to_delete:
            self.graph.remove_node(k)

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]


class OntologyMondo(OntologyObo):

    def __init__(
            self,
            **kwargs
    ):
        super().__init__(obo="http://purl.obolibrary.org/obo/mondo.obo")

        # Clean up nodes:
        nodes_to_delete = []
        for k, v in self.graph.nodes.items():
            if "name" not in v.keys():
                nodes_to_delete.append(k)
        for k in nodes_to_delete:
            self.graph.remove_node(k)

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]


class OntologyCellosaurus(OntologyExtendedObo):

    def __init__(
            self,
            **kwargs
    ):
        download_link = "https://ftp.expasy.org/databases/cellosaurus/cellosaurus.obo"

        if os.name == "nt":  # if running on windows, do not download obo file, but rather pass url directly to obonet
            super().__init__(obo=download_link)
        else:
            # Identify cache:
            folder = FILE_PATH.split(os.sep)[:-4]
            folder.insert(1, os.sep)
            ontology_cache_dir = os.path.join(*folder, "cache", "ontologies", "cellosaurus")
            fn = "cellosaurus.obo"
            obofile = os.path.join(ontology_cache_dir, fn)
            # Download if necessary:
            if not os.path.isfile(obofile):
                def download_cl():
                    print(f"Downloading: {fn}")
                    if not os.path.exists(ontology_cache_dir):
                        os.makedirs(ontology_cache_dir)
                    r = requests.get(download_link, allow_redirects=True)
                    open(obofile, 'wb').write(r.content)
                download_cl()
            super().__init__(obo=obofile)

        # Clean up nodes:
        # edge_types = ["derived_from", "originate_from_same_individual_as"]
        nodes_to_delete = []
        for k, v in self.graph.nodes.items():
            if "name" not in v.keys():
                nodes_to_delete.append(k)
        for k in nodes_to_delete:
            self.graph.remove_node(k)

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]


class OntologySinglecellLibraryConstruction(OntologyEbi):

    def __init__(
            self,
            ontology: str = "efo",
            root_term: str = "EFO_0010183",
    ):
        super().__init__(
            ontology=ontology,
            root_term=root_term,
            additional_terms={
                "microwell-seq": {"name": "microwell-seq"}
            }
        )
