import abc
import networkx
import numpy as np
import obonet
import os
# import owlready2
import pickle

import pandas as pd
import requests
from functools import lru_cache
from typing import Dict, List, Tuple, Union


from sfaira import settings
from sfaira.versions.metadata.extensions import EFO_TERMS

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


def cached_load_file(url, ontology_cache_dir, ontology_cache_fn, recache: bool = False):
    if os.name == "nt":  # if running on windows, do not download obo file, but rather pass url directly to obonet
        # TODO add caching option.
        obofile = url
    else:
        ontology_cache_dir = os.path.join(settings.cachedir_ontologies, ontology_cache_dir)
        obofile = os.path.join(ontology_cache_dir, ontology_cache_fn)
        # Download if necessary:
        if not os.path.exists(obofile) or recache:
            os.makedirs(name=ontology_cache_dir, exist_ok=True)

            def download_file():
                print(f"Downloading: {ontology_cache_fn} to {ontology_cache_dir}")
                if not os.path.exists(ontology_cache_dir):
                    os.makedirs(ontology_cache_dir)
                r = requests.get(url, allow_redirects=True)
                # if url.startswith("https://raw.githubusercontent.com"):
                #     open(obofile, 'wb').write(r.text)
                # else:
                open(obofile, 'wb').write(r.content)

            download_file()
    return obofile


def cached_load_ebi(ontology_cache_dir, ontology_cache_fn, recache: bool = False) -> (networkx.MultiDiGraph, os.PathLike):
    """
    Load pickled graph object if available.

    :param ontology_cache_dir:
    :param ontology_cache_fn:
    :param recache:
    :return:
    """
    ontology_cache_dir = os.path.join(settings.cachedir_ontologies, ontology_cache_dir)
    picklefile = os.path.join(ontology_cache_dir, ontology_cache_fn)
    if os.path.isfile(picklefile) and not recache:
        with open(picklefile, 'rb') as f:
            graph = pickle.load(f)
    else:
        os.makedirs(name=ontology_cache_dir, exist_ok=True)
        graph = None
    return graph, picklefile


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
    Basic unordered ontology container.

    Node IDs and names are the same.
    """
    nodes: Union[dict, list]

    def __init__(
            self,
            terms: Union[List[Union[str, bool, int]], Dict[str, dict]],
            **kwargs
    ):
        self.nodes = terms
        self._dict_mode = isinstance(terms, dict)

    @property
    def node_names(self) -> List[str]:
        if self._dict_mode:
            return [v["name"] for v in self.nodes.values()]
        else:
            return self.nodes

    @property
    def node_ids(self) -> List[str]:
        if self._dict_mode:
            return list(self.nodes.keys())
        else:
            return self.nodes

    @property
    def leaves(self) -> List[str]:
        return self.node_ids

    @property
    def n_leaves(self) -> int:
        return len(self.node_ids)

    def map_to_leaves(
            self,
            node: str,
            return_type: str = "ids",
    ) -> Union[List[str], np.ndarray]:
        """
        Map a given node to leave nodes.

        Note that for a list ontology, this maps each node to itself.

        :param node: Node(s) to map as symbol(s) or ID(s).
        :param return_type:

            "ids": IDs of mapped leave nodes. This is the query node itself for OntologyList.
            "idx": indices in leave note list of mapped leave nodes. This is the index of query node itself for
                    OntologyList.
        :return:
        """
        node = self.__convert_to_id_cached(node)
        if return_type == "ids":
            return [node]
        elif return_type == "idx":
            return np.where(self.node_ids == node)[0]
        else:
            raise ValueError(f"return_type {return_type} not recognized")

    def prepare_maps_to_leaves(
            self,
            include_self: bool = True
    ) -> Dict[str, np.ndarray]:
        """
        Precomputes all maps of nodes to their leave nodes.

        Note that for a list ontology, this maps each node to itself.

        :param include_self: whether to include node itself
        :return: Dictionary of index vectors of leave node matches for each node (key).
        """
        if include_self:
            return dict([(x, np.array([self.leaves.index(x)])) for x in self.leaves])
        else:
            return dict([(x, np.array([])) for x in self.leaves])

    def is_a_node_id(self, x: str) -> bool:
        return x in self.node_ids

    def is_a_node_name(self, x: str) -> bool:
        return x in self.node_names

    def convert_to_name(self, x: Union[str, List[str]]) -> Union[str, List[str]]:
        if self._dict_mode:
            was_str = isinstance(x, str)
            if was_str:
                x = [x]
            if self.is_a_node_id(x[0]):
                self.__validate_node_ids(x=x)
                x = [self.nodes[k]["name"] for k in x]
            elif self.is_a_node_name(x[0]):
                self.__validate_node_names(x=x)
            else:
                raise ValueError(f"node {x[0]} not recognized")
            if was_str:
                x = x[0]
        return x

    def convert_to_id(self, x: Union[str, List[str]]) -> Union[str, List[str]]:
        if self._dict_mode:
            was_str = isinstance(x, str)
            if was_str:
                x = [x]
            if self.is_a_node_id(x[0]):
                self.__validate_node_ids(x=x)
            elif self.is_a_node_name(x[0]):
                self.__validate_node_names(x=x)
                x = [
                    [k for k, v in self.nodes.items() if v["name"] == z][0]
                    for z in x
                ]
            else:
                raise ValueError(f"node {x[0]} not recognized")
            if was_str:
                x = x[0]
        return x

    @lru_cache(maxsize=None)
    def __convert_to_id_cached(self, x: str) -> str:
        return self.convert_to_id(x)

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

    def get_ancestors(self, node: str) -> List[str]:
        return []

    def get_descendants(self, node: str) -> List[str]:
        return []

    def __validate_node_ids(self, x: Union[str, List[str]]):
        if isinstance(x, str):
            x = [x]
        node_ids = self.node_ids
        if np.any([y not in node_ids for y in x]):
            raise ValueError(f"queried node id {x} are not all in graph")

    def __validate_node_names(self, x: Union[str, List[str]]):
        if isinstance(x, str):
            x = [x]
        node_names = self.node_names
        if np.any([y not in node_names for y in x]):
            raise ValueError(f"queried node names {x} are not all in graph")


class OntologyHierarchical(Ontology, abc.ABC):
    """
    Basic ordered ontology container
    """
    _graph: networkx.MultiDiGraph
    _node_names: Union[None, List[str]] = None

    def _clear_caches(self):
        self.get_ancestors.cache_clear()
        self.get_descendants.cache_clear()
        self._node_names = None

    @property
    def graph(self) -> networkx.MultiDiGraph:
        return self._graph

    @graph.setter
    def graph(self, graph: networkx.MultiDiGraph):
        self._graph = graph
        self._clear_caches()

    def _check_graph(self, verbose=0):
        if not networkx.is_directed_acyclic_graph(self.graph) and verbose > 0:
            print(f"Ontology {type(self)} is not a DAG, treat child-parent reasoning with care.")

    def __validate_node_ids(self, x: Union[str, List[str]]):
        if isinstance(x, str):
            x = [x]
        node_ids = self.node_ids
        if np.any([y not in node_ids for y in x]):
            raise ValueError(f"queried node id {x} are not all in graph")

    def __validate_node_names(self, x: Union[str, List[str]]):
        if isinstance(x, str):
            x = [x]
        node_names = self.node_names
        if np.any([y not in node_names for y in x]):
            raise ValueError(f"queried node names {x} are not all in graph")

    @property
    def nodes(self) -> List[Tuple[str, dict]]:
        return list(self.graph.nodes.items())

    @nodes.setter
    def nodes(self, x: List[str]):
        """
        Sets new nodes-space for graph.

        This clips nodes that are not in target set.
        :param x: New set of nodes, identified as IDs.
        """
        x = self.convert_to_id(x)
        nodes_to_remove = []
        for y in self.node_ids:
            if y not in x:
                nodes_to_remove.append(y)
        self.graph.remove_nodes_from(nodes_to_remove)
        self._clear_caches()

    @property
    def nodes_dict(self) -> dict:
        return dict(list(self.graph.nodes.items()))

    @property
    def node_names(self) -> List[str]:
        if self._node_names is None:
            self._node_names = [v["name"] for v in self.graph.nodes.values()]
        return self._node_names

    @property
    def node_ids(self) -> List[str]:
        return list(self.graph.nodes())

    def is_a_node_id(self, x: str) -> bool:
        return x in self.node_ids

    def is_a_node_name(self, x: str) -> bool:
        return x in self.node_names

    def convert_to_name(self, x: Union[str, List[str]]) -> Union[str, List[str]]:
        was_str = isinstance(x, str)
        if was_str:
            x = [x]
        if self.is_a_node_id(x[0]):
            self.__validate_node_ids(x=x)
            x = [self.graph.nodes[k]["name"] for k in x]
        elif self.is_a_node_name(x[0]):
            self.__validate_node_names(x=x)
        else:
            raise ValueError(f"node {x[0]} not recognized")
        if was_str:
            return x[0]
        else:
            return x

    def convert_to_id(self, x: Union[str, List[str]]) -> Union[str, List[str]]:
        was_str = isinstance(x, str)
        if was_str:
            x = [x]
        if self.is_a_node_id(x[0]):
            self.__validate_node_ids(x=x)
        elif self.is_a_node_name(x[0]):
            self.__validate_node_names(x=x)
            x = [
                [k for k, v in self.graph.nodes.items() if v["name"] == z][0]
                for z in x
            ]
        else:
            raise ValueError(f"node {x[0]} not recognized")
        if was_str:
            return x[0]
        else:
            return x

    @lru_cache(maxsize=None)
    def __convert_to_id_cached(self, x: str) -> str:
        return self.convert_to_id(x)

    @property
    def leaves(self) -> List[str]:
        return [x for x in self.graph.nodes() if self.graph.in_degree(x) == 0]

    @leaves.setter
    def leaves(self, x: List[str]):
        """
        Sets new leaf-space for graph.

        This clips nodes that are not upstream of defined leaves.
        :param x: New set of leaves nodes, identified as IDs.
        """
        x = self.convert_to_id(x)
        nodes_to_remove = []
        for y in self.graph.nodes():
            if not np.any([self.is_a(query=z, reference=y, convert_to_id=False) for z in x]):
                nodes_to_remove.append(y)
        self.graph.remove_nodes_from(nodes_to_remove)
        self._clear_caches()

    @property
    def n_leaves(self) -> int:
        return len(self.leaves)

    def get_effective_leaves(self, x: List[str]) -> List[str]:
        """
        Get effective leaves in ontology given set of observed nodes.

        The effective leaves are the minimal set of nodes such that all nodes in x are ancestors of this set, ie the
        observed nodes which represent leaves of a sub-DAG of the ontology DAG, which captures all observed nodes.

        :param x: Observed node IDs.
        :return: Effective leaves.
        """
        if isinstance(x, str):
            x = [x]
        if isinstance(x, np.ndarray):
            x = x.tolist()
        assert isinstance(x, list), "supply either list or str to get_effective_leaves"
        if len(x) == 0:
            raise ValueError("x was empty list, get_effective_leaves cannot be called on empty list")
        x = list(np.unique(x))
        x = self.convert_to_id(x=x)
        leaves = []
        for y in x:
            if not np.any([self.is_a(query=z, reference=y, convert_to_id=False) for z in list(set(x) - {y})]):
                leaves.append(y)
        return leaves

    @lru_cache(maxsize=None)
    def get_ancestors(self, node: str) -> List[str]:
        node = self.__convert_to_id_cached(node)
        return list(networkx.ancestors(self.graph, node))

    @lru_cache(maxsize=None)
    def get_descendants(self, node: str) -> List[str]:
        node = self.__convert_to_id_cached(node)
        return list(networkx.descendants(self.graph, node))

    def is_a(self, query: str, reference: str, convert_to_id: bool = True) -> bool:
        """
        Checks if query node is reference node or an ancestor thereof.

        :param query: Query node name. Node ID or name.
        :param reference: Reference node name. Node ID or name.
        :param convert_to_id: Whether to call self.convert_to_id on `query` and `reference` input arguments
        :return: If query node is reference node or an ancestor thereof.
        """
        if convert_to_id:
            query = self.__convert_to_id_cached(query)
            reference = self.__convert_to_id_cached(reference)
        return query in self.get_ancestors(node=reference) or query == reference

    def map_to_leaves(
            self,
            node: str,
            return_type: str = "ids",
            include_self: bool = True
    ) -> Union[List[str], np.ndarray]:
        """
        Map a given node to leave nodes.

        :param node: Node(s) to map as symbol(s) or ID(s).
        :param return_type:

            "ids": IDs of mapped leave nodes
            "idx": indices in leave note list of mapped leave nodes
        :param include_self: DEPRECEATED.
        :return:
        """
        node = self.__convert_to_id_cached(node)
        ancestors = self.get_ancestors(node)
        # Add node itself to list of ancestors.
        ancestors = ancestors + [node]
        leaves = self.leaves
        if return_type == "ids":
            return [x for x in leaves if x in ancestors]
        elif return_type == "idx":
            return np.sort([i for i, x in enumerate(leaves) if x in ancestors])
        else:
            raise ValueError(f"return_type {return_type} not recognized")

    def prepare_maps_to_leaves(
            self,
            include_self: bool = True
    ) -> Dict[str, np.ndarray]:
        """
        Precomputes all maps of nodes to their leave nodes.

        :param include_self: whether to include node itself
        :return: Dictionary of index vectors of leave node matches for each node (key).
        """
        nodes = self.node_ids
        maps = {}
        import time
        t0 = time.time()
        for x in nodes:
            maps[x] = self.map_to_leaves(node=x, return_type="idx", include_self=include_self)
        print(f"time for precomputing ancestors: {time.time()-t0}")
        return maps

    def reset_root(self, root: Union[str, List[str]]):
        if isinstance(root, str):
            root = [root]
        new_nodes = self.convert_to_id(x=root)
        for x in root:
            new_nodes += self.get_ancestors(node=x)
        self.graph = self.graph.subgraph(nodes=new_nodes).copy()

    @abc.abstractmethod
    def synonym_node_properties(self) -> List[str]:
        pass


class OntologyEbi(OntologyHierarchical):
    """
    Recursively assembles ontology by querying EBI web interface.

    Not recommended for large ontologies because of the iterative query of the web API.
    """

    def __init__(
            self,
            ontology: str,
            root_term: str,
            additional_terms: dict,
            additional_edges: List[Tuple[str, str]],
            ontology_cache_fn: str,
            recache: bool,
            **kwargs
    ):
        # Note on base URL: EBI OLS points to different resources depending on the ontology used, this needs to be
        # accounted for here.
        if ontology == "hancestro":
            base_url = f"https://www.ebi.ac.uk/ols/api/ontologies/{ontology}/terms/" \
                       f"http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F"
        elif ontology == "efo":
            base_url = f"https://www.ebi.ac.uk/ols/api/ontologies/{ontology}/terms/" \
                       f"http%253A%252F%252Fwww.ebi.ac.uk%252F{ontology}%252F"
        else:
            assert False

        def get_url_self(iri):
            return f"{base_url}{iri}"

        def get_url_children(iri):
            return f"{base_url}{iri}/children"

        def get_iri_from_node(x):
            return x["iri"].split("/")[-1]

        def get_id_from_iri(x):
            x = ":".join(x.split("_"))
            return x

        def get_id_from_node(x):
            x = get_iri_from_node(x)
            x = get_id_from_iri(x)
            return x

        def recursive_search(iri):
            """
            This function queries all nodes that are children of a given node at one time. This is faster than querying
            the characteristics of each node separately but leads to slightly awkward code, the root node has to be
            queried separately for example below.

            :param iri: Root node IRI.
            :return: Tuple of

                - nodes (dictionaries of node ID and node values) and
                - edges (node ID of parent and child).
            """
            terms_children = requests.get(get_url_children(iri=iri)).json()["_embedded"]["terms"]
            nodes_new = {}
            edges_new = []
            direct_children = []
            k_self = get_id_from_iri(iri)
            # Define root node if this is the first iteration, this node is otherwise not defined through values.
            if k_self == ":".join(root_term.split("_")):
                terms_self = requests.get(get_url_self(iri=iri)).json()
                nodes_new[k_self] = {
                    "name": terms_self["label"],
                    "description": terms_self["description"],
                    "synonyms": terms_self["synonyms"],
                    "has_children": terms_self["has_children"],
                }
            for c in terms_children:
                k_c = get_id_from_node(c)
                nodes_new[k_c] = {
                    "name": c["label"],
                    "description": c["description"],
                    "synonyms": c["synonyms"],
                    "has_children": c["has_children"],
                }
                direct_children.append(k_c)
                if c["has_children"]:
                    nodes_x, edges_x = recursive_search(iri=get_iri_from_node(c))
                    nodes_new.update(nodes_x)
                    # Update nested edges of between children:
                    edges_new.extend(edges_x)
            # Update edges to children:
            edges_new.extend([(k_self, k_c) for k_c in direct_children])
            return nodes_new, edges_new

        graph, picklefile = cached_load_ebi(ontology_cache_dir=ontology, ontology_cache_fn=ontology_cache_fn,
                                            recache=recache)
        if graph is None:
            self.graph = networkx.MultiDiGraph()
            nodes, edges = recursive_search(iri=root_term)
            nodes.update(additional_terms)
            edges.extend(additional_edges)
            for k, v in nodes.items():
                self.graph.add_node(node_for_adding=k, **v)
            for x in edges:
                parent, child = x
                self.graph.add_edge(child, parent)
            with open(picklefile, 'wb') as f:
                pickle.dump(obj=self.graph, file=f)
        else:
            self.graph = graph

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
            for k, v in self.graph.nodes.items()
        ])
        # Suggest top n_suggest hits by string match:
        return [self.node_names[i] for i in np.argsort(scores)[-n_suggest:]][::-1]

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonyms"]


class OntologyOwl(OntologyHierarchical, abc.ABC):

    # onto_owl = owlready2.Ontology

    def __init__(
            self,
            owl: str,
            **kwargs
    ):
        # self.onto_owl = owlready2.get_ontology(owl)
        # self.onto_owl.load()
        self.graph = None


class OntologyObo(OntologyHierarchical, abc.ABC):

    def __init__(
            self,
            obo: str,
            **kwargs
    ):
        self.graph = obonet.read_obo(obo, ignore_obsolete=True)

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


class OntologyExtended:
    """
    Basic ontology extended by additional nodes and edges without breaking DAG.
    """

    def add_children(self, dict_ontology: Dict[str, Dict[str, dict]]):
        """
        Extend ontology by additional edges and children nodes defined in a dictionary.

        Checks that DAG is not broken after graph assembly.

        :param dict_ontology: Dictionary of nodes and edges to add to ontology. Parsing:

            - keys: parent nodes (which must be in ontology)
            - values: children nodes (which can be in ontology), must be given as a dictionary in which keys are
                ontology IDs and values are node values.
                If these are in the ontology, an edge is added, otherwise, an edge and the node are added.
        :return:
        """
        for k, v in dict_ontology.items():
            assert isinstance(v, dict), "dictionary values should be dictionaries"
            # Check that parent node is present:
            if k not in self.node_ids:
                raise ValueError(f"key {k} was not in reference ontology")
            # Check if edge is added only, or edge and node.
            for child_node_k, child_node_v in v.items():
                if child_node_k not in self.node_ids:  # Add node
                    self.graph.add_node(node_for_adding=child_node_k, **child_node_v)
                # Add edge.
                self.graph.add_edge(child_node_k, k)
        # Check that DAG was not broken:
        self._check_graph()

    def add_parents(self, dict_ontology: Dict[str, List[str]]):
        """
        Extend ontology by additional edges and parent nodes defined in a dictionary.

        Checks that DAG is not broken after graph assembly.

        :param dict_ontology: Dictionary of nodes and edges to add to ontology. Parsing:

            - keys: parent nodes (which are not in ontology)
            - values: children nodes (which must be in ontology)
        :return:
        """
        for k, v in dict_ontology.items():
            # Add parent:
            self.graph.add_node(node_for_adding=k, **{"name": k})
            # Check if edge is added only, or edge and node.
            v = self.convert_to_id(x=v)
            for child_node_k in v:
                if child_node_k not in self.node_ids:  # Add node
                    raise ValueError(f"key {child_node_k} was not in reference ontology")
                # Add edge.
                self.graph.add_edge(child_node_k, k)
        # Check that DAG was not broken:
        self._check_graph(verbose=0)

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]


class OntologyUberon(OntologyObo, OntologyExtended):

    def __init__(
            self,
            branch: str,
            recache: bool = False,
            **kwargs
    ):
        obofile = cached_load_file(
            url=f"https://raw.githubusercontent.com/obophenotype/uberon/{branch}/uberon.obo",
            ontology_cache_dir="uberon",
            ontology_cache_fn=f"uberon_{branch}.obo",
            recache=recache,
        )
        super().__init__(obo=obofile)

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
        edges_to_delete = []
        for i, x in enumerate(self.graph.edges):
            if x[2] not in [
                "is_a",
                "located_in",
                "part_of",
            ]:
                edges_to_delete.append((x[0], x[1]))
        for x in edges_to_delete:
            self.graph.remove_edge(u=x[0], v=x[1])
        self._check_graph(verbose=0)

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym", "latin term", "has relational adjective"]


class OntologyUberonLifecyclestage(OntologyUberon):

    """
    Subset of UBERON for generic life cycle stages that can be used for organism not covered by specific developmental
    ontologies.
    """

    def __init__(
            self,
            branch: str,
            recache: bool = False,
            **kwargs
    ):
        super().__init__(branch=branch, recache=recache, **kwargs)
        self.reset_root(root="UBERON:0000105")


class OntologyCl(OntologyObo, OntologyExtended):

    def __init__(
            self,
            branch: str,
            use_developmental_relationships: bool = False,
            recache: bool = False,
            **kwargs
    ):
        """

        Developmental edges are not desired in all interactions with this ontology, double-negative thymocytes are for
        example not an intuitive parent node for a fine grained T cell label in a non-thymic tissue.
        :param branch:
        :param use_developmental_relationships: Whether to keep developmental relationships.
        :param kwargs:
        """
        obofile = cached_load_file(
            url=f"https://raw.github.com/obophenotype/cell-ontology/{branch}/cl.obo",
            ontology_cache_dir="cl",
            ontology_cache_fn=f"cl_{branch}.obo",
            recache=recache,
        )
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
        if use_developmental_relationships:
            edges_allowed = ["is_a", "develops_from"]
        else:
            edges_allowed = ["is_a"]
        for i, x in enumerate(self.graph.edges):
            if x[2] not in edge_types:
                print(f"NON-CRITICAL WARNING: cl edge type {x[2]} not in reference list yet")
            if x[2] not in edges_allowed:
                edges_to_delete.append((x[0], x[1]))
        for x in edges_to_delete:
            self.graph.remove_edge(u=x[0], v=x[1])
        self._check_graph(verbose=0)

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]


class OntologyOboCustom(OntologyObo, OntologyExtended):

    def __init__(
            self,
            obo: str,
            **kwargs
    ):
        super().__init__(obo=obo, **kwargs)


# use OWL for OntologyHancestro


class OntologyHsapdv(OntologyObo):

    def __init__(
            self,
            branch: str,
            recache: bool = False,
            **kwargs
    ):
        # Note on URL: berkeleybop is the ontology supported by cellxgene, GitHub seems to be older but versioned.
        obofile = cached_load_file(
            url="http://ontologies.berkeleybop.org/hsapdv.obo",
            # url=f"https://raw.githubusercontent.com/obophenotype/developmental-stage-ontologies/{branch}/src/hsapdv/hsapdv.obo",
            ontology_cache_dir="hsapdv",
            ontology_cache_fn=f"hsapdv_{branch}.obo",
            recache=recache,
        )
        super().__init__(obo=obofile)

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


class OntologyMmusdv(OntologyObo, OntologyExtended):

    def __init__(
            self,
            branch: str,
            recache: bool = False,
            **kwargs
    ):
        # URL for releases, not used here yet because versioning with respect to releases below is not consistent yet.
        # url=f"https://raw.githubusercontent.com/obophenotype/developmental-stage-ontologies/{branch}/src/mmusdv/mmusdv.obo"
        obofile = cached_load_file(
            url="http://ontologies.berkeleybop.org/mmusdv.obo",
            ontology_cache_dir="mmusdv",
            ontology_cache_fn="mmusdv.obo",
            recache=recache,
        )
        super().__init__(obo=obofile)

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


class OntologyMondo(OntologyObo, OntologyExtended):

    def __init__(
            self,
            branch: str,
            recache: bool = False,
            **kwargs
    ):
        # TODO Need to switch to owl, the newer releases of MONDO only push the owl release to GitHub, so versioned
        #  ontology access is only possible with owl.
        obofile = cached_load_file(
            url=f"http://purl.obolibrary.org/obo/mondo.obo",
            ontology_cache_dir="mondo",
            ontology_cache_fn=f"mondo_{branch}.obo",
            recache=recache,
        )
        super().__init__(obo=obofile)

        # Clean up nodes:
        nodes_to_delete = []
        for k, v in self.graph.nodes.items():
            if "name" not in v.keys():
                nodes_to_delete.append(k)
        for k in nodes_to_delete:
            self.graph.remove_node(k)

        # add healthy property
        # Add node "healthy" under root node "MONDO:0000001": "quality".
        # We use a PATO node for this label: PATO:0000461.
        self.add_children(dict_ontology={
            "MONDO:0000001": {
                "PATO:0000461": {"name": "healthy"}
            },
        })

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]


class OntologyCellosaurus(OntologyObo, OntologyExtended):

    def __init__(
            self,
            recache: bool = False,
            **kwargs
    ):
        obofile = cached_load_file(
            url="https://ftp.expasy.org/databases/cellosaurus/cellosaurus.obo",
            ontology_cache_dir="cellosaurus",
            ontology_cache_fn="cellosaurus.obo",
            recache=recache,
        )
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


class OntologyHancestro(OntologyEbi, OntologyExtended):

    """
    TODO move this to .owl backend once available.
    TODO root term: No term HANCESTRO_0001 ("Thing"?) accessible through EBI interface, because of that country-related
        higher order terms are not available as they are parallel to HANCESTRO_0004. Maybe fix with .owl backend?
    """

    def __init__(self, recache: bool = False):
        super().__init__(
            ontology="hancestro",
            root_term="HANCESTRO_0004",
            additional_terms={},
            additional_edges=[],
            ontology_cache_fn="hancestro.pickle",
            recache=recache,
        )

        # Add additional nodes under root node.
        self.add_children(dict_ontology={
            "HANCESTRO:0004": {"multiethnic": {"name": "multiethnic"}},
        })


class OntologyTaxon(OntologyObo, OntologyExtended):

    """
    Note on ontology: The same repo also contains ncbitaxon.obo, the full ontology which is ~500MB large and
    takes multiple minutes to load. We are using a reduced version, taxslim, here.

    See also https://github.com/obophenotype/ncbitaxon/releases/download/{branch}/ncbitaxon.obo.
    """

    def __init__(
            self,
            branch: str,
            recache: bool = False,
            **kwargs
    ):
        obofile = cached_load_file(
            url=f"https://github.com/obophenotype/ncbitaxon/releases/download/{branch}/taxslim.obo",
            ontology_cache_dir="ncbitaxon",
            ontology_cache_fn=f"ncbitaxon_{branch}.obo",
            recache=recache,
        )
        super().__init__(obo=obofile)

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


class OntologyEfo(OntologyObo, OntologyExtended):

    def __init__(
            self,
            recache: bool = False,
            **kwargs
    ):
        obofile = cached_load_file(
            url="https://www.ebi.ac.uk/efo/efo.obo",
            ontology_cache_dir="efo",
            ontology_cache_fn="efo.obo",
            recache=recache,
        )
        super().__init__(obo=obofile)
        # Subset EFO to relevant sub trees:
        self.reset_root(root=["OBI:0001686", "EFO:0010183", "EFO:0002772"])

        # Clean up nodes:
        nodes_to_delete = []
        for k, v in self.graph.nodes.items():
            if "name" not in v.keys():
                nodes_to_delete.append(k)
        for k in nodes_to_delete:
            self.graph.remove_node(k)

        self.add_parents(dict_ontology=EFO_TERMS)

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]


class OntologySex(OntologyObo, OntologyExtended):

    """
    Sex is defined based on a subset of the PATO ontology.
    """

    def __init__(
            self,
            branch: str,
            recache: bool = False,
            **kwargs
    ):
        obofile = cached_load_file(
            url=f"https://raw.githubusercontent.com/pato-ontology/pato/{branch}/pato-base.obo",
            ontology_cache_dir="pato",
            ontology_cache_fn=f"pato_{branch}.obo",
            recache=recache,
        )
        super().__init__(obo=obofile)
        nodes_to_delete = []
        # Subset ontology: see also:
        # https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/2.0.0/schema.md#sex_ontology_term_id
        targets = self.get_ancestors("PATO:0001894")
        for k, v in self.graph.nodes.items():
            if k not in targets:
                nodes_to_delete.append(k)
        # Clean up nodes:
        for k, v in self.graph.nodes.items():
            if "name" not in v.keys():
                nodes_to_delete.append(k)
        nodes_to_delete = np.unique(nodes_to_delete)
        for k in nodes_to_delete:
            self.graph.remove_node(k)

    @property
    def synonym_node_properties(self) -> List[str]:
        return ["synonym"]
