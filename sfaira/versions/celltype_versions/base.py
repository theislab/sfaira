import abc
import networkx
import numpy as np
import obonet
import pandas as pd
from typing import Dict, List, Tuple, Union
import warnings

from sfaira.versions.celltype_versions.extensions import ONTOLOGIY_EXTENSION_HUMAN, ONTOLOGIY_EXTENSION_MOUSE


class OntologyBase:
    leaves: list

    @abc.abstractmethod
    def set_leaves(self, nodes: list = None):
        pass

    @abc.abstractmethod
    def get_ancestors(self, node: str) -> List[str]:
        pass

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


class OntologyDict(OntologyBase):

    def __init__(self, onto: dict):
        self.onto = onto

    def set_leaves(self, nodes: list = None):
        self.leaves = nodes

    def get_ancestors(self, node: str) -> List[str]:
        return self.onto[node] if node in self.onto.keys() else [node]


class OntologyObo(OntologyBase):

    graph: networkx.MultiDiGraph

    def __init__(
            self,
            obo: str = "http://purl.obolibrary.org/obo/cl.obo",
            namespace_id: str = "cell",
            **kwargs
    ):
        self.graph = obonet.read_obo(obo)

        # Clean up nodes:
        nodes_to_delete = []
        for k, v in self.graph.nodes.items():
            if "namespace" not in v.keys() or v["namespace"] != namespace_id:
                nodes_to_delete.append(k)
            elif "name" not in v.keys():
                nodes_to_delete.append(k)
        for k in nodes_to_delete:
            self.graph.remove_node(k)

        # Clean up edges:
        # The graph object can hold different types of edges, assert that they are "is_a" edges:
        # All edge types (based on previous download, assert below that this is not extended):
        edge_types = [
            'is_a',  # nomenclature DAG -> include because of annotation coarseness differences
            'develops_from',  # developmental DAG -> include because of developmental differences
            'has_part',  # ?
            'develops_into',  # inverse developmental DAG -> do not include
            'RO:0002120',  # ?
            'RO:0002103',  # ?
            'lacks_plasma_membrane_part'  # ?
        ]
        edges_to_delete = []
        for i, x in enumerate(self.graph.edges):
            assert x[2] in edge_types, x
            if x[2] not in ["is_a", "develops_from"]:
                edges_to_delete.append((x[0], x[1]))
        for x in edges_to_delete:
            self.graph.remove_edge(u=x[0], v=x[1])
        self._check_graph()

    def _check_graph(self):
        if not networkx.is_directed_acyclic_graph(self.graph):
            warnings.warn("DAG was broken")

    @property
    def nodes(self):
        return self.graph.nodes()

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
        return list(networkx.ancestors(self.graph, node))

    def map_class_to_id(self, x):
        """
        Map ontology class to ID.
        :param x:
        :return:
        """
        assert False  # ToDo

    def map_id_to_class(self, x):
        """
        Map ontology ID to class.
        :param x:
        :return:
        """
        assert False  # ToDo

    def find_nodes_fuzzy(
            self,
            source,
            match_only: bool = False,
            include_synonyms: bool = True,
            constrain_by_anatomy: Union[str, None] = None,
            omit_list: list = [],
            n_suggest: int = 10,
    ) -> pd.DataFrame:
        """
        Map free text node names to ontology node names via fuzzy string matching.

        If this function does not yield good matches, consider querying this web interface:
        https://www.ebi.ac.uk/ols/index

        :param source: Free text node labels which are to be matched to ontology nodes.
        :param match_only: Whether to include strict matches only in output.
        :param include_synonyms: Whether to include synonyms of nodes in string search.
        :param constrain_by_anatomy: Whether to require suggestions to be within a target anatomy defined within UBERON.
        :param omit_list: Free text node labels to omit in map.
        :param n_suggest: Number of cell types to suggest.
        :return: Table with source and target node names. Columns: "source", "target"
        """
        from fuzzywuzzy import fuzz
        matches = []
        nodes = [(k, v) for k, v in self.graph.nodes.items()]
        include = []
        if isinstance(source, pd.DataFrame):
            source = list(zip(source.iloc[:, 0].values, source.iloc[:, 1].values))
        for x in source:
            if not isinstance(x, list) and not isinstance(x, tuple):
                x = [x, "nan"]
            scores = np.array([
                np.max([
                    fuzz.ratio(x[0].lower().strip("'").strip("\""), y[1]["name"].lower())
                ] + [
                    fuzz.ratio(x[0].lower().strip("'").strip("\"").strip("]").strip("["), yy.lower())
                    for yy in y[1]["synonym"]
                ]) if "synonym" in y[1].keys() and include_synonyms else
                np.max([
                    fuzz.ratio(x[0].lower().strip("'").strip("\""), y[1]["name"].lower())
                ])
                for y in nodes
            ])
            include.append(x[0].lower().strip("'").strip("\"") not in omit_list)
            if match_only:
                matches.append(np.any(scores == 100))  # perfect match
            else:
                if np.any(scores == 100):
                    matches.append([nodes[i][1]["name"] for i in np.where(scores == 100)[0]])
                else:
                    if constrain_by_anatomy is not None:
                        # Select up to 5 nodes which match the anatomical constraint:
                        # ToDo: need to check if constrain_by_anatomy is a parent or child node of:
                        # v['relationship'] = ['part_of UBERON:0001885']
                        # This then implies that the suggested cell type occurs in
                        # a) parent -> a more general setting across anatomies from which one was sampled
                        # b) child -> a sub anatomy of the sampled tissue.
                        matchesi = [
                            nodes[i][1]["name"]
                            for i in np.argsort(scores)
                            if "def" in nodes[i][1].keys() and constrain_by_anatomy in nodes[i][1]["def"]
                        ][-5:][::-1]
                        # Select best remaining matches until n_suggests:
                        matchesi = matchesi + [
                            nodes[i][1]["name"]
                            for i in np.argsort(scores)
                            if nodes[i][1]["name"] not in matchesi
                        ][-np.max(n_suggest - len(matchesi), 0):][::-1]
                    else:
                        # Suggest top 10 hits by string match:
                        matchesi = [nodes[i][1]["name"] for i in np.argsort(scores)[-n_suggest:]][::-1]
                    matches.append(matchesi)
        tab = pd.DataFrame({
            "source": source,
            "target": [":".join(z) for z in matches]
        })
        return tab.loc[include]


class OntologyExtendedObo(OntologyObo):
    """
    Basic .obo ontology extended by additional nodes and edges without breaking DAG.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add_extension(dict_ontology=ONTOLOGIY_EXTENSION_HUMAN)  # ToDo distinguish here

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


class CelltypeUniverse:
    """
    Cell type universe (list) and ontology (hierarchy) container class.


    Basic checks on the organ specific instance are performed in the constructor.
    """
    ontology: OntologyBase
    _target_universe: Union[List[str], None]

    def __init__(self, organism: str, **kwargs):
        """

        :param organism: Organism, defines ontology extension used.
        :param kwargs:
        """
        self.onto = OntologyExtendedObo(**kwargs)
        self._target_universe = None
        self._set_extension(organism=organism)

    def _set_extension(self, organism):
        """

        :param organism: Organism, defines ontology extension used.
        """
        if organism == "human":
            self.onto.add_extension(ONTOLOGIY_EXTENSION_HUMAN)
        elif organism == "mouse":
            self.onto.add_extension(ONTOLOGIY_EXTENSION_MOUSE)
        else:
            raise ValueError(f"organism {organism} not found")

    @property
    def target_universe(self):
        """
        Ontology classes of target universe (understandable cell type names).

        :return:
        """
        return self._target_universe

    @target_universe.setter
    def target_universe(self, x: List[str]):
        # Check that all nodes are valid:
        for xx in x:
            if xx not in self.onto.nodes:
                raise ValueError(f"cell type {xx} was not in ontology")
        # Default universe is the full set of leave nodes of ontology:
        self.target_universe = self.onto.leaves
        self.onto.set_leaves(self.target_universe)

    @property
    def target_universe_ids(self):
        """
        Ontology IDs of target universe (codified cell type names).

        :return:
        """
        return [self.onto.map_class_to_id(x) for x in self._target_universe]

    @property
    def ntypes(self):
        """
        Number of different cell types in target universe.
        """
        return len(self.target_universe)

    def __validate_target_universe_table(self, tab: pd.DataFrame):
        assert len(tab.columns) == 2
        assert tab.columns[0] == "name" and tab.columns[1] == "id"

    def load_target_universe(self, organ):
        """

        :param organ: Anatomic structure to load target universe for.
        :return:
        """
        # ToDo: Use pydoc based query of universes stored in ./target_universes/..
        tab = None
        self.__validate_target_universe_table(tab=tab)
        self.target_universe = None  # ToDo

    def read_target_universe_csv(self, fn):
        """

        :param fn: File containing target universe.
        :return:
        """
        tab = pd.read_csv(fn)
        self.__validate_target_universe_table(tab=tab)
        self.target_universe = tab["name"].values

    def map_to_target_leaves(
            self,
            nodes: List[str],
            return_type: str = "elements"
    ):
        """
        Map a given list of nodes to leave nodes defined for this ontology.
        :param nodes:
        :param return_type:

            "elements": names of mapped leave nodes
            "idx": indices in leave note list of of mapped leave nodes
        :return:
        """
        return [self.onto.map_to_leaves(x, return_type=return_type) for x in nodes]
