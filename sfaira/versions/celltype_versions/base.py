import abc
import networkx
import numpy as np
import obonet
import pandas as pd
from typing import Dict, List, Tuple, Union
import warnings

from sfaira.versions.celltype_versions.extensions import ONTOLOGIY_EXTENSION_HUMAN, ONTOLOGIY_EXTENSION_MOUSE


class OntologyObo:

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
    def nodes(self):
        return list(self.graph.nodes.items())

    @property
    def nodes_dict(self):
        return self.graph.nodes.items()

    @property
    def node_names(self):
        return [x["name"] for x in self.graph.nodes.values()]

    @property
    def node_ids(self):
        return list(self.graph.nodes())

    def id_from_name(self, x: str):
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
        return list(networkx.ancestors(self.graph, node))

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

    def validate_node(self, x: str):
        if x not in self.node_names:
            suggestions = self.map_node_suggestion(x=x, include_synonyms=False)
            raise ValueError(f"Node label {x} not found. Did you mean any of {suggestions}?")

    @abc.abstractmethod
    def synonym_node_properties(self) -> List[str]:
        pass


class OntologyExtendedObo(OntologyObo):
    """
    Basic .obo ontology extended by additional nodes and edges without breaking DAG.
    """

    def __init__(self, obo, **kwargs):
        super().__init__(obo=obo, **kwargs)
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
            **kwargs
    ):
        super().__init__(obo="http://purl.obolibrary.org/obo/cl.obo")

        # Clean up nodes:
        nodes_to_delete = []
        for k, v in self.graph.nodes.items():
            if "namespace" not in v.keys() or v["namespace"] != "cell":
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


class CelltypeUniverse:
    """
    Cell type universe (list) and ontology (hierarchy) container class.


    Basic checks on the organ specific instance are performed in the constructor.
    """
    ontology: OntologyCelltypes
    _target_universe: Union[List[str], None]

    def __init__(self, organism: str, **kwargs):
        """

        :param organism: Organism, defines ontology extension used.
        :param kwargs:
        """
        self.onto_cl = OntologyCelltypes(**kwargs)
        self.onto_anatomy = OntologyUberon(**kwargs)
        self._target_universe = None
        self._set_extension(organism=organism)

    def _set_extension(self, organism):
        """

        :param organism: Organism, defines ontology extension used.
        """
        if organism == "human":
            self.onto_cl.add_extension(ONTOLOGIY_EXTENSION_HUMAN)
        elif organism == "mouse":
            self.onto_cl.add_extension(ONTOLOGIY_EXTENSION_MOUSE)
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
            if xx not in self.onto_cl.nodes:
                raise ValueError(f"cell type {xx} was not in ontology")
        # Default universe is the full set of leave nodes of ontology:
        self.target_universe = self.onto_cl.leaves
        self.onto_cl.set_leaves(self.target_universe)

    @property
    def target_universe_ids(self):
        """
        Ontology IDs of target universe (codified cell type names).

        :return:
        """
        return [self.onto_cl.map_class_to_id(x) for x in self._target_universe]

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
        return [self.onto_cl.map_to_leaves(x, return_type=return_type) for x in nodes]

    def prepare_celltype_map_fuzzy(
            self,
            source,
            match_only: bool = False,
            include_synonyms: bool = True,
            anatomical_constraint: Union[str, None] = None,
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
        :param anatomical_constraint: Whether to require suggestions to be within a target anatomy defined within UBERON.
        :param omit_list: Free text node labels to omit in map.
        :param n_suggest: Number of cell types to suggest.
        :return: Table with source and target node names. Columns: "source", "target"
        """
        from fuzzywuzzy import fuzz
        matches = []
        nodes = self.onto_cl.nodes
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
                    if anatomical_constraint is not None:
                        # Check that anatomical constraint is a term in UBERON and get UBERON ID:
                        anatomical_constraint_id = self.onto_anatomy.id_from_name(anatomical_constraint)
                        # Select up to 5 nodes which match the anatomical constraint:
                        # The entries look as follows:
                        # node.value['relationship'] = ['part_of UBERON:0001885']
                        # Find nodes that can be matched to UBERON:
                        anatomical_subselection = [
                            "relationship" in y[1].keys() and
                            np.any(["part_of UBERON" in yy for yy in y[1]["relationship"]]) and
                            np.any([
                                yy.split("part_of ")[-1] in self.onto_anatomy.node_ids
                                for yy in y[1]["relationship"]
                            ])
                            for y in nodes
                        ]
                        uberon_ids = [
                            y[1]["relationship"][
                                np.where(["part_of UBERON" in yy for yy in y[1]["relationship"]])[0][0]
                            ].split("part_of ")[1]
                            if z else None
                            for y, z in zip(nodes, anatomical_subselection)
                        ]
                        # Check relationship in UBERON. Select for:
                        # a) parent -> a more general setting across anatomies from which one was sampled
                        # b) child -> a sub anatomy of the sampled tissue.
                        # Check this by checking if one is an ancestor of the other:
                        anatomical_subselection = [
                            z and (
                                anatomical_constraint_id in self.onto_anatomy.get_ancestors(node=y) or
                                y in self.onto_anatomy.get_ancestors(node=anatomical_constraint_id)
                            )
                            for y, z in zip(uberon_ids, anatomical_subselection)
                        ]
                        # Iterate over nodes sorted by string match score and masked by constraint:
                        matchesi = [
                            nodes[i][1]["name"]
                            for i in np.argsort(scores)
                            if anatomical_subselection[i]
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
