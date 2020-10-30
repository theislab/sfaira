import abc
import networkx
import numpy as np
import obonet
import pandas as pd
from typing import List, Tuple, Union


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

    def __init__(self, obo: str = "http://purl.obolibrary.org/obo/cl.obo"):
        self.graph = obonet.read_obo(obo)
        assert networkx.is_directed_acyclic_graph(self.graph)

    def set_leaves(self, nodes: list = None):
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

    def fuzzymatch_nodes(
            self,
            source,
            match_only: bool = False,
            include_old: bool = False,
            include_synonyms: bool = True,
            remove: list = []
    ):
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
                ] + ([
                    fuzz.ratio(x[0].lower().strip("'").strip("\"").strip("]").strip("["), yy.lower())
                    for yy in y[1]["synonym"]
                ] if "synonym" in y[1].keys() and include_synonyms else []))
                for y in nodes
            ])
            include.append(x[0].lower().strip("'").strip("\"") not in remove)
            if match_only:
                matches.append(np.any(scores == 100))  # perfect match
            else:
                if np.any(scores == 100):
                    matches.append([(nodes[i][1]["name"], nodes[i][0]) for i in np.where(scores == 100)[0]])
                else:
                    matchesi = [(
                        nodes[i][1]["name"] + "[" + ";".join([
                            yy.strip("'").strip("\"").strip("]").strip("[")
                            for yy in nodes[i][1]["synonym"]
                        ]) + "}"
                        if "synonym" in nodes[i][1].keys() and include_synonyms else nodes[i][1]["name"],
                        nodes[i][0]
                    ) for i in np.argsort(scores)[-10:]]
                    if include_old:
                        matchesi = matchesi + [(x[0].upper(), x[1])]
                    matches.append(matchesi)
        if match_only:
            tab = pd.DataFrame({"name": source, "matched": matches})
        else:
            tab = pd.DataFrame({"name,id": [" ".join([",".join(zz) for zz in z]) for z in matches]})
        return tab.loc[include]


class CelltypeVersionsBase:
    """
    Versioned cell type universe (list) and ontology (hierarchy) container class.

    This class is subclassed once for each anatomical structure (organ).
    Cell type versions take the form x.y:
        - x is a major version and is incremented if the cell type identities or number changes.
        - y is a minor version and is incremented if cell type annotation such as names are altered without changing
            the described cell type or adding cell types.

    Basic checks on the organ specific instance are performed in the constructor.
    """
    celltype_universe: dict
    ontology: dict
    version: str

    def __init__(self, **kwargs):
        # Check that versions are consistent.
        if not list(self.celltype_universe.keys()) == list(self.ontology.keys()):
            raise ValueError(
                "error in matching versions of cell type universe and ontology in %s" %
                type(self)
            )
        # Check that ontology terms are unique also between ontologies
        if np.sum([len(x) for x in self.ontology.values()]) != \
            len(np.unique(np.array([list(x) for x in self.ontology.values()]))):
            raise ValueError(
                "duplicated ontology terms found between ontologies in %s" %
                type(self)
            )

    @property
    def versions(self):
        """
        Available cell type universe versions loaded in this instance.

        :return:
        """
        return self.celltype_universe.keys()

    def _check_version(self, version: str):
        if version not in self.celltype_universe.keys():
            raise ValueError("Version %s not found. Check self.version for available versions." % version)

    def set_version(
            self,
            version: str
    ):
        """
        Set a cell type universe version for this instance.

        :param version: Full version string "a.b.c" or celltype version "a".
        :return:
        """
        if len(version.split(".")) == 3:
            version = version.split(".")[0]
            self._check_version(version=version)
            self.version = version
        elif len(version.split(".")) == 1:
            self._check_version(version=version)
            self.version = version
        else:
            raise ValueError("version supplied should be either in format `a.b.c` or `a`")


    @property
    def ids(self):
        """
        List of all human understandable cell type names of this instance.

        :return:
        """
        return self._ids(self.version)

    def _ids(self, version: str):
        return [x[0] for x in self.celltype_universe[version]]

    @property
    def ontology_ids(self):
        """
        List of all cell type IDs (based on an ontology ID scheme) of this instance.

        :return:
        """
        return self._ontology_ids(self.version)

    def _ontology_ids(self, version: str):
        return [x[1] for x in self.celltype_universe[version]]

    @property
    def ntypes(self):
        """
        Number of different cell types in this instance.

        :return:
        """
        return self._ntypes(self.version)

    def _ntypes(self, version: str):
        return len(self.celltype_universe[version])

    def read_csv(self, fn) -> List[Tuple[str, str]]:
        tab = pd.read_csv(fn)
        return [(x, y) for x, y in zip(tab["name"].values, tab["id"].values)]

    def map_to_leaves(
            self,
            nodes: List[str],
            ontology: str = "custom",
            ontology_id: str = "names",
            return_type: str = "elements"
    ):
        """
        Map a given list of nodes to leave nodes defined for this ontology
        :param nodes:
        :param ontology:
        :param return_type:

            "elements": names of mapped leave nodes
            "idx": indicies in leave note list of of mapped leave nodes
        :return:
        """
        if ontology == "custom":
            onto = OntologyDict(self.ontology[self.version][ontology_id])
        elif ontology == "cl":
            onto = OntologyObo()
        else:
            assert False
        onto.set_leaves(self.celltype_universe[self.version])
        return [onto.map_to_leaves(x, return_type=return_type) for x in nodes]
