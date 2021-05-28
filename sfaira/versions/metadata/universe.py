import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Union

from sfaira.versions.metadata import OntologyCl, OntologyUberon
from sfaira.versions.metadata.extensions import ONTOLOGIY_EXTENSION

TARGET_UNIVERSE_KEY_NAME = "name"
TARGET_UNIVERSE_KEY_ID = "id"


class CelltypeUniverse:
    """
    Cell type universe (list) and ontology (hierarchy) container class.


    Basic checks on the organ specific instance are performed in the constructor.
    """
    onto_cl: OntologyCl
    onto_uberon: OntologyUberon
    _target_universe: Union[List[str], None]

    def __init__(self, cl: OntologyCl, uberon: OntologyUberon, **kwargs):
        self.onto_cl = cl
        self.onto_uberon = uberon
        self._target_universe = None
        self._set_extension()

    def _set_extension(self):
        self.onto_cl.add_extension(ONTOLOGIY_EXTENSION)

    def __validate_target_universe_table(self, tab: pd.DataFrame):
        assert len(tab.columns) == 2
        assert tab.columns[0] == "name" and tab.columns[1] == "id"

    def load_target_universe(self, fn):
        """

        :param fn: .csv file containing target universe.
        :return:
        """
        tab = pd.read_csv(fn, sep="\t", index_col=None)
        self.__validate_target_universe_table(tab=tab)
        self.onto_cl.leaves = tab["name"].values

    def write_target_universe(
            self,
            fn,
            x: List[str],
    ):
        """

        :param fn: .csv file containing target universe.
        :param x: Nodes that make up target universe.
        :return:
        """
        tab = pd.DataFrame({
            TARGET_UNIVERSE_KEY_NAME: self.onto_cl.convert_to_name(x),
            TARGET_UNIVERSE_KEY_ID: self.onto_cl.convert_to_id(x),
        })
        self.__validate_target_universe_table(tab=tab)
        tab.to_csv(path_or_buf=fn, sep="\t", index=False)

    def prepare_celltype_map_fuzzy(
            self,
            source,
            match_only: bool = False,
            include_synonyms: bool = True,
            anatomical_constraint: Union[str, None] = None,
            choices_for_perfect_match: bool = True,
            omit_list: list = [],
            omit_target_list: list = ["cell"],
            n_suggest: int = 4,
            threshold_for_partial_matching: float = 90.,
    ) -> Tuple[
        List[Dict[str, Union[List[str], str]]],
        List[bool]
    ]:
        """
        Map free text node names to ontology node names via fuzzy string matching and return as list

        If this function does not yield good matches, consider querying this web interface:
        https://www.ebi.ac.uk/ols/index

        Search strategies:

            - exact_match: Only exact string matches to name or synonym in ontology. This is the only strategy that is
                enabled if match_only is True.
            - lenient_match: Fuzzy string matches to name or synonym in ontology based on ratio of match errors
                ((fuzz.ratio).
            - very_lenient_match: Fuzzy string matches to name or synonym in ontology based on ratio of matches
                characters from query (fuzz.partial_ratio)

        Search strategies with anatomical constraints:
        An anatomic constraint is a name of an anatomical structure that can be mapped to UBERON.

            - anatomic_onotolgy_match:
                We select cell types expected in this UBERON clade based on the link between CL and UBERON.
            - anatomic_string_match:
                We perform an additional fuzzy string matching with the anatomical structure added to the proposed
                label. This is often beneficial because analysts do not always prefix such extension (e.g. pancreatic)
                to the free text cell type labels if the entire sample consists only of cells from this anatomical
                structure. Note that if the maps from 1) were perfect, this would not be necessary. In practice, we
                find this to still recover some hits that are otherwise missed.

        Note that matches are shadowed in lower priorty strategies, ie a perfect match will not show up in the list
        of hits of any other strategy.

        :param source: Free text node labels which are to be matched to ontology nodes.
        :param match_only: Whether to include strict matches only in output.
        :param include_synonyms: Whether to include synonyms of nodes in string search.
        :param anatomical_constraint: Whether to require suggestions to be within a target anatomy defined
            within UBERON.
        :param choices_for_perfect_match: Whether to give additional matches if a perfect match was found and an
            anatomical_constraint is not not defined. This is overridden by match_only.
        :param omit_list: Free text node labels to omit in map.
        :param omit_target_list: Ontology nodes to not match to.
        :param n_suggest: Number of cell types to suggest per search strategy.
        :param threshold_for_partial_matching: Maximum fuzzy match score below which lenient matching (ratio) is
            extended through partial_ratio.
        :return: Tuple

            - List with matches for each source, each entry is a dictionary,
                of lists of search strategies named by strategy name. If a search strategy yields perfect matches, it
                does not return a list of strings but just a single string.
            - List with boolean indicator whether or not this output should be reported.
        """
        from fuzzywuzzy import fuzz
        matches = []
        nodes = self.onto_cl.nodes
        nodes = [x for x in nodes if x[1]["name"] not in omit_target_list]
        include_terms = []
        if isinstance(source, pd.DataFrame):
            source = list(zip(source.iloc[:, 0].values, source.iloc[:, 1].values))
        for x in source:
            if not isinstance(x, list) and not isinstance(x, tuple):
                x = [x, "nan"]
            term = x[0].lower().strip("'").strip("\"").strip("'").strip("\"").strip("]").strip("[")
            # Test for perfect string matching:
            scores_strict = np.array([
                np.max(
                    [
                        100 if term == y[1]["name"].lower() else 0
                    ] + [
                        100 if term == yy.lower() else 0
                        for yy in y[1]["synonym"]
                    ]
                ) if "synonym" in y[1].keys() and include_synonyms else 100 if term == y[1]["name"].lower() else 0
                for y in nodes
            ])
            # Test for partial string matching:
            # fuzz ratio and partial_ratio capture different types of matches well, we use both here and decide below
            # which scores are used in which scenario defined through the user input.
            # Formatting of synonyms: These are richly annotated, we strip references following after either:
            # BROAD, EXACT
            # in the synonym string and characters: "'

            def synonym_string_processing(y):
                return y.lower().split("broad")[0].split("exact")[0].lower().strip("'").strip("\"").split("\" ")[0]

            scores_lenient = np.array([
                np.max([fuzz.ratio(term, y[1]["name"].lower())] + [
                    fuzz.ratio(term, synonym_string_processing(yy))
                    for yy in y[1]["synonym"]
                ]) if "synonym" in y[1].keys() and include_synonyms else
                fuzz.ratio(term, y[1]["name"].lower())
                for y in nodes
            ])
            scores_very_lenient = np.array([
                np.max([fuzz.partial_ratio(term, y[1]["name"].lower())] + [
                    fuzz.partial_ratio(term, synonym_string_processing(yy))
                    for yy in y[1]["synonym"]
                ]) if "synonym" in y[1].keys() and include_synonyms else
                fuzz.partial_ratio(term, y[1]["name"].lower())
                for y in nodes
            ])
            include_terms.append(term not in omit_list)
            if match_only and not anatomical_constraint:
                # Explicitly trying to report perfect matches (match_only is True).
                matches.append({"perfect_match": [nodes[i][1]["name"] for i in np.where(scores_strict == 100)[0]][0]})
            else:
                matches_i = {}
                if np.any(scores_strict == 100) and not anatomical_constraint:
                    # Perfect match and not additional information through anatomical_constraint, ie no reason to assume
                    # that the user is not looking for this hit.
                    matches_i.update({
                        "perfect_match": [nodes[i][1]["name"] for i in np.where(scores_strict == 100)[0]][0]
                    })
                    if choices_for_perfect_match:
                        matches_i.update({
                            "lenient_match": [
                                nodes[i][1]["name"] for i in np.argsort(scores_lenient)[::-1]
                                if not np.any([nodes[i][1]["name"] in v for v in matches_i.values()])
                            ][:n_suggest]
                        })
                        if np.max(scores_lenient) < threshold_for_partial_matching:
                            matches_i.update({
                                "very_lenient_match": [
                                    nodes[i][1]["name"]
                                    for i in np.argsort(scores_very_lenient)[::-1]
                                    if not np.any([nodes[i][1]["name"] in v for v in matches_i.values()])
                                ][:n_suggest]
                            })
                else:
                    if anatomical_constraint is not None:
                        # Use anatomical constraints two fold:
                        # 1. Select cell types that are in the correct ontology.
                        # 2. Run a second string matching with the anatomical word included.

                        # 1. Select cell types that are in the correct ontology.
                        # Check that anatomical constraint is a term in UBERON and get UBERON ID:
                        anatomical_constraint_id = self.onto_uberon.convert_to_id(anatomical_constraint)
                        # Select up to 5 nodes which match the anatomical constraint:
                        # The entries look as follows:
                        # node.value['relationship'] = ['part_of UBERON:0001885']
                        # Find nodes that can be matched to UBERON:
                        anatomical_subselection = [
                            "relationship" in y[1].keys() and
                            np.any(["part_of UBERON" in yy for yy in y[1]["relationship"]]) and
                            np.any([
                                yy.split("part_of ")[-1] in self.onto_uberon.node_ids
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
                                anatomical_constraint_id in self.onto_uberon.get_ancestors(node=y) or
                                y in self.onto_uberon.get_ancestors(node=anatomical_constraint_id)
                            )
                            for y, z in zip(uberon_ids, anatomical_subselection)
                        ]
                        # Iterate over nodes sorted by string match score and masked by constraint:
                        matches_i.update({
                            "anatomic_onotolgy_match": [
                                nodes[i][1]["name"]
                                for i in np.argsort(scores_lenient)
                                if anatomical_subselection[i] and not
                                np.any([nodes[i][1]["name"] in v for v in matches_i.values()])
                            ][-n_suggest:][::-1]
                        })

                        # 2. Run a second string matching with the anatomical word included.
                        modified_term = anatomical_constraint + " " + x[0].lower().strip("'").strip("\"").strip("]"). \
                            strip("[")
                        scores_anatomy = np.array([
                            np.max([
                                fuzz.partial_ratio(modified_term, y[1]["name"].lower())
                            ] + [
                                fuzz.partial_ratio(modified_term, synonym_string_processing(yy))
                                for yy in y[1]["synonym"]
                            ]) if "synonym" in y[1].keys() and include_synonyms else
                            fuzz.partial_ratio(modified_term, y[1]["name"].lower())
                            for y in nodes
                        ])
                        matches_i.update({
                            "anatomic_string_match": [
                                nodes[i][1]["name"] for i in np.argsort(scores_anatomy)
                                if nodes[i][1]["name"] and not
                                np.any([nodes[i][1]["name"] in v for v in matches_i.values()])
                            ][-n_suggest:][::-1]
                        })

                        # Select best overall matches based on lenient and strict matching:
                        matches_i.update({
                            "perfect_match": [
                                nodes[i][1]["name"]
                                for i in np.argsort(scores_strict)[::-1]
                            ][:n_suggest]
                        })
                        matches_i.update({
                            "lenient_match": [
                                nodes[i][1]["name"]
                                for i in np.argsort(scores_lenient)[::-1]
                                if not np.any([nodes[i][1]["name"] in v for v in matches_i.values()])
                            ][:n_suggest]
                        })
                        if np.max(scores_lenient) < threshold_for_partial_matching:
                            matches_i.update({
                                "very_lenient_match": [
                                    nodes[i][1]["name"]
                                    for i in np.argsort(scores_very_lenient)[::-1]
                                    if not np.any([nodes[i][1]["name"] in v for v in matches_i.values()])
                                ][:n_suggest]
                            })
                    else:
                        # Suggest top hits by string match:
                        matches_i.update({
                            "lenient_match": [
                                nodes[i][1]["name"] for i in np.argsort(scores_lenient)[::-1]
                            ][:n_suggest]
                        })
                        if np.max(scores_lenient) < threshold_for_partial_matching:
                            matches_i.update({
                                "very_lenient_match": [
                                    nodes[i][1]["name"]
                                    for i in np.argsort(scores_very_lenient)[::-1]
                                    if not np.any([nodes[i][1]["name"] in v for v in matches_i.values()])
                                ][:n_suggest]
                            })
                matches.append(matches_i)
        return matches, include_terms

    def prepare_celltype_map_tab(
            self,
            source,
            match_only: bool = False,
            include_synonyms: bool = True,
            anatomical_constraint: Union[str, None] = None,
            omit_list: list = [],
            n_suggest: int = 10,
            separator_suggestions: str = ":",
            separator_groups: str = ":|||:",
    ) -> pd.DataFrame:
        """
        Map free text node names to ontology node names via fuzzy string matching and return as matching table.

        :param source: Free text node labels which are to be matched to ontology nodes.
        :param match_only: Whether to include strict matches only in output.
        :param include_synonyms: Whether to include synonyms of nodes in string search.
        :param anatomical_constraint: Whether to require suggestions to be within a target anatomy defined within
            UBERON.
        :param omit_list: Free text node labels to omit in map.
        :param n_suggest: Number of cell types to suggest per search strategy.
        :param separator_suggestions: String separator for matches of a single strategy in output target column.
        :param separator_groups: String separator for search strategy grouped matches in output target column.
        :return: Table with source and target node names. Columns: "source", "target"
        """
        matches, include_terms = self.prepare_celltype_map_fuzzy(
            source=source,
            match_only=match_only,
            include_synonyms=include_synonyms,
            anatomical_constraint=anatomical_constraint,
            choices_for_perfect_match=False,
            omit_list=omit_list,
            n_suggest=n_suggest,
        )
        tab = pd.DataFrame({
            "source": source,
            "target": [
                separator_groups.join([
                    separator_suggestions.join(v)
                    if isinstance(v, list) else v
                    for v in x.values()
                ])
                for x in matches
            ]
        })
        return tab.loc[include_terms]
