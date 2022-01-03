import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Union

from sfaira.versions.metadata import OntologyCl, OntologyList, OntologyUberon

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

    def __init__(self, cl: Union[OntologyCl, OntologyList], uberon: OntologyUberon, **kwargs):
        self.onto_cl = cl
        self.onto_uberon = uberon
        self._target_universe = None

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
