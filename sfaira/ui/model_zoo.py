import numpy as np
import pandas as pd
from typing import List, Union

from sfaira.versions.metadata import CelltypeUniverse
from sfaira.consts import OCS
from sfaira.versions.topologies import TopologyContainer, TOPOLOGIES


class ModelZoo:

    """
    Model zoo class.
    """

    _model_id: Union[str, None]
    available_model_ids: Union[list, None]
    celltypes: Union[CelltypeUniverse, None]
    topology_container: Union[None, TopologyContainer]
    zoo: Union[dict, None]

    TOPOLOGIES = TOPOLOGIES
    TOPOLOGY_CONTAINER_CLASS = TopologyContainer

    def __init__(
            self,
            model_lookuptable: Union[None, pd.DataFrame] = None,
            model_class: Union[str, None] = None,
    ):
        """
        :param model_lookuptable: model_lookuptable.
        :param model_class: Model class to subset to.
        """
        self._ontology_container_sfaira = OCS
        self._model_id = None
        self.topology_container = None

        if model_lookuptable is not None:  # check if models in repository
            self._load_model_ids(model_ids=model_lookuptable['model_id'].values, model_class=model_class)
            self._construct_zoo_from_model_ids()
        else:
            self.zoo = None
            self.available_model_ids = None

    def _load_model_ids(
            self,
            model_ids,
            model_class: Union[str, None] = None,
    ):
        """
        Load model ids based on models available in model lookup tables.

        :param model_ids: Table listing all available model_ids.
        :param model_class: Model class to subset to
        """
        self.available_model_ids = [x for x in model_ids if (x.split('_')[0] == model_class or model_class is None)]

    def _construct_zoo_from_model_ids(self):
        """
        Load model zoo based on models available model_ids.
        """
        id_df = pd.DataFrame(
            [i.split('_')[1:3] for i in self.available_model_ids],
            columns=['name', 'organisation']
        )
        orgs = np.unique(id_df['organisation'])
        zoo = dict.fromkeys(orgs)
        for o in orgs:
            id_df_o = id_df[id_df['organisation'] == o]
            name = np.unique(id_df_o['name'])
            zoo[o] = dict.fromkeys(name)
        self.zoo = zoo

    @staticmethod
    def _order_versions(
            versions: List[str]
    ):
        """
        Order list of versions of the form 'vx.y.z' from newest to oldest.

        :param versions: Unordered list of version IDs.
        :return: Ordered list of versions.
        """
        versions.sort(key=lambda s: [int(u) for u in s.split('.')])
        return versions

    def topology(
            self,
            model_type: str,
            organisation: str
    ) -> List[str]:
        """
        Return list of available model topologies that trained by a given organisation, and a given model

        :param model_type: Identifier of model_type to show versions for.
        :param organisation: Identifier of organisation to show versions for.
        :return: List of versions available.
        """
        assert organisation in self.zoo.keys(), "organisation requested was not found in zoo"
        assert model_type in self.zoo[organisation].keys(), \
            "model_type requested was not found in zoo"
        return self.zoo[organisation][model_type]

    def versions(
            self,
            model_type: str,
            organisation: str,
            model_topology: str
    ) -> List[str]:
        """
        Return list of available model versions of a given organisation for a given organism and organ and model.

        :param model_type: Identifier of model_type to show versions for.
        :param organisation: Identifier of organisation to show versions for.
        :param model_topology: Identifier of model_topology to show versions for.
        :return: List of versions available.
        """
        assert organisation in self.zoo.keys(), "organisation requested was not found in zoo"
        assert model_type in self.zoo[organisation].keys(), \
            "model_type requested was not found in zoo"
        assert model_topology in self.zoo[organisation][model_type].keys(), \
            "model_topology requested was not found in zoo"
        return self.zoo[organisation][model_type][model_topology]

    @property
    def model_hyperparameters(self) -> dict:
        assert self.topology_container is not None, "set model_id first"
        return self.topology_container.topology["hyper_parameters"]

    @property
    def celltypes(self):
        assert self.topology_container is not None, "set model_id first"
        return self.topology_container.topology["output"]["targets"]

    @celltypes.setter
    def celltypes(self, x: List):
        assert self.topology_container is not None, "set model_id first"
        self.topology_container.topology["output"]["targets"] = x

    @property
    def model_id(self):
        return self._model_id

    @model_id.setter
    def model_id(self, x: str):
        """
        Set model ID to a manually supplied ID and automatically set topology container.

        :param x: Model ID to set. Format: modelclass_organism-organ-modeltype-topology-version_organisation
        """
        assert self.available_model_ids is None or x in self.available_model_ids,\
            f"{x} not found in available_model_ids, please check available models using ModelZoo.available_model_ids"
        assert len(x.split('_')) == 3, f'model_id {x} is invalid'
        self._model_id = x
        self.topology_container = self.TOPOLOGY_CONTAINER_CLASS(
            topology=self.TOPOLOGIES[self.model_organism][self.model_class][self.model_type][self.model_topology],
            topology_id=self.model_version
        )

    @property
    def model_class(self):
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[0]

    @property
    def model_name(self):
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[1]

    @property
    def model_organism(self):
        # TODO: this relies on theislab model_name formatting
        assert self.model_id is not None, "set model_id first"
        return {
            "homosapiens": "Homo sapiens",
            "musmusculus": "Mus musculus",
            "human": "Homo sapiens",  # necessary for old sfaira model uploads
            "mouse": "Mus musculus",  # necessary for old sfaira model uploads
        }[self.model_id.split('_')[1].split("-")[0]]

    @property
    def model_organ(self):
        # TODO: this relies on theislab model_name formatting
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[1].split("-")[1]

    @property
    def model_type(self):
        # TODO: this relies on theislab model_name formatting
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[1].split("-")[2]

    @property
    def model_topology(self):
        # TODO: this relies on theislab model_name formatting
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[1].split("-")[3]

    @property
    def model_version(self):
        # TODO: this relies on theislab model_name formatting
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[1].split("-")[4]

    @property
    def organisation(self):
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[2]
