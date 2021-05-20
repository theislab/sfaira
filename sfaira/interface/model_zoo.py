import abc
import numpy as np
import pandas as pd
from typing import List, Union

from sfaira.versions.metadata import CelltypeUniverse
from sfaira.consts import OntologyContainerSfaira
from sfaira.versions.topologies import TopologyContainer, TOPOLOGIES


class ModelZoo(abc.ABC):
    """
    Model ontology base class.
    """
    topology_container: TopologyContainer
    ontology: dict
    _model_id: Union[str, None]
    celltypes: Union[CelltypeUniverse, None]

    def __init__(
            self,
            model_lookuptable: Union[None, pd.DataFrame] = None,
            model_class: Union[str, None] = None,
    ):
        """
        :param model_lookuptable: model_lookuptable.
        :param model_class: Model class to subset to.
        """
        self._ontology_container_sfaira = OntologyContainerSfaira()
        if model_lookuptable is not None:  # check if models in repository
            self.ontology = self.load_ontology_from_model_ids(model_ids=model_lookuptable['model_id'].values,
                                                              model_class=model_class)
        self._model_id = None
        self.celltypes = None

    @staticmethod
    def load_ontology_from_model_ids(
            model_ids,
            model_class: Union[str, None] = None,
    ) -> dict:
        """
        Load model ontology based on models available in model lookup tables.

        :param model_ids: Table listing all available model_ids.
        :param model_class: Model class to subset to.
        :return: Dictionary formatted ontology.
        """

        ids = [x for x in model_ids if (x.split('_')[0] == model_class or model_class is None)]
        id_df = pd.DataFrame(
            [i.split('_')[1:3] for i in ids],
            columns=['name', 'organisation']
        )
        model = np.unique(id_df['name'])
        ontology = dict.fromkeys(model)
        for m in model:
            id_df_m = id_df[id_df['name'] == m]
            orga = np.unique(id_df_m['organisation'])
            ontology[m] = dict.fromkeys(orga)
        return ontology

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
        assert model_type in self.ontology.keys(), "model_type requested was not found in ontology"
        assert organisation in self.ontology[model_type].keys(), \
            "organisation requested was not found in ontology"
        return self.ontology[model_type][organisation]

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
        assert model_type in self.ontology.keys(), "model_type requested was not found in ontology"
        assert organisation in self.ontology[model_type].keys(), \
            "organisation requested was not found in ontology"
        assert model_topology in self.ontology[model_type][organisation].keys(), \
            "model_topology requested was not found in ontology"
        return self.ontology[model_type][organisation][model_topology]

    @property
    def model_hyperparameters(self) -> dict:
        assert self.topology_container is not None
        return self.topology_container.topology["hyper_parameters"]

    @property
    def topology_container(self) -> TopologyContainer:
        # TODO: this ID decomposition to organism is custom to the topologies handled in this package.
        organism = self.model_name.split("-")[0]
        return TopologyContainer(
            topology=TOPOLOGIES[organism][self.model_class][self.model_type][self.model_topology],
            topology_id=self.model_version
        )

    @property
    def model_id(self):
        return self._model_id

    @model_id.setter
    def model_id(self, x: str):
        """
        Set model ID to a manually supplied ID.

        :param x: Model ID to set. Format: modelclass_organism-organ-modeltype-topology-version_organisation
        """
        assert len(x.split('_')) == 3, f'model_id {x} is invalid'
        self._model_id = x

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
        # TODO: this is a custom name ontology
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[1].split("-")[0]

    @property
    def model_organ(self):
        # TODO: this is a custom name ontology
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[1].split("-")[1]

    @property
    def model_type(self):
        # TODO: this is a custom name ontology
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[1].split("-")[2]

    @property
    def model_topology(self):
        # TODO: this is a custom name ontology
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[1].split("-")[3]

    @property
    def model_version(self):
        # TODO: this is a custom name ontology
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[1].split("-")[4]

    @property
    def organisation(self):
        assert self.model_id is not None, "set model_id first"
        return self.model_id.split('_')[2]
