import abc
try:
    import kipoi
except ImportError:
    kipoi = None
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
    model_id: Union[str, None]
    model_class: Union[str, None]
    model_class: Union[str, None]
    model_type: Union[str, None]
    model_topology: Union[str, None]
    model_version: Union[str, None]
    celltypes: Union[CelltypeUniverse, None]

    def __init__(
            self,
            model_lookuptable: Union[None, pd.DataFrame] = None
    ):
        """
        :param model_lookuptable: model_lookuptable.
        """
        self._ontology_container_sfaira = OntologyContainerSfaira()
        if model_lookuptable is not None:  # check if models in repository
            self.ontology = self.load_ontology_from_model_ids(model_lookuptable['model_id'].values)
        self.model_id = None
        self.model_class = None
        self.model_type = None
        self.organisation = None
        self.model_topology = None
        self.model_version = None
        self.celltypes = None

    @abc.abstractmethod
    def load_ontology_from_model_ids(
            self,
            model_ids
    ):
        pass

    def _order_versions(
            self,
            versions: List[str]
    ):
        """
        Order list of versions of the form 'vx.y.z' from newest to oldest.

        :param versions: Unordered list of version IDs.
        :return: Ordered list of versions.
        """
        versions.sort(key=lambda s: [int(u) for u in s.split('.')])

        return versions

    def set_model_id(
            self,
            model_id: str
    ):
        """
        Set model ID to a manually supplied ID.

        :param model_id: Model ID to set. Format: pipeline_genome_organ_model_organisation_topology_version
        """
        if len(model_id.split('_')) < 6:
            raise RuntimeError(f'Model ID {model_id} is invalid!')
        self.model_id = model_id
        ixs = self.model_id.split('_')
        self.model_class = ixs[0]
        self.model_id = ixs[1]
        self.model_type = ixs[2]
        self.organisation = ixs[3]
        self.model_topology = ixs[4]
        self.model_version = ixs[5]

    def save_weights_to_remote(self, path=None):
        """
        Saves model weights to repository XY.
        Increments 3rd digit of version number.
        Adds model_id to the text file, updates model_index
        """
        raise NotImplementedError()

    def save_weights_to_public(self):
        """
        Saves model weights to cloud under an organization name.
        Increments 2nd digit of version number.
        Adds model_id to the text file, updates model_index
        """
        raise NotImplementedError()

    def call_kipoi(self):
        """
        Returns kipoi_experimental model call from remote directly on local data using kipoi_experimental.

        Runs model defined in self.model_id.
        For this, the remote server associated with the model_id has to be identified via find_remote().

        :return: Predictions
        """
        raise NotImplementedError()

    def models(self) -> List[str]:
        """
        Return list of available models.

        :return: List of models available.
        """
        return self.ontology.keys()

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
        assert self.model_id is not None, "set model_id first"
        # TODO: this ID decomposition to organism is custom to the topologies handled in this package.
        organism = self.model_id.split("-")[0]
        return TopologyContainer(
            topology=TOPOLOGIES[organism][self.model_class][self.model_type][self.model_version],
            topology_id=self.model_version
        )


class ModelZooEmbedding(ModelZoo):

    def load_ontology_from_model_ids(
            self,
            model_ids
    ) -> dict:
        """
        Load model ontology based on models available in model lookup tables.

        :param model_ids: Table listing all available model_ids.
        :return: Dictionary formatted ontology.
        """

        ids = [i for i in model_ids if i.split('_')[0] == 'embedding']
        id_df = pd.DataFrame(
            [i.split('_')[1:6] for i in ids],
            columns=['id', 'model_type', 'organisation', 'model_topology', 'model_version']
        )
        model = np.unique(id_df['model_type'])
        ontology = dict.fromkeys(model)
        for m in model:
            id_df_m = id_df[id_df.model_type == m]
            orga = np.unique(id_df_m['organisation'])
            ontology[m] = dict.fromkeys(orga)
            for org in orga:
                id_df_org = id_df_m[id_df_m.organisation == org]
                topo = np.unique(id_df_org['model_topology'])
                ontology[m][org] = dict.fromkeys(topo)
                for t in topo:
                    id_df_t = id_df_org[id_df_org.model_topology == t]
                    ontology[m][org][t] = id_df_t.model_version.tolist()

        return ontology

    def set_latest(
            self,
            model_type: str,
            organisation: str,
            model_topology: str
    ):
        """
        Set model ID to latest model in given ontology group.

        :param model_type: Identifier of model_type to select.
        :param organisation: Identifier of organisation to select.
        :param model_topology: Identifier of model_topology to select
        :return:
        """
        assert model_type in self.ontology.keys(), "model_type requested was not found in ontology"
        assert organisation in self.ontology[model_type].keys(), \
            "organisation requested was not found in ontology"
        assert model_topology in self.ontology[model_type][organisation].keys(), \
            "model_topology requested was not found in ontology"

        versions = self.versions(
            model_type=model_type,
            organisation=organisation,
            model_topology=model_topology
        )
        self.model_type = model_type
        self.organisation = organisation
        self.model_topology = model_topology  # set to model for now, could be organism/organ specific later

        self.model_version = self._order_versions(versions=versions)[0]
        self.model_id = '_'.join([
            'embedding',
            self.id,
            self.model_type,
            self.organisation,
            self.model_topology,
            self.model_version
        ])


class ModelZooCelltype(ModelZoo):
    """
    Note on topology id: The topology ID is x.y.z, x is the major cell type version and y.z is the cell type model
    topology. Cell type model ontologies do not include the output size as this is set by the cell type version.
    """

    def load_ontology_from_model_ids(
            self,
            model_ids
    ) -> dict:
        """
        Load model ontology based on models available in model lookup tables.

        :param model_ids: Table listing all available model_ids.
        :return: Dictionary formatted ontology.
        """

        ids = [i for i in model_ids if i.split('_')[0] == 'celltype']
        id_df = pd.DataFrame(
            [i.split('_')[1:6] for i in ids],
            columns=['id', 'model_type', 'organisation', 'model_topology', 'model_version']
        )
        model = np.unique(id_df['model_type'])
        ontology = dict.fromkeys(model)
        for m in model:
            id_df_m = id_df[id_df.model_type == m]
            orga = np.unique(id_df_m['organisation'])
            ontology[m] = dict.fromkeys(orga)
            for org in orga:
                id_df_org = id_df_m[id_df_m.organisation == org]
                topo = np.unique(id_df_org['model_topology'])
                ontology[m][org] = dict.fromkeys(topo)
                for t in topo:
                    id_df_t = id_df_org[id_df_org.model_topology == t]
                    ontology[m][org][t] = id_df_t.model_version.tolist()

        return ontology

    def set_latest(
            self,
            model_type: str,
            organisation: str,
            model_topology: str
    ):
        """
        Set model ID to latest model in given ontology group.

        :param organism: Identifier of organism to select.
        :param organ: Identifier of organ to select.
        :param model_type: Identifier of model_type to select.
        :param organisation: Identifier of organisation to select.
        :param model_topology: Identifier of model_topology to select
        :return:
        """
        assert model_type in self.ontology.keys(), "model_type requested was not found in ontology"
        assert organisation in self.ontology[model_type].keys(), \
            "organisation requested was not found in ontology"
        assert model_topology in self.ontology[model_type][organisation].keys(), \
            "model_topology requested was not found in ontology"

        versions = self.versions(
            model_type=model_type,
            organisation=organisation,
            model_topology=model_topology
        )

        self.model_type = model_type
        self.organisation = organisation
        self.model_topology = model_topology  # set to model for now, could be organism/organ specific later

        self.model_version = self._order_versions(versions=versions)[0]
        self.model_id = '_'.join([
            'celltype',
            self.id,
            self.model_type,
            self.organisation,
            self.model_topology,
            self.model_version
        ])
