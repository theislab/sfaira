import abc
try:
    import kipoi
except ImportError:
    kipoi = None
import numpy as np
import pandas as pd
from typing import List, Union

from .external import celltype_versions, Topologies


class ModelZoo(abc.ABC):
    """
    Model ontology base class.
    """
    topology_container: Topologies
    ontology: dict
    model_id: Union[str, None]
    model_class: Union[str, None]
    species: Union[str, None]
    organ: Union[str, None]
    model_class: Union[str, None]
    model_type: Union[str, None]
    model_topology: Union[str, None]
    model_version: Union[str, None]
    celltypes: Union[List, None]

    def __init__(
            self,
            model_lookuptable: Union[None, pd.DataFrame] = None
    ):
        """
        :param model_lookuptable: model_lookuptable.
        """
        if model_lookuptable is not None:  # check if models in repository
            self.ontology = self.load_ontology_from_model_ids(model_lookuptable['model_id'].values)
        self.model_id = None
        self.model_class = None
        self.species = None
        self.organ = None
        self.model_type = None
        self.organisation = None
        self.model_topology = None
        self.model_version = None
        self.topology_container = None
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
        :return:
        """
        self.model_id = model_id
        ixs = self.model_id.split('_')
        self.model_class = ixs[0]
        self.species = ixs[1]
        self.organ = ixs[2]
        self.model_type = ixs[3]
        self.organisation = ixs[4]
        self.model_topology = ixs[5]
        self.model_version = ixs[6]

        self.topology_container = Topologies(
            species=self.species,
            model_class=self.model_class,
            model_type=self.model_type,
            topology_id=self.model_topology
        )

    def save_weights_to_remote(self, path=None):
        """
        Saves model weights to repository XY.
        Increments 3rd digit of version number.
        Adds model_id to the text file, updates model_index

        :return:
        """
        raise NotImplementedError()

    def save_weights_to_public(self):
        """
        Saves model weights to cloud under an organization name.
        Increments 2nd digit of version number.
        Adds model_id to the text file, updates model_index

        :return:
        """
        raise NotImplementedError()

    def call_kipoi(self):
        """
        Returns kipoi_experimental model call from remote directly on local data using kipoi_experimental.

        Runs model defined in self.model_id.
        For this, the remote server associated with the model_id has to be identified via find_remote().

        :return: Predictions
        """
        return kipoi.get_model(
            self.model_id,
            source='kipoi_experimental',
            with_dataloader=True
        )  # TODO make sure that this is in line with kipoi_experimental model names
        # alternatively:
        #return kipoi_experimental.get_model("https://github.com/kipoi/models/tree/7d3ea7800184de414aac16811deba6c8eefef2b6/pwm_HOCOMOCO/human/CTCF", source='github-permalink')

    def species(self) -> List[str]:
        """
        Return list of available species.

        :return: List of species available.
        """
        return self.ontology.keys()

    def organs(
            self,
            species: str
    ) -> List[str]:
        """
        Return list of available organs for a given species.

        :param species: Identifier of species to show organs for.
        :return: List of organs available.
        """
        assert species in self.ontology.keys(), "species requested was not found in ontology"
        return self.ontology[species].keys()

    def models(
            self,
            species: str,
            organ: str
    ) -> List[str]:
        """
        Return list of available models for a given species, organ.

        :param species: Identifier of species to show organs for.
        :param organ: Identifier of organ to show versions for.
        :return: List of models available.
        """
        assert species in self.ontology.keys(), "species requested was not found in ontology"
        assert organ in self.ontology[species].keys(), "organ requested was not found in ontology"
        return self.ontology[species][organ].keys()

    def organisation(
            self,
            species: str,
            organ: str,
            model_type: str
    ) -> List[str]:
        """
        Return list of available organisation that trained a given model for a given species and organ

        :param species: Identifier of species to show versions for.
        :param organ: Identifier of organ to show versions for.
        :param model_type: Identifier of model to show versions for.
        :return: List of versions available.
        """
        assert species in self.ontology.keys(), "species requested was not found in ontology"
        assert organ in self.ontology[species].keys(), "organ requested was not found in ontology"
        assert model_type in self.ontology[species][organ].keys(), "model_type requested was not found in ontology"
        return self.ontology[species][organ][model_type]

    def topology(
            self,
            species: str,
            organ: str,
            model_type: str,
            organisation: str
    ) -> List[str]:
        """
        Return list of available model topologies that trained by a given organisation,
        a given model for a given species and organ

        :param species: Identifier of species to show versions for.
        :param organ: Identifier of organ to show versions for.
        :param model_type: Identifier of model_type to show versions for.
        :param organisation: Identifier of organisation to show versions for.
        :return: List of versions available.
        """
        assert species in self.ontology.keys(), "species requested was not found in ontology"
        assert organ in self.ontology[species].keys(), "organ requested was not found in ontology"
        assert model_type in self.ontology[species][organ].keys(), "model_type requested was not found in ontology"
        assert organisation in self.ontology[species][organ][model_type].keys(), \
            "organisation requested was not found in ontology"
        return self.ontology[species][organ][model_type][organisation]

    def versions(
            self,
            species: str,
            organ: str,
            model_type: str,
            organisation: str,
            model_topology: str
    ) -> List[str]:
        """
        Return list of available model versions of a given organisation for a given species and organ and model.

        :param species: Identifier of species to show versions for.
        :param organ: Identifier of organ to show versions for.
        :param model_type: Identifier of model_type to show versions for.
        :param organisation: Identifier of organisation to show versions for.
        :param model_topology: Identifier of model_topology to show versions for.
        :return: List of versions available.
        """
        assert species in self.ontology.keys(), "species requested was not found in ontology"
        assert organ in self.ontology[species].keys(), "organ requested was not found in ontology"
        assert model_type in self.ontology[species][organ].keys(), "model_type requested was not found in ontology"
        assert organisation in self.ontology[species][organ][model_type].keys(), \
            "organisation requested was not found in ontology"
        assert model_topology in self.ontology[species][organ][model_type][organisation].keys(), \
            "model_topology requested was not found in ontology"
        return self.ontology[species][organ][model_type][organisation][model_topology]

    @property
    def genome(self):
        return self.model_hyperparameters["genome"]

    @property
    def gene_names(self):
        return self.topology_container.genome_container.names

    @property
    def ensemble_names(self):
        return self.topology_container.genome_container.ensembl

    @property
    def model_hyperparameters(self) -> dict:
        assert self.topology_container is not None
        return self.topology_container.topology["hyper_parameters"]


class ModelZooEmbedding(ModelZoo):
    """
    The supported model ontology is:

        species -> organ -> model -> organisation -> topology -> version -> ID

    Maybe: include experimental protocol? Ie droplet, full-length, single-nuclei.
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

        ids = [i for i in model_ids if i.split('_')[0] == 'embedding']
        id_df = pd.DataFrame(
            [i.split('_')[1:7] for i in ids],
            columns=['species', 'organ', 'model_type', 'organisation', 'model_topology', 'model_version']
        )
        species = np.unique(id_df['species'])
        ontology = dict.fromkeys(species)
        for g in species:
            id_df_g = id_df[id_df.species == g]
            organ = np.unique(id_df_g['organ'])
            ontology[g] = dict.fromkeys(organ)
            for o in organ:
                id_df_o = id_df_g[id_df_g.organ == o]
                model = np.unique(id_df_o['model_type'])
                ontology[g][o] = dict.fromkeys(model)
                for m in model:
                    id_df_m = id_df_o[id_df_o.model_type == m]
                    orga = np.unique(id_df_m['organisation'])
                    ontology[g][o][m] = dict.fromkeys(orga)
                    for org in orga:
                        id_df_org = id_df_m[id_df_m.organisation == org]
                        topo = np.unique(id_df_org['model_topology'])
                        ontology[g][o][m][org] = dict.fromkeys(topo)
                        for t in topo:
                            id_df_t = id_df_org[id_df_org.model_topology == t]
                            ontology[g][o][m][org][t] = id_df_t.model_version.tolist()

        return ontology

    def set_latest(
            self,
            species: str,
            organ: str,
            model_type: str,
            organisation: str,
            model_topology: str
    ):
        """
        Set model ID to latest model in given ontology group.

        :param species: Identifier of species to select.
        :param organ: Identifier of organ to select.
        :param model_type: Identifier of model_type to select.
        :param organisation: Identifier of organisation to select.
        :param model_topology: Identifier of model_topology to select
        :return:
        """
        assert species in self.ontology.keys(), "species requested was not found in ontology"
        assert organ in self.ontology[species].keys(), "organ requested was not found in ontology"
        assert model_type in self.ontology[species][organ].keys(), "model_type requested was not found in ontology"
        assert organisation in self.ontology[species][organ][model_type].keys(), \
            "organisation requested was not found in ontology"
        assert model_topology in self.ontology[species][organ][model_type][organisation].keys(), \
            "model_topology requested was not found in ontology"

        versions = self.versions(
            species=species,
            organ=organ,
            model_type=model_type,
            organisation=organisation,
            model_topology=model_topology
        )
        self.species = species
        self.organ = organ
        self.model_type = model_type
        self.organisation = organisation
        self.model_topology = model_topology  # set to model for now, could be species/organ specific later

        self.model_version = self._order_versions(versions=versions)[0]
        self.model_id = '_'.join([
            'embedding',
            self.species,
            self.organ,
            self.model_type,
            self.organisation,
            self.model_topology,
            self.model_version
        ])
        self.topology_container = Topologies(
            species=self.species,
            model_class="embedding",
            model_type=self.model_type,
            topology_id=self.model_topology
        )


class ModelZooCelltype(ModelZoo):
    """
    The supported model ontology is:

        species -> organ -> model -> organisation -> topology -> version -> ID

    Maybe: include experimental protocol? Ie droplet, full-length, single-nuclei.

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
            [i.split('_')[1:7] for i in ids],
            columns=['species', 'organ', 'model_type', 'organisation', 'model_topology', 'model_version']
        )
        species = np.unique(id_df['species'])
        ontology = dict.fromkeys(species)
        for g in species:
            id_df_g = id_df[id_df.species == g]
            organ = np.unique(id_df_g['organ'])
            ontology[g] = dict.fromkeys(organ)
            for o in organ:
                id_df_o = id_df_g[id_df_g.organ == o]
                model = np.unique(id_df_o['model_type'])
                ontology[g][o] = dict.fromkeys(model)
                for m in model:
                    id_df_m = id_df_o[id_df_o.model_type == m]
                    orga = np.unique(id_df_m['organisation'])
                    ontology[g][o][m] = dict.fromkeys(orga)
                    for org in orga:
                        id_df_org = id_df_m[id_df_m.organisation == org]
                        topo = np.unique(id_df_org['model_topology'])
                        ontology[g][o][m][org] = dict.fromkeys(topo)
                        for t in topo:
                            id_df_t = id_df_org[id_df_org.model_topology == t]
                            ontology[g][o][m][org][t] = id_df_t.model_version.tolist()

        return ontology

    def set_latest(
            self,
            species: str,
            organ: str,
            model_type: str,
            organisation: str,
            model_topology: str
    ):
        """
        Set model ID to latest model in given ontology group.

        :param species: Identifier of species to select.
        :param organ: Identifier of organ to select.
        :param model_type: Identifier of model_type to select.
        :param organisation: Identifier of organisation to select.
        :param model_topology: Identifier of model_topology to select
        :return:
        """
        assert species in self.ontology.keys(), "species requested was not found in ontology"
        assert organ in self.ontology[species].keys(), "organ requested was not found in ontology"
        assert model_type in self.ontology[species][organ].keys(), "model_type requested was not found in ontology"
        assert organisation in self.ontology[species][organ][model_type].keys(), \
            "organisation requested was not found in ontology"
        assert model_topology in self.ontology[species][organ][model_type][organisation].keys(), \
            "model_topology requested was not found in ontology"

        versions = self.versions(
            species=species,
            organ=organ,
            model_type=model_type,
            organisation=organisation,
            model_topology=model_topology
        )

        self.species = species
        self.organ = organ
        self.model_type = model_type
        self.organisation = organisation
        self.model_topology = model_topology  # set to model for now, could be species/organ specific later

        self.model_version = self._order_versions(versions=versions)[0]
        self.model_id = '_'.join([
            'celltype',
            self.species,
            self.organ,
            self.model_type,
            self.organisation,
            self.model_topology,
            self.model_version
        ])
        self.topology_container = Topologies(
            species=self.species,
            model_class="celltype",
            model_type=self.model_type,
            topology_id=self.model_topology
        )
        self.celltypes = celltype_versions.SPECIES_DICT[self.species][self.organ].celltype_universe[self.model_version.split(".")[0]]
