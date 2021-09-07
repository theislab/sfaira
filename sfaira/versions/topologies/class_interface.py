from typing import Union

from sfaira.versions.genomes.genomes import GenomeContainer


class TopologyContainer:

    """
    Class interface for a YAML-style defined model topology that loads a genome container tailored to the model.
    """

    def __init__(
            self,
            topology: dict,
            topology_id: str,
            custom_genome_constainer: Union[GenomeContainer, None] = None,
    ):
        self.topology = topology
        if custom_genome_constainer is None:
            self.gc = GenomeContainer(assembly=self.topology["input"]["genome"])
        else:
            assert isinstance(custom_genome_constainer, GenomeContainer)
            self.gc = custom_genome_constainer
        self.gc.subset(**dict([tuple(self.topology["input"]["genes"])]))
        self.topology_id = topology_id

    @property
    def model_type(self):
        return self.topology["model_type"]

    @property
    def output(self):
        return self.topology["output"]

    @property
    def n_var(self):
        return self.gc.n_var

    @property
    def organism(self):
        return self.gc.organism
