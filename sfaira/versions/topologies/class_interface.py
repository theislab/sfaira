from sfaira.versions.genomes import GenomeContainer


class TopologyContainer:

    """
    Class interface for a YAML-style defined model topology that loads a genome container tailored to the model.
    """

    def __init__(
            self,
            topology: dict,
            topology_id: str,
    ):
        self.topology = topology
        self.gc = GenomeContainer(assembly=self.topology["input"]["genome"])
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
