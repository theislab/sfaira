from typing import Union

from sfaira.versions.genomes.genomes import GenomeContainer, ReactiveFeatureContainer


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
            if self.topology["input"]["genome"] is not None:
                self.gc = GenomeContainer(
                    organism=" ".join(self.topology["input"]["genome"].split(".")[0].split("_")),
                    release=self.topology["input"]["genome"].split(".")[-1],
                )
            else:
                self.gc = ReactiveFeatureContainer()
        else:
            assert isinstance(custom_genome_constainer, GenomeContainer)
            self.gc = custom_genome_constainer
        self.gc.set(**dict([tuple(self.topology["input"]["genes"])]))
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
