import abc
import pandas


class SuperGenomeContainer:
    _cache_tab: pandas.DataFrame
    genome: str
    species: str

    def __init__(
            self,
            species: str,
            genome: str
    ):
        self.species = species
        if self.species == "human":
            try:
                from sfaira_extension.versions.genome_versions.human import GenomeContainer
                if genome not in GenomeContainer.available_genomes:
                    from .human import GenomeContainer
                    if genome not in GenomeContainer.available_genomes:
                        raise ValueError(f"Genome {genome} not recognised.")
            except ImportError:
                from .human import GenomeContainer
                if genome not in GenomeContainer.available_genomes:
                    raise ValueError(f"Genome {genome} not recognised.")
        elif self.species == "mouse":
            try:
                from sfaira_extension.versions.genome_versions.mouse import GenomeContainer
                if genome not in GenomeContainer.available_genomes:
                    from .mouse import GenomeContainer
                    if genome not in GenomeContainer.available_genomes:
                        raise ValueError(f"Genome {genome} not recognised.")
            except ImportError:
                from .mouse import GenomeContainer
                if genome not in GenomeContainer.available_genomes:
                    raise ValueError(f"Genome {genome} not recognised.")
        else:
            raise ValueError(f"Species {species} not recognised.")

        self.gc = GenomeContainer()
        self.set_genome(genome=genome)

    @property
    def cache_tab(self):
        return self._cache_tab

    def set_genome(self, genome):
        self.genome = genome
        self._cache_tab = self.gc.read_local_csv(genome=genome)
        assert self.gc.genome_sizes[self.genome][0] == self.cache_tab.shape[0]

    def show_genomes(self):
        return list(self.gc.genomes.keys())

    @property
    def names(self):
        return self.cache_tab["name"].values.tolist()

    @property
    def ensembl(self):
        return self.cache_tab["ensg"].values.tolist()

    @property
    def type(self):
        return self.cache_tab["type"].values.tolist()

    @property
    def ngenes(self):
        return self.gc.genome_sizes[self.genome][0]

    @property
    def names_to_id_dict(self):
        return dict(zip(self.cache_tab["name"].values.tolist(), self.cache_tab["ensg"].values.tolist()))

    @property
    def id_to_names_dict(self):
        return dict(zip(self.cache_tab["ensg"].values.tolist(), self.cache_tab["name"].values.tolist()))

    @property
    def strippednames_to_id_dict(self):
        return dict(zip([i.split(".")[0] for i in self.cache_tab["name"]], self.cache_tab["ensg"].values.tolist()))
