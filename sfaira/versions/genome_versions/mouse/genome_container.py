import os
import pandas

from .genome_sizes import GENOME_SIZE_DICT


class GenomeContainer:
    available_genomes = ["Mus_musculus_GRCm38_97"]

    def __init__(self):
        self.genomes = {
            "Mus_musculus_GRCm38_97": "Mus_musculus_GRCm38_97.csv"
        }
        self.genome_sizes = {
            "Mus_musculus_GRCm38_97": GENOME_SIZE_DICT["Mus_musculus_GRCm38_97"]
        }

    def read_local_csv(self, genome):
        return pandas.read_csv(os.path.join(str(os.path.dirname(__file__)), self.genomes[genome]))