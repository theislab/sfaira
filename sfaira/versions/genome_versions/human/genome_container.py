import os
import pandas

from .genome_sizes import GENOME_SIZE_DICT


class GenomeContainer:

    def __init__(self):
        self.genomes = {
            "Homo_sapiens_GRCh38_97": "Homo_sapiens_GRCh38_97.csv"
        }
        self.genome_sizes = {
            "Homo_sapiens_GRCh38_97": GENOME_SIZE_DICT["Homo_sapiens_GRCh38_97"]
        }

    def read_local_csv(self, genome):
        return pandas.read_csv(os.path.join(str(os.path.dirname(__file__)), self.genomes[genome]))