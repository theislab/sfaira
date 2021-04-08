"""
Functionalities to interact with gene sets defined in an assembly and gene-annotation (such as protein-coding).
"""

import gzip
import os
from typing import Union
import pandas
import urllib.request

KEY_SYMBOL = "gene_name"
KEY_ID = "gene_id"
KEY_TYPE = "gene_biotype"
VALUE_GTF_GENE = "gene"
KEY_GTF_REGION_TYPE = 2
KEY_GTF_REGION_DETAIL_FIELD = 8
IDX_GTF_REGION_DETAIL_FIELD_ID = 0
IDX_GTF_REGION_DETAIL_FIELD_SYMBOL = 2
IDX_GTF_REGION_DETAIL_FIELD_TYPE = 4


class GtfInterface:

    def __init__(self, assembly: str):
        self.assembly = assembly

    @property
    def cache_dir(self):
        """
        The cache dir is in a cache directory in the sfaira installation that is excempt from git versioning.
        """
        cache_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "cache", "genomes")
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)
        return cache_dir

    @property
    def cache_fn(self):
        return os.path.join(self.cache_dir, self.assembly + ".csv")

    @property
    def release(self) -> str:
        return self.assembly.split(".")[-1]

    @property
    def organism(self):
        return self.assembly.split(".")[0].lower()

    @property
    def url_ensembl_ftp(self):
        return f"ftp://ftp.ensembl.org/pub/release-{self.release}/gtf/{self.organism}/{self.assembly}.gtf.gz"

    def download_gtf_ensembl(self):
        """
        Download .gtf file from ensembl FTP server and turn into reduced, gene-centric cache .csv.
        """
        temp_file = os.path.join(self.cache_dir, self.assembly + ".gtf.gz")
        print(f"downloading {self.url_ensembl_ftp} into a temporary file {temp_file}")
        _ = urllib.request.urlretrieve(url=self.url_ensembl_ftp, filename=temp_file)
        with gzip.open(temp_file) as f:
            tab = pandas.read_csv(f, sep="\t", comment="#", header=None)
        os.remove(temp_file)  # Delete temporary file .gtf.gz.
        tab = tab.loc[tab[KEY_GTF_REGION_TYPE].values == VALUE_GTF_GENE, :]
        conversion_tab = pandas.DataFrame({
            "gene_id": [
                x.split(";")[IDX_GTF_REGION_DETAIL_FIELD_ID].split(" ")[-1].strip("\"")
                for x in tab[KEY_GTF_REGION_DETAIL_FIELD].values],
            "gene_name": [
                x.split(";")[IDX_GTF_REGION_DETAIL_FIELD_SYMBOL].split(" ")[-1].strip("\"")
                for x in tab[KEY_GTF_REGION_DETAIL_FIELD].values],
            "gene_biotype": [
                x.split(";")[IDX_GTF_REGION_DETAIL_FIELD_TYPE].split(" ")[-1].strip("\"")
                for x in tab[KEY_GTF_REGION_DETAIL_FIELD].values],
        }).sort_values("gene_id")
        conversion_tab.to_csv(self.cache_fn)

    @property
    def cache(self) -> pandas.DataFrame:
        if not os.path.exists(self.cache_fn):
            self.download_gtf_ensembl()
        return pandas.read_csv(self.cache_fn)


class GenomeContainer:
    genome_tab: pandas.DataFrame
    assembly: str
    organism: str

    def __init__(
            self,
            organism: str,
            assembly: Union[None, str],
    ):
        self.organism = organism
        # Set defaults:
        if self.organism == "human":
            self.assembly = assembly if assembly is not None else "Homo_sapiens.GRCh38.102"
        elif self.organism == "mouse":
            self.assembly = assembly if assembly is not None else "Mus_musculus.GRCm38.102"
        else:
            raise ValueError(f"organism {organism} not found")
        self.gc = GtfInterface(assembly=self.assembly)
        self.load_genome()

    def load_genome(self):
        self.genome_tab = self.gc.cache

    def subset(self, gene_biotype: str):
        self.genome_tab = self.genome_tab.loc[self.genome_tab[KEY_TYPE].values == gene_biotype, :].copy()

    @property
    def names(self):
        return self.genome_tab[KEY_SYMBOL].values.tolist()

    @property
    def ensembl(self):
        return self.genome_tab[KEY_ID].values.tolist()

    @property
    def type(self):
        return self.genome_tab[KEY_TYPE].values.tolist()

    @property
    def ngenes(self) -> int:
        return self.genome_tab.shape[0]

    @property
    def names_to_id_dict(self):
        return dict(zip(self.genome_tab[KEY_SYMBOL].values.tolist(), self.genome_tab[KEY_ID].values.tolist()))

    @property
    def id_to_names_dict(self):
        return dict(zip(self.genome_tab[KEY_ID].values.tolist(), self.genome_tab[KEY_SYMBOL].values.tolist()))

    @property
    def strippednames_to_id_dict(self):
        return dict(zip([i.split(".")[0] for i in self.genome_tab[KEY_SYMBOL]],
                        self.genome_tab[KEY_ID].values.tolist()))
