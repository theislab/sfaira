"""
Functionalities to interact with gene sets defined in an assembly and gene-annotation (such as protein-coding).
"""

import gzip
import numpy as np
import os
from typing import Union
import pandas
import pathlib
import urllib.error
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
        cache_dir_path = pathlib.Path(cache_dir)
        cache_dir_path.mkdir(parents=True, exist_ok=True)
        return cache_dir

    @property
    def cache_fn(self):
        return os.path.join(self.cache_dir, self.assembly + ".csv")

    @property
    def release(self) -> str:
        return self.assembly.split(".")[-1]

    @property
    def organism(self) -> str:
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
        try:
            _ = urllib.request.urlretrieve(url=self.url_ensembl_ftp, filename=temp_file)
        except urllib.error.URLError as e:
            raise ValueError(f"Could not download gtf from {self.url_ensembl_ftp} with urllib.error.URLError: {e}")
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

    def __init__(
            self,
            organism: Union[None, str] = None,
            assembly: Union[None, str] = None,
    ):
        if assembly is None:
            # Set defaults based on organism if assembly is not given.
            if organism is None:
                raise ValueError("Supply either organism or assembly to GenomeContainer().")
            if organism == "human":
                self.assembly = "Homo_sapiens.GRCh38.102"
            elif organism == "mouse":
                self.assembly = "Mus_musculus.GRCm38.102"
            else:
                raise ValueError(f"organism {organism} not found")
        else:
            self.assembly = assembly
        self.gtfi = GtfInterface(assembly=self.assembly)
        self.load_genome()

    @property
    def organism(self):
        return self.gtfi.organism

    def load_genome(self):
        self.genome_tab = self.gtfi.cache

    def subset(
            self,
            biotype: Union[None, str] = None,
            symbols: Union[None, str] = None,
            ensg: Union[None, str] = None,
    ):
        """
        Subset by gene biotype or to gene list defined by identifiers (symbol or ensemble ID).

        Can subset by multiple factors at the same time.

        :param biotype: Gene biotype(s) of gene(s) to subset genome to. Separate in string via "," if choosing multiple.
        :param symbols: Gene symbol(s) of gene(s) to subset genome to. Separate in string via "," if choosing multiple.
        :param ensg: Ensemble gene ID(s) of gene(s) to subset genome to. Separate in string via "," if choosing
            multiple.
        """
        subset = np.ones((self.n_var,), "int") == 1
        if biotype is not None:
            if not isinstance(biotype, str):
                raise ValueError("Supply biotype as string, see also function annotation.")
            biotype = biotype.split(",")
            subset = np.logical_and(
                subset,
                [x in biotype for x in self.genome_tab[KEY_TYPE].values]
            )
        if symbols is not None:
            if not isinstance(symbols, str):
                raise ValueError("Supply symbols as string, see also function annotation.")
            symbols = symbols.split(",")
            subset = np.logical_and(
                subset,
                [x in symbols for x in self.genome_tab[KEY_SYMBOL].values]
            )
        if ensg is not None:
            if not isinstance(ensg, str):
                raise ValueError("Supply ensg as string, see also function annotation.")
            ensg = ensg.split(",")
            subset = np.logical_and(
                subset,
                [x in ensg for x in self.genome_tab[KEY_ID].values]
            )
        self.genome_tab = self.genome_tab.loc[subset, :].copy()

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
    def n_var(self) -> int:
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
