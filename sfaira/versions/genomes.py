"""
Functionalities to interact with gene sets defined in an assembly and gene-annotation (such as protein-coding).
"""

import gzip
import numpy as np
import os
from typing import List, Union
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
            KEY_ID: [
                x.split(";")[IDX_GTF_REGION_DETAIL_FIELD_ID].split(" ")[-1].strip("\"")
                for x in tab[KEY_GTF_REGION_DETAIL_FIELD].values],
            KEY_SYMBOL: [
                x.split(";")[IDX_GTF_REGION_DETAIL_FIELD_SYMBOL].split(" ")[-1].strip("\"")
                for x in tab[KEY_GTF_REGION_DETAIL_FIELD].values],
            KEY_TYPE: [
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
            assembly: str = None,
    ):
        if not isinstance(assembly, str):
            raise ValueError(f"supplied assembly {assembly} was not a string")
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
            biotype: Union[None, str, List[str]] = None,
            symbols: Union[None, str, List[str]] = None,
            ensg: Union[None, str, List[str]] = None,
    ):
        """
        Subset by gene biotype or to gene list defined by identifiers (symbol or ensemble ID).

        Will subset by multiple factors if more than one parameter is not None.

        :param biotype: Gene biotype(s) of gene(s) to subset genome to. Elements have to appear in genome.
            Separate in string via "," if choosing multiple or supply as list of string.
        :param symbols: Gene symbol(s) of gene(s) to subset genome to. Elements have to appear in genome.
            Separate in string via "," if choosing multiple or supply as list of string.
        :param ensg: Ensemble gene ID(s) of gene(s) to subset genome to. Elements have to appear in genome.
            Separate in string via "," if choosing multiple or supply as list of string.
        """
        subset = np.ones((self.n_var,), "int") == 1
        if biotype is not None:
            if isinstance(biotype, list):
                pass
            elif isinstance(biotype, str):
                biotype = biotype.split(",")
            else:
                raise ValueError(f"Supply biotype as string, see also function annotation. Supplied {biotype}.")
            self.__validate_types(x=biotype)
            subset = np.logical_and(
                subset,
                [x in biotype for x in self.genome_tab[KEY_TYPE].values]
            )
        if symbols is not None:
            if isinstance(symbols, list):
                pass
            elif isinstance(symbols, str):
                symbols = symbols.split(",")
            else:
                raise ValueError(f"Supply symbols as string, see also function annotation. Supplied {symbols}.")
            self.__validate_symbols(x=symbols)
            subset = np.logical_and(
                subset,
                [x in symbols for x in self.genome_tab[KEY_SYMBOL].values]
            )
        if ensg is not None:
            if isinstance(ensg, list):
                pass
            elif isinstance(ensg, str):
                ensg = ensg.split(",")
            else:
                raise ValueError(f"Supply ensg as string, see also function annotation. Supplied {ensg}.")
            self.__validate_ensembl(x=ensg)
            subset = np.logical_and(
                subset,
                [x in ensg for x in self.genome_tab[KEY_ID].values]
            )
        self.genome_tab = self.genome_tab.loc[subset, :].copy()

    @property
    def symbols(self):
        return self.genome_tab[KEY_SYMBOL].values.tolist()

    @property
    def ensembl(self):
        return self.genome_tab[KEY_ID].values.tolist()

    @property
    def biotype(self):
        return self.genome_tab[KEY_TYPE].values.tolist()

    def __validate_ensembl(self, x: List[str]):
        not_found = [y for y in x if y not in self.ensembl]
        if len(not_found) > 0:
            raise ValueError(f"Could not find ensembl: {not_found}")

    def __validate_symbols(self, x: List[str]):
        not_found = [y for y in x if y not in self.symbols]
        if len(not_found) > 0:
            raise ValueError(f"Could not find names: {not_found}")

    def __validate_types(self, x: List[str]):
        not_found = [y for y in x if y not in self.biotype]
        if len(not_found) > 0:
            raise ValueError(f"Could not find type: {not_found}")

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


class CustomFeatureContainer(GenomeContainer):

    def __init__(
            self,
            genome_tab: pandas.DataFrame,
    ):
        """

        :param genome_tab: Table characterising feature space. Must be a data frame with 3 columns:

            - "gene_name": Name of features.
            - "gene_id": ID of features, can be the same as values of "gene_name"
            - "gene_biotype": Types of features, can be arbitrary like "embedding"
        """
        self.assembly = "custom"
        assert len(genome_tab.columns) == 3
        assert KEY_SYMBOL in genome_tab.columns
        assert KEY_ID in genome_tab.columns
        assert KEY_TYPE in genome_tab.columns
        self.genome_tab = genome_tab
