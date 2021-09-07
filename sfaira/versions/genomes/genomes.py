"""
Functionalities to interact with gene sets defined in an assembly and gene-annotation (such as protein-coding).
"""

import gzip
import numpy as np
import os
from typing import Iterable, List, Union
import pandas
import pathlib
import urllib.error
import urllib.request

from sfaira.consts.directories import CACHE_DIR_GENOMES

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
        cache_dir_path = pathlib.Path(CACHE_DIR_GENOMES)
        cache_dir_path.mkdir(parents=True, exist_ok=True)
        return CACHE_DIR_GENOMES

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
        try:
            _ = urllib.request.urlretrieve(url=self.url_ensembl_ftp, filename=temp_file)
        except urllib.error.URLError as e:
            raise ValueError(f"Could not download gtf from {self.url_ensembl_ftp} with urllib.error.URLError: {e}, "
                             f"check if assembly name '{self.assembly}' corresponds to an actual assembly.")
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
    """
    Container class for a genome annotation for a specific release.

    This class can be used to translate between symbols and ENSEMBL IDs for a specific assembly, to store specific gene
    subsets of an assembly, and to subselect genes by biotypes in an assembly.
    """

    genome_tab: pandas.DataFrame
    assembly: str

    def __init__(
            self,
            assembly: str = None,
    ):
        """
        Are you not sure which assembly to use?

            - You could use the newest one for example, check the ENSEMBL site regularly for updates:
                http://ftp.ensembl.org/pub/
            - You could use one used by a specific aligner, the assemblies used by 10x cellranger are described here
                for example: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build

        :param assembly: The full name of the genome assembly, e.g. Homo_sapiens.GRCh38.102.
        """
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
    def symbols(self) -> List[str]:
        """
        List of symbols of genes in genome container.
        """
        return self.genome_tab[KEY_SYMBOL].values.tolist()

    @property
    def ensembl(self) -> List[str]:
        """
        List of ENSEMBL IDs of genes in genome container.
        """
        return self.genome_tab[KEY_ID].values.tolist()

    @property
    def biotype(self) -> List[str]:
        """
        List of biotypes of genes in genome container.
        """
        return self.genome_tab[KEY_TYPE].values.tolist()

    def __validate_ensembl(self, x: List[str], enforce_captitalization: bool = True):
        if enforce_captitalization:
            not_found = [y for y in x if y not in self.ensembl]
        else:
            ensembl_upper = [y.upper() for y in self.ensembl]
            not_found = [y for y in x if y.upper() not in ensembl_upper]
        if len(not_found) > 0:
            raise ValueError(f"Could not find ENSEMBL ID: {not_found}")

    def __validate_symbols(self, x: List[str], enforce_captitalization: bool = True):
        if enforce_captitalization:
            not_found = [y for y in x if y not in self.symbols]
        else:
            symbols_upper = [y.upper() for y in self.symbols]
            not_found = [y for y in x if y.upper() not in symbols_upper]
        if len(not_found) > 0:
            raise ValueError(f"Could not find symbol: {not_found}")

    def __validate_types(self, x: List[str]):
        not_found = [y for y in x if y not in self.biotype]
        if len(not_found) > 0:
            raise ValueError(f"Could not find type: {not_found}")

    @property
    def n_var(self) -> int:
        """
        Number of genes in genome container.
        """
        return self.genome_tab.shape[0]

    @property
    def symbol_to_id_dict(self):
        """
        Dictionary-formatted map of gene symbols to ENSEMBL IDs.
        """
        return dict(zip(self.genome_tab[KEY_SYMBOL].values.tolist(), self.genome_tab[KEY_ID].values.tolist()))

    @property
    def id_to_symbols_dict(self):
        """
        Dictionary-formatted map of ENSEMBL IDs to gene symbols.
        """
        return dict(zip(self.genome_tab[KEY_ID].values.tolist(), self.genome_tab[KEY_SYMBOL].values.tolist()))

    def translate_symbols_to_id(self, x: Union[str, Iterable[str]]) -> Union[str, List[str]]:
        """
        Translate gene symbols to ENSEMBL IDs.

        :param x: Symbol(s) to translate.
        :return: ENSEMBL IDs
        """
        if isinstance(x, str):
            x = [x]
        self.__validate_symbols(x=x, enforce_captitalization=False)
        map_dict = dict([(k.upper(), v) for k, v in self.symbol_to_id_dict.items()])
        y = [map_dict[xx.upper()] for xx in x]
        if len(y) == 1:
            y = y[0]
        return y

    def translate_id_to_symbols(self, x: Union[str, Iterable[str]]) -> Union[str, List[str]]:
        """
        Translate ENSEMBL IDs to gene symbols.

        :param x: ENSEMBL ID(s) to translate.
        :return: Gene symbols.
        """
        if isinstance(x, str):
            x = [x]
        self.__validate_ensembl(x=x, enforce_captitalization=False)
        map_dict = dict([(k.upper(), v) for k, v in self.id_to_symbols_dict.items()])
        y = [map_dict[xx.upper()] for xx in x]
        if len(y) == 1:
            y = y[0]
        return y

    @property
    def strippednames_to_id_dict(self):
        return dict(zip([i.split(".")[0] for i in self.genome_tab[KEY_SYMBOL]],
                        self.genome_tab[KEY_ID].values.tolist()))


class CustomFeatureContainer(GenomeContainer):

    def __init__(
            self,
            genome_tab: pandas.DataFrame,
            organism: str,
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
        self._organism = organism

    @property
    def organism(self):
        return self._organism
