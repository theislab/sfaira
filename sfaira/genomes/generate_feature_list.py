import numpy as np
import pandas
from typing import Union


class ExtractFeatureList:
    gene_table: Union[None, pandas.DataFrame]
    organism: Union[None, str]
    release: Union[None, str]

    def __init__(self):
        self.organism = None
        self.release = None
        self.gene_table = None

    def reduce_types(self, types: list):
        """
        Reduce entries in gene table to genes of a specific type.

        :param types: List of types to reduce to
        :return:
        """
        assert self.gene_table is not None
        self.gene_table = self.gene_table.loc[[x in types for x in self.gene_table["type"].values], :]


class ExtractFeatureListEnsemble(ExtractFeatureList):

    def from_ensemble_gtf(
            self,
            fn: str,
            skiprows: int = 5
    ):
        """
        Read ENSEMBL formatted .gtf and return HGNC and ENS gene names.

        :param fn: File name of .gtf to read.
        :param skiprows: Number of commen rows at top of file.
        :return:
        """
        gtf_name = fn.split("/")[-1]
        self.organism = gtf_name.split(".")[0]
        self.release = "_".join(gtf_name.split(".")[1:-1])

        tab = pandas.read_table(
            fn,
            sep="\t",
            skipinitialspace=True,
            skiprows=skiprows,
            header=None
        )
        tab_gene = tab.iloc[np.where([x == "gene" for x in tab.iloc[:, 2].values])[0], :]
        self.gene_table = pandas.DataFrame({
            "hgnc": [x.split(";")[2].split("gene_name ")[1].strip("\"").strip(" ")
                     for x in tab_gene.iloc[:, 8].values],
            "ensg": [x.split(";")[0].split("gene_id ")[1].strip("\"").strip(" ")
                     for x in tab_gene.iloc[:, 8].values],
            "type": [x.split(";")[4].split("gene_biotype ")[1].strip("\"").strip(" ")
                     for x in tab_gene.iloc[:, 8].values]
        }).sort_values(by="ensg")

    def reduce_types_protein_coding(self):
        self.reduce_types(types=["protein_coding"])

    def write_gene_table_to_csv(self, path):
        self.gene_table.to_csv(path_or_buf=path + self.organism + "_" + self.release + ".csv")
