from typing import Iterable, List, Union

from sfaira.versions.genomes import GenomeContainer


def translate_symbols_to_id(x: Union[str, Iterable[str]], assembly: str) -> Union[str, List[str]]:
    """
    Translate gene symbols to ENSEMBL IDs.

    Input captitalization is ignored but the output capitalisation matches the ENSEMBL .gtf files.

    Are you not sure which assembly to use?

    - You could use the newest one for example, check the ENSEMBL site regularly for updates:
        http://ftp.ensembl.org/pub/
    - You could use one used by a specific aligner, the assemblies used by 10x cellranger are described here
        for example: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build

    :param x: Symbol(s) to translate.
    :param assembly: The full name of the genome assembly, e.g. "Homo_sapiens.GRCh38.102".
    :return: ENSEMBL IDs
    """
    return GenomeContainer(assembly=assembly).translate_symbols_to_id(x=x)


def translate_id_to_symbols(x: Union[str, Iterable[str]], assembly: str) -> Union[str, List[str]]:
    """
    Translate ENSEMBL IDs to gene symbols.

    Input captitalization is ignored but the output capitalisation matches the ENSEMBL .gtf files.

    Are you not sure which assembly to use?

    - You could use the newest one for example, check the ENSEMBL site regularly for updates:
        http://ftp.ensembl.org/pub/
    - You could use one used by a specific aligner, the assemblies used by 10x cellranger are described here
        for example: https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build

    :param x: ENSEMBL ID(s) to translate.
    :param assembly: The full name of the genome assembly, e.g. "Homo_sapiens.GRCh38.102".
    :return: Gene symbols.
    """
    return GenomeContainer(assembly=assembly).translate_id_to_symbols(x=x)
