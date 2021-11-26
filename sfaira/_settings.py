"""
Settings class which for example holds paths to cache directories used throughout the code.
"""

import os


SFAIRA_REPO_URL = "https://zenodo.org/record/4836517/files/"


class SfairaConfig:
    """\
    Config manager for sfaira.
    """

    def __init__(self):
        self.sfaira_repo_url = SFAIRA_REPO_URL
        self._cachedir_base = os.path.join(os.path.expanduser("~"), ".cache", "sfaira")
        self._cachedir_databases = os.path.join(self._cachedir_base, "dataset_meta")
        self._cachedir_databases_cellxgene = os.path.join(self._cachedir_databases, "cellxgene")
        self._cachedir_genomes = os.path.join(self._cachedir_base, "genomes")
        self._cachedir_ontologies = os.path.join(self._cachedir_base, "ontologies")

    @property
    def cachedir_base(self) -> str:
        os.makedirs(self._cachedir_base, exist_ok=True)
        return self._cachedir_base

    @cachedir_base.setter
    def cachedir_base(self, cachedir_base):
        if not isinstance(cachedir_base, str):
            raise ValueError(f"cachedir_base needs to be provided as a string, was {type(cachedir_base)}")
        if cachedir_base == "repo":
            cachedir_base = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "cache")
        self._cachedir_base = cachedir_base

    @property
    def cachedir_databases(self) -> str:
        os.makedirs(self._cachedir_databases, exist_ok=True)
        return self._cachedir_databases

    @cachedir_databases.setter
    def cachedir_databases(self, cachedir_databases):
        raise ValueError("cachedir_databases cannot be set manually as it is defined as a subdirectory of"
                         " cachedir_base. please modify cachedir_base instead")

    @property
    def cachedir_databases_cellxgene(self) -> str:
        os.makedirs(self._cachedir_databases_cellxgene, exist_ok=True)
        return self._cachedir_databases_cellxgene

    @cachedir_databases_cellxgene.setter
    def cachedir_databases_cellxgene(self, cachedir_databases_cellxgene):
        raise ValueError("cachedir_databases_cellxgene cannot be set manually as it is defined as a subdirectory"
                         " of cachedir_base. please modify cachedir_base instead")

    @property
    def cachedir_genomes(self) -> str:
        os.makedirs(self._cachedir_genomes, exist_ok=True)
        return self._cachedir_genomes

    @cachedir_genomes.setter
    def cachedir_genomes(self, cachedir_genomes):
        raise ValueError("cachedir_genomes cannot be set manually as it is defined as a subdirectory of cachedir_base."
                         "please modify cachedir_base instead")

    @property
    def cachedir_ontologies(self) -> str:
        os.makedirs(self._cachedir_ontologies, exist_ok=True)
        return self._cachedir_ontologies

    @cachedir_ontologies.setter
    def cachedir_ontologies(self, cachedir_ontologies):
        raise ValueError("cachedir_ontologies cannot be set manually as it is defined as a subdirectory of cachedir_base. please modify cachedir_base instead")


settings = SfairaConfig()
