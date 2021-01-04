import numpy as np


class CelltypeVersionsBase:
    """
    Versioned cell type universe (list) and ontology (hierarchy) container class.

    This class is subclassed once for each anatomical structure (organ).
    Cell type versions take the form x.y:
        - x is a major version and is incremented if the cell type identities or number changes.
        - y is a minor version and is incremented if cell type annotation such as names are altered without changing
            the described cell type or adding cell types.

    Basic checks on the organ specific instance are performed in the constructor.
    """

    celltype_universe: dict
    ontology: dict
    version: str

    def __init__(self, **kwargs):
        # Check that versions are consistent.
        if not list(self.celltype_universe.keys()) == list(self.ontology.keys()):
            raise ValueError(
                "error in matching versions of cell type universe and ontology in %s" %
                type(self)
            )
        # Check that ontology terms are unique also between ontologies
        if np.sum([len(x) for x in self.ontology.values()]) != \
            len(np.unique(np.array([list(x) for x in self.ontology.values()]))):
            raise ValueError(
                "duplicated ontology terms found between ontologies in %s" %
                type(self)
            )

    @property
    def versions(self):
        """
        Available cell type universe versions loaded in this instance.

        :return:
        """
        return self.celltype_universe.keys()

    def _check_version(self, version: str):
        if version not in self.celltype_universe.keys():
            raise ValueError("Version %s not found. Check self.version for available versions." % version)

    def set_version(
            self,
            version: str
    ):
        """
        Set a cell type universe version for this instance.

        :param version: Full version string "a.b.c" or celltype version "a".
        :return:
        """
        if len(version.split(".")) == 3:
            version = version.split(".")[0]
            self._check_version(version=version)
            self.version = version
        elif len(version.split(".")) == 1:
            self._check_version(version=version)
            self.version = version
        else:
            raise ValueError("version supplied should be either in format `a.b.c` or `a`")


    @property
    def ids(self):
        """
        List of all loaders understandable cell type names of this instance.

        :return:
        """
        return self._ids(self.version)

    def _ids(self, version: str):
        return [x[0] for x in self.celltype_universe[version]]

    @property
    def ontology_ids(self):
        """
        List of all cell type IDs (based on an ontology ID scheme) of this instance.

        :return:
        """
        return self._ontology_ids(self.version)

    def _ontology_ids(self, version: str):
        return [x[1] for x in self.celltype_universe[version]]

    @property
    def ntypes(self):
        """
        Number of different cell types in this instance.

        :return:
        """
        return self._ntypes(self.version)

    def _ntypes(self, version: str):
        return len(self.celltype_universe[version])

    def to_csv(
            self,
            fn: str
    ):
        pass
