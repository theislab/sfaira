from .base import CelltypeVersionsBase
from . import mouse
from . import human


mouse = mouse.ORGAN_DICT
human = human.ORGAN_DICT

# Load versions from extension if available:
try:
    from sfaira_extension.versions.celltype_versions import SPECIES_DICT as SPECIES_DICT_EXTENSION

    for organ in mouse.keys():
        if organ in SPECIES_DICT_EXTENSION["mouse"].keys():
            for v in SPECIES_DICT_EXTENSION["mouse"][organ].versions:
                if v in mouse[organ].celltype_universe.keys():
                    mouse[organ].celltype_universe[v] = {
                        **mouse[organ].celltype_universe[v],
                        **SPECIES_DICT_EXTENSION["mouse"][organ].celltype_universe[v]
                    }
                    mouse[organ].ontology[v] = {
                        **mouse[organ].ontology[v],
                        **SPECIES_DICT_EXTENSION["mouse"][organ].ontology[v]
                    }
                else:
                    mouse[organ].celltype_universe[v] = SPECIES_DICT_EXTENSION["mouse"][organ].celltype_universe[v]
                    mouse[organ].ontology[v] = SPECIES_DICT_EXTENSION["mouse"][organ].ontology[v]

    for organ in human.keys():
        if organ in SPECIES_DICT_EXTENSION["human"].keys():
            for v in SPECIES_DICT_EXTENSION["human"][organ].versions:
                if v in human[organ].celltype_universe.keys():
                    human[organ].celltype_universe[v] = {
                        **human[organ].celltype_universe[v],
                        **SPECIES_DICT_EXTENSION["human"][organ].celltype_universe[v]
                    }
                    human[organ].ontology[v] = {
                        **human[organ].ontology[v],
                        **SPECIES_DICT_EXTENSION["human"][organ].ontology[v]
                    }
                else:
                    human[organ].celltype_universe[v] = SPECIES_DICT_EXTENSION["human"][organ].celltype_universe[v]
                    human[organ].ontology[v] = SPECIES_DICT_EXTENSION["human"][organ].ontology[v]
except ImportError:
    pass

SPECIES_DICT = {
    "mouse": mouse,
    "human": human
}