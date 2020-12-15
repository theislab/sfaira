from .base import CelltypeVersionsBase
from . import mouse
from . import human


mouse = mouse.ORGAN_DICT
human = human.ORGAN_DICT

# Load versions from extension if available:
try:
    from sfaira_extension.versions.celltype_versions import ORGANISM_DICT as ORGANISM_DICT_EXTENSION

    for organ in mouse.keys():
        if organ in ORGANISM_DICT_EXTENSION["mouse"].keys():
            for v in ORGANISM_DICT_EXTENSION["mouse"][organ].versions:
                if v in mouse[organ].celltype_universe.keys():
                    raise ValueError(f'Celltype version {v} already defined for mouse organ {organ} in base sfaira. '
                                     f'Please define a new version in sfaira_extension.')
                else:
                    mouse[organ].celltype_universe[v] = ORGANISM_DICT_EXTENSION["mouse"][organ].celltype_universe[v]
                    mouse[organ].ontology[v] = ORGANISM_DICT_EXTENSION["mouse"][organ].ontology[v]

    for organ in human.keys():
        if organ in ORGANISM_DICT_EXTENSION["human"].keys():
            for v in ORGANISM_DICT_EXTENSION["human"][organ].versions:
                if v in human[organ].celltype_universe.keys():
                    raise ValueError(f'Celltype version {v} already defined for human organ {organ} in base sfaira. '
                                     f'Please define a new version in sfaira_extension.')
                else:
                    human[organ].celltype_universe[v] = ORGANISM_DICT_EXTENSION["human"][organ].celltype_universe[v]
                    human[organ].ontology[v] = ORGANISM_DICT_EXTENSION["human"][organ].ontology[v]
except ImportError:
    pass

ORGANISM_DICT = {
    "mouse": mouse,
    "human": human
}