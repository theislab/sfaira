from .base import CelltypeVersionsBase
from . import mouse
from . import human


mouse = mouse.ORGAN_DICT
human = human.ORGAN_DICT

# Load versions from extension if available:
try:
    from sfaira_extension.versions.celltype_versions import SPECIES_DICT as SPECIES_DICT_EXTENSION
    mouse_e = SPECIES_DICT_EXTENSION["mouse"]
    human_e = SPECIES_DICT_EXTENSION["human"]
    for k in mouse.keys():
        if k in mouse_e.keys():
            mouse[k].celltype_universe.update(mouse_e[k])
            mouse[k].ontology.update(mouse_e[k])
        if k in mouse_e.keys():
            human[k].celltype_universe.update(human_e[k])
            human[k].ontology.update(human_e[k])
except ImportError:
    pass

SPECIES_DICT = {
    "mouse": mouse,
    "human": human
}