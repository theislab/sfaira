from .bladder import CelltypeVersionsMouseBladder
from .brain import CelltypeVersionsMouseBrain
from .diaphragm import CelltypeVersionsMouseDiaphragm
from .adipose import CelltypeVersionsMouseAdipose
from .heart import CelltypeVersionsMouseHeart
from .kidney import CelltypeVersionsMouseKidney
from .colon import CelltypeVersionsMouseColon
from .muscle import CelltypeVersionsMouseMuscle
from .liver import CelltypeVersionsMouseLiver
from .lung import CelltypeVersionsMouseLung
from .mammarygland import CelltypeVersionsMouseMammarygland
from .bone import CelltypeVersionsMouseBone
from .ovary import CelltypeVersionsMouseOvary
from .blood import CelltypeVersionsMouseBlood
from .placenta import CelltypeVersionsMousePlacenta
from .pancreas import CelltypeVersionsMousePancreas
from .prostate import CelltypeVersionsMouseProstate
from .rib import CelltypeVersionsMouseRib
from .skin import CelltypeVersionsMouseSkin
from .ileum import CelltypeVersionsMouseIleum
from .spleen import CelltypeVersionsMouseSpleen
from .stomach import CelltypeVersionsMouseStomach
from .malegonad import CelltypeVersionsMouseMalegonad
from .thymus import CelltypeVersionsMouseThymus
from .tongue import CelltypeVersionsMouseTongue
from .trachea import CelltypeVersionsMouseTrachea
from .uterus import CelltypeVersionsMouseUterus

ORGAN_DICT = {
    "bladder": CelltypeVersionsMouseBladder(),
    "brain": CelltypeVersionsMouseBrain(),
    "diaphragm": CelltypeVersionsMouseDiaphragm(),
    "adipose": CelltypeVersionsMouseAdipose(),
    "heart": CelltypeVersionsMouseHeart(),
    "kidney": CelltypeVersionsMouseKidney(),
    "colon": CelltypeVersionsMouseColon(),
    "muscle": CelltypeVersionsMouseMuscle(),
    "liver": CelltypeVersionsMouseLiver(),
    "lung": CelltypeVersionsMouseLung(),
    "mammarygland": CelltypeVersionsMouseMammarygland(),
    "bone": CelltypeVersionsMouseBone(),
    "ovary": CelltypeVersionsMouseOvary(),
    "blood": CelltypeVersionsMouseBlood(),
    "placenta": CelltypeVersionsMousePlacenta(),
    "pancreas": CelltypeVersionsMousePancreas(),
    "prostate": CelltypeVersionsMouseProstate(),
    "rib": CelltypeVersionsMouseRib(),
    "skin": CelltypeVersionsMouseSkin(),
    "ileum": CelltypeVersionsMouseIleum(),
    "spleen": CelltypeVersionsMouseSpleen(),
    "stomach": CelltypeVersionsMouseStomach(),
    "malegonad": CelltypeVersionsMouseMalegonad(),
    "thymus": CelltypeVersionsMouseThymus(),
    "tongue": CelltypeVersionsMouseTongue(),
    "trachea": CelltypeVersionsMouseTrachea(),
    "uterus": CelltypeVersionsMouseUterus()
}
