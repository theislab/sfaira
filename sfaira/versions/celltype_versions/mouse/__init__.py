from .bladder import CelltypeVersionsMouseBladder
from .brain import CelltypeVersionsMouseBrain
from .diaphragm import CelltypeVersionsMouseDiaphragm
from .fat import CelltypeVersionsMouseFat
from .heart import CelltypeVersionsMouseHeart
from .kidney import CelltypeVersionsMouseKidney
from .large_intestine import CelltypeVersionsMouseLargeintestine
from .limb_muscle import CelltypeVersionsMouseLimbmuscle
from .liver import CelltypeVersionsMouseLiver
from .lung import CelltypeVersionsMouseLung
from .mammary_gland import CelltypeVersionsMouseMammarygland
from .marrow import CelltypeVersionsMouseMarrow
from .ovary import CelltypeVersionsMouseOvary
from .peripheral_blood import CelltypeVersionsMousePeripheralblood
from .placenta import CelltypeVersionsMousePlacenta
from .pancreas import CelltypeVersionsMousePancreas
from .prostate import CelltypeVersionsMouseProstate
from .rib import CelltypeVersionsMouseRib
from .skin import CelltypeVersionsMouseSkin
from .small_intestine import CelltypeVersionsMouseSmallintestine
from .spleen import CelltypeVersionsMouseSpleen
from .stomach import CelltypeVersionsMouseStomach
from .testis import CelltypeVersionsMouseTestis
from .thymus import CelltypeVersionsMouseThymus
from .tongue import CelltypeVersionsMouseTongue
from .trachae import CelltypeVersionsMouseTrachae
from .uterus import CelltypeVersionsMouseUterus

ORGAN_DICT = {
    "bladder": CelltypeVersionsMouseBladder(),
    "brain": CelltypeVersionsMouseBrain(),
    "diaphragm": CelltypeVersionsMouseDiaphragm(),
    "fat": CelltypeVersionsMouseFat(),
    "heart": CelltypeVersionsMouseHeart(),
    "kidney": CelltypeVersionsMouseKidney(),
    "largeintestine": CelltypeVersionsMouseLargeintestine(),
    "limbmuscle": CelltypeVersionsMouseLimbmuscle(),
    "liver": CelltypeVersionsMouseLiver(),
    "lung": CelltypeVersionsMouseLung(),
    "mammarygland": CelltypeVersionsMouseMammarygland(),
    "marrow": CelltypeVersionsMouseMarrow(),
    "ovary": CelltypeVersionsMouseOvary(),
    "peripheralblood": CelltypeVersionsMousePeripheralblood(),
    "placenta": CelltypeVersionsMousePlacenta(),
    "pancreas": CelltypeVersionsMousePancreas(),
    "prostate": CelltypeVersionsMouseProstate(),
    "rib": CelltypeVersionsMouseRib(),
    "skin": CelltypeVersionsMouseSkin(),
    "smallintestine": CelltypeVersionsMouseSmallintestine(),
    "spleen": CelltypeVersionsMouseSpleen(),
    "stomach": CelltypeVersionsMouseStomach(),
    "testis": CelltypeVersionsMouseTestis(),
    "thymus": CelltypeVersionsMouseThymus(),
    "tongue": CelltypeVersionsMouseTongue(),
    "trachae": CelltypeVersionsMouseTrachae(),
    "uterus": CelltypeVersionsMouseUterus()
}
