import os
import sys
import tensorflow as tf

# Any data loader here to extract path:
from sfaira.data import DatasetSuperGroupSfaira

print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

config_path = str(sys.argv[1])


def clean(s):
    if s is not None:
        s = s.replace(' ', '').replace('-', '').replace('_', '').replace("'", '').lower()
    return s


configs_to_write = {
    "human": [
        "adipose tissue",
        "adrenal gland",
        "artery",
        "blood",
        "bone marrow",
        "brain",
        "chorionic villus",
        "diaphragm",
        "esophagus",
        "eye",
        "gall bladder",
        "heart",
        "intestine",
        "kidney",
        "liver",
        "lung",
        "muscle organ",
        "ovary",
        "pancreas",
        "placenta",
        "pleura",
        "prostate gland",
        "rib"
        "skeleton",
        "skin of body",
        "spinal cord",
        "spleen",
        "stomach",
        "testis",
        "tongue",
        "thymus",
        "thyroid gland",
        "trachea",
        "ureter",
        "urinary bladder",
        "uterine cervix",
        "uterus",
        "vault of skull",
    ],
    "mouse": [
        "adipose tissue",
        "blood",
        "bone marrow",
        "brain",
        "diaphragm",
        "heart",
        "intestine",
        "kidney",
        "liver",
        "lung",
        "mammary gland",
        "muscle organ",
        "ovary",
        "pancreas",
        "placenta",
        "prostate gland",
        "skin of body",
        "spleen",
        "stomach",
        "testis",
        "thymus",
        "tongue",
        "trachea",
        "urinary bladder",
        "uterus",
    ]
}

for organism, organs in configs_to_write.items():
    for organ in organs:
        dsgs = DatasetSuperGroupSfaira(
            data_path=".",
            meta_path=".",
            cache_path="."
        )
        dsgs.subset(key="organism", values=[organism])
        dsgs.subset(key="organ", values=[organ])
        dsgs.write_config(os.path.join(config_path, f"config_{clean(organism)}_{clean(organ)}.csv"))
