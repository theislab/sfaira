import os
import sys
import tensorflow as tf

# Any data loader here to extract path:
from sfaira.data import DistributedStore, clean_string

print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

store_path = str(sys.argv[1])
config_path = str(sys.argv[2])


configs_to_write = {
    "human": [
        "adipose tissue",
        "adrenal gland",
        "artery",
        "blood",
        "bone marrow",
        "brain",
        "chorionic villus",
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
        "rib",
        "skeleton",
        "skin of body",
        "spinal cord",
        "spleen",
        "stomach",
        "testis",
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
    ],
}

for organism, organs in configs_to_write.items():
    for organ in organs:
        print(f"Writing {organism} {organ}")
        store = DistributedStore(cache_path=store_path)
        store.subset(attr_key="organism", values=[organism])
        store.subset(attr_key="organ", values=[organ])
        store.write_config(os.path.join(config_path, f"config_{clean_string(organism)}_{clean_string(organ)}"))
