import sys
import tensorflow as tf
from sfaira.data import human


print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

path = str(sys.argv[1])
path_meta = str(sys.argv[2])

ds_dict = {
    "adipose": human.DatasetGroupAdipose(path=path, meta_path=path_meta),
    "adrenalgland": human.DatasetGroupAdrenalgland(path=path, meta_path=path_meta),
    "mixed": human.DatasetGroupMixed(path=path, meta_path=path_meta),
    "artery": human.DatasetGroupArtery(path=path, meta_path=path_meta),
    "bladder": human.DatasetGroupBladder(path=path, meta_path=path_meta),
    "blood": human.DatasetGroupBlood(path=path, meta_path=path_meta),
    "bone": human.DatasetGroupBone(path=path, meta_path=path_meta),
    "brain": human.DatasetGroupBrain(path=path, meta_path=path_meta),
    "calvaria": human.DatasetGroupCalvaria(path=path, meta_path=path_meta),
    "cervix": human.DatasetGroupCervix(path=path, meta_path=path_meta),
    "chorionicvillus": human.DatasetGroupChorionicvillus(path=path, meta_path=path_meta),
    "colon": human.DatasetGroupColon(path=path, meta_path=path_meta),
    "duodenum": human.DatasetGroupDuodenum(path=path, meta_path=path_meta),
    "epityphlon": human.DatasetGroupEpityphlon(path=path, meta_path=path_meta),
    "esophagus": human.DatasetGroupEsophagus(path=path, meta_path=path_meta),
    "eye": human.DatasetGroupEye(path=path, meta_path=path_meta),
    "fallopiantube": human.DatasetGroupFallopiantube(path=path, meta_path=path_meta),
    "femalegonad": human.DatasetGroupFemalegonad(path=path, meta_path=path_meta),
    "gallbladder": human.DatasetGroupGallbladder(path=path, meta_path=path_meta),
    "heart": human.DatasetGroupHeart(path=path, meta_path=path_meta),
    "hesc": human.DatasetGroupHesc(path=path, meta_path=path_meta),
    "ileum": human.DatasetGroupIleum(path=path, meta_path=path_meta),
    "jejunum": human.DatasetGroupJejunum(path=path, meta_path=path_meta),
    "kidney": human.DatasetGroupKidney(path=path, meta_path=path_meta),
    "liver": human.DatasetGroupLiver(path=path, meta_path=path_meta),
    "lung": human.DatasetGroupLung(path=path, meta_path=path_meta),
    "malegonad": human.DatasetGroupMalegonad(path=path, meta_path=path_meta),
    "muscle": human.DatasetGroupMuscle(path=path, meta_path=path_meta),
    "omentum": human.DatasetGroupOmentum(path=path, meta_path=path_meta),
    "pancreas": human.DatasetGroupPancreas(path=path, meta_path=path_meta),
    "placenta": human.DatasetGroupPlacenta(path=path, meta_path=path_meta),
    "pleura": human.DatasetGroupPleura(path=path, meta_path=path_meta),
    "prostate": human.DatasetGroupProstate(path=path, meta_path=path_meta),
    "rectum": human.DatasetGroupRectum(path=path, meta_path=path_meta),
    "rib": human.DatasetGroupRib(path=path, meta_path=path_meta),
    "skin": human.DatasetGroupSkin(path=path, meta_path=path_meta),
    "spinalcord": human.DatasetGroupSpinalcord(path=path, meta_path=path_meta),
    "spleen": human.DatasetGroupSpleen(path=path, meta_path=path_meta),
    "stomach": human.DatasetGroupStomach(path=path, meta_path=path_meta),
    "thymus": human.DatasetGroupThymus(path=path, meta_path=path_meta),
    "thyroid": human.DatasetGroupThyroid(path=path, meta_path=path_meta),
    "trachea": human.DatasetGroupTrachea(path=path, meta_path=path_meta),
    "ureter": human.DatasetGroupUreter(path=path, meta_path=path_meta),
    "uterus": human.DatasetGroupUterus(path=path, meta_path=path_meta),
}
for k in list(ds_dict.keys()):
    for kk in ds_dict[k].ids:
        ds_dict[k].datasets[kk].write_meta(dir_out=path_meta)
