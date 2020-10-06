import sys
import tensorflow as tf
import sfaira.api as sfaira
import os

print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

path_base = str(sys.argv[1])
fn = str(sys.argv[2])
genome = str(sys.argv[3])

path = os.path.join(path_base, "data")
path_meta = os.path.join(path_base, "meta_leander")

from sfaira.data import human
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
ds = sfaira.data.DatasetSuperGroup(
    dataset_groups=[ds_dict[k] for k in list(ds_dict.keys())]
)
ds.load_all_tobacked(
    fn_backed=fn,
    genome=genome,
    shuffled=False,
    as_dense=False,
    annotated_only=False
)
