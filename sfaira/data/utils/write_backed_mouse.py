import sys
import tensorflow as tf
import sfaira
import os

from sfaira.data import mouse


print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

path = str(sys.argv[1])
fn = str(sys.argv[2])
genome = str(sys.argv[3])

path_meta = os.path.join(path, "meta")

ds_dict = {
    "bladder": mouse.DatasetGroupBladder(path=path, meta_path=path_meta, cache_path=path_meta),
    "brain": mouse.DatasetGroupBrain(path=path, meta_path=path_meta, cache_path=path_meta),
    "diaphragm": mouse.DatasetGroupDiaphragm(path=path, meta_path=path_meta, cache_path=path_meta),
    "adipose": mouse.DatasetGroupAdipose(path=path, meta_path=path_meta, cache_path=path_meta),
    "heart": mouse.DatasetGroupHeart(path=path, meta_path=path_meta, cache_path=path_meta),
    "kidney": mouse.DatasetGroupKidney(path=path, meta_path=path_meta, cache_path=path_meta),
    "colon": mouse.DatasetGroupColon(path=path, meta_path=path_meta, cache_path=path_meta),
    "muscle": mouse.DatasetGroupMuscle(path=path, meta_path=path_meta, cache_path=path_meta),
    "liver": mouse.DatasetGroupLiver(path=path, meta_path=path_meta, cache_path=path_meta),
    "lung": mouse.DatasetGroupLung(path=path, meta_path=path_meta, cache_path=path_meta),
    "mammarygland": mouse.DatasetGroupMammaryGland(path=path, meta_path=path_meta, cache_path=path_meta),
    "bone": mouse.DatasetGroupBone(path=path, meta_path=path_meta, cache_path=path_meta),
    "femalegonad": mouse.DatasetGroupFemalegonad(path=path, meta_path=path_meta, cache_path=path_meta),
    "pancreas": mouse.DatasetGroupPancreas(path=path, meta_path=path_meta, cache_path=path_meta),
    "blood": mouse.DatasetGroupBlood(path=path, meta_path=path_meta, cache_path=path_meta),
    "placenta": mouse.DatasetGroupPlacenta(path=path, meta_path=path_meta, cache_path=path_meta),
    "prostate": mouse.DatasetGroupProstate(path=path, meta_path=path_meta, cache_path=path_meta),
    "rib": mouse.DatasetGroupRib(path=path, meta_path=path_meta, cache_path=path_meta),
    "skin": mouse.DatasetGroupSkin(path=path, meta_path=path_meta, cache_path=path_meta),
    "ileum": mouse.DatasetGroupIleum(path=path, meta_path=path_meta, cache_path=path_meta),
    "spleen": mouse.DatasetGroupSpleen(path=path, meta_path=path_meta, cache_path=path_meta),
    "stomach": mouse.DatasetGroupStomach(path=path, meta_path=path_meta, cache_path=path_meta),
    "malegonad": mouse.DatasetGroupMalegonad(path=path, meta_path=path_meta, cache_path=path_meta),
    "thymus": mouse.DatasetGroupThymus(path=path, meta_path=path_meta, cache_path=path_meta),
    "tongue": mouse.DatasetGroupTongue(path=path, meta_path=path_meta, cache_path=path_meta),
    "trachea": mouse.DatasetGroupTrachea(path=path, meta_path=path_meta, cache_path=path_meta),
    "uterus": mouse.DatasetGroupUterus(path=path, meta_path=path_meta, cache_path=path_meta),
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
