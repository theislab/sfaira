import sys
import tensorflow as tf
from sfaira.data import human


print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

path = str(sys.argv[1])
path_meta = str(sys.argv[2])

ds_dict = {
    "bladder": mouse.DatasetGroupBladder(path=path, meta_path=path_meta),
    "brain": mouse.DatasetGroupBrain(path=path, meta_path=path_meta),
    "diaphragm": mouse.DatasetGroupDiaphragm(path=path, meta_path=path_meta),
    "fat": mouse.DatasetGroupFat(path=path, meta_path=path_meta),
    "heart": mouse.DatasetGroupHeart(path=path, meta_path=path_meta),
    "kidney": mouse.DatasetGroupKidney(path=path, meta_path=path_meta),
    "largeintestine": mouse.DatasetGroupLargeintestine(path=path, meta_path=path_meta),
    "limbmuscle": mouse.DatasetGroupLimbmuscle(path=path, meta_path=path_meta),
    "liver": mouse.DatasetGroupLiver(path=path, meta_path=path_meta),
    "lung": mouse.DatasetGroupLung(path=path, meta_path=path_meta),
    "mammarygland": mouse.DatasetGroupMammaryGland(path=path, meta_path=path_meta),
    "marrow": mouse.DatasetGroupMarrow(path=path, meta_path=path_meta),
    "ovary": mouse.DatasetGroupOvary(path=path, meta_path=path_meta),
    "pancreas": mouse.DatasetGroupPancreas(path=path, meta_path=path_meta),
    "peripheralblood": mouse.DatasetGroupPeripheralBlood(path=path, meta_path=path_meta),
    "placenta": mouse.DatasetGroupPlacenta(path=path, meta_path=path_meta),
    "prostate": mouse.DatasetGroupProstate(path=path, meta_path=path_meta),
    "rib": mouse.DatasetGroupRib(path=path, meta_path=path_meta),
    "skin": mouse.DatasetGroupSkin(path=path, meta_path=path_meta),
    "smallintestine": mouse.DatasetGroupSmallintestine(path=path, meta_path=path_meta),
    "spleen": mouse.DatasetGroupSpleen(path=path, meta_path=path_meta),
    "stomach": mouse.DatasetGroupStomach(path=path, meta_path=path_meta),
    "testis": mouse.DatasetGroupTestis(path=path, meta_path=path_meta),
    "thymus": mouse.DatasetGroupThymus(path=path, meta_path=path_meta),
    "tongue": mouse.DatasetGroupTongue(path=path, meta_path=path_meta),
    "trachae": mouse.DatasetGroupTrachea(path=path, meta_path=path_meta),
    "uterus": mouse.DatasetGroupUterus(path=path, meta_path=path_meta)
}
for k in list(ds_dict.keys()):
    for kk in ds_dict[k].ids:
        ds_dict[k].datasets[kk].write_meta(dir_out=path_meta)
