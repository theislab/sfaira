import os
import pydoc
import sfaira
import sys
import tensorflow as tf

# Any data loader here to extract path:
from sfaira.data.dataloaders.loaders.d10_1016_j_cmet_2019_01_021 import FILE_PATH

print(tf.__version__)

# Set global variables.
print("sys.argv", sys.argv)

path = str(sys.argv[1])
path_meta = str(sys.argv[2])
path_cache = str(sys.argv[3])
processes = int(str(sys.argv[4]))

# Data loaders with one file per data set rather than one summary file:
# ToDo: not supported yet.
# TMS d10_1101_661728
# MCA d10_1016_j_cell_2018_02_001
# HCL d10_1038_s41586_020_2157_4
studys_separate_csvs = []
#    "d10_1101_661728",
#    "d10_1016_j_cell_2018_02_001",
#    "d10_1038_s41586_020_2157_4"
# ]

dir_sfaira_dataloaders = "/" + str(os.path.join(*str(os.path.dirname(FILE_PATH)).split("/")[:-1]))

dir_prefix = "d"
dir_exlcude = []
for dir_study in os.listdir(dir_sfaira_dataloaders):
    print(dir_study)
    if os.path.isdir(os.path.join(dir_sfaira_dataloaders, dir_study)):  # only directories
        # Narrow down to data set directories:
        if dir_study[:len(dir_prefix)] == dir_prefix and dir_study not in dir_exlcude:
            for f_dataset in os.listdir(os.path.join(dir_sfaira_dataloaders, dir_study)):
                if os.path.isfile(os.path.join(dir_sfaira_dataloaders, dir_study, f_dataset)):  # only files
                    print(f_dataset)
                    # Narrow down to data set files:
                    if f_dataset.split(".")[-1] == "py" and \
                            f_dataset.split(".")[0] not in ["__init__", "base", "group"]:
                        file_module = ".".join(str(f_dataset).split(".")[:-1])
                        DatasetFound = pydoc.locate(
                            "sfaira.data.dataloaders.loaders." + dir_study + "." + file_module + ".Dataset")
                        # Check if global objects are available:
                        # - SAMPLE_FNS: for DatasetBaseGroupLoadingManyFiles
                        # - SAMPLE_IDS: for DatasetBaseGroupLoadingOneFile
                        sample_fns = pydoc.locate(
                            "sfaira.data.dataloaders.loaders." + dir_study + "." + file_module + ".SAMPLE_FNS")
                        sample_ids = pydoc.locate(
                            "sfaira.data.dataloaders.loaders." + dir_study + "." + file_module + ".SAMPLE_IDS")
                        if sample_fns is not None and sample_ids is None:
                            # DatasetBaseGroupLoadingManyFiles:
                            datasets_f = [
                                DatasetFound(
                                    sample_fn=x,
                                    path=path,
                                    meta_path=path_meta,
                                    cache_path=path_cache
                                )
                                for x in sample_fns
                            ]
                        elif sample_fns is None and sample_ids is not None:
                            # DatasetBaseGroupLoadingManyFiles:
                            datasets_f = [
                                DatasetFound(
                                    sample_id=x,
                                    path=path,
                                    meta_path=path_meta,
                                    cache_path=path_cache
                                )
                                for x in sample_ids
                            ]
                        elif sample_fns is not None and sample_ids is not None:
                            raise ValueError(f"sample_fns and sample_ids both found for {f_dataset}")
                        else:
                            datasets_f = [DatasetFound(
                                path=path,
                                meta_path=path_meta,
                                cache_path=path_cache
                            )]
                        dsg_f = sfaira.data.DatasetGroup(datasets=dict([(x.id, x) for x in datasets_f]))
                        dsg_f.load(
                            load_raw=False,
                            allow_caching=True
                        )
                        if str(dir_study) in studys_separate_csvs:
                            pass
                        else:
                            dsg_f.write_ontology_class_map(
                                fn=os.path.join(dir_sfaira_dataloaders, dir_study, file_module + ".csv"),
                                protected_writing=True,
                                n_suggest=4,
                            )
