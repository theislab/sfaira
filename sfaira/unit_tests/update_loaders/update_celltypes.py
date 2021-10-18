import numpy as np
import os
import pydoc
import pytest

from sfaira.data.dataloaders.loaders import DatasetSuperGroupLoaders
from sfaira.data import Universe

from sfaira.unit_tests.directories import DIR_SFAIRA_LOADERS


def test_celltype_update():
    """
    Warning: changes files in package if ontology mismatches are found!
    """
    # Directory choice hyperparamters:
    dir_prefix = "d"
    dir_exclude = []
    dssg = DatasetSuperGroupLoaders()
    for dsg in dssg.dataset_groups:
        collection = dsg.collection_id
        for k, ds in dsg.datasets.items():
            fn_map = str(pydoc.locate(f"sfaira.data.dataloaders.loaders.{collection}.{k}.tsv"))
            if os.path.exists(fn_map):
                # Access reading and value protection mechanisms from first data set loaded in group.
                ds.read_ontology_class_map(fn=fn_map)
                ds.write_ontology_class_map(fn=fn_map, protected_writing=False)
            else:
                print(f"did not find map for {collection}, {k}: {fn_map}")

#for d in os.listdir(DIR_SFAIRA_LOADERS):
#    if os.path.isdir(os.path.join(DIR_SFAIRA_LOADERS, d)):  # Narrow down to directories.
#        if d[:len(dir_prefix)] == dir_prefix and d not in dir_exclude:  # Narrow down to data set directories.
#            path_dsg = str(pydoc.locate(f"sfaira.data.dataloaders.loaders.{f}.FILE_PATH"))
#for f in os.listdir(d):  # Loop over files in target directory.
#    if os.path.isfile(os.path.join(d, f)):  # Narrow down to files.
#        if f.split(".")[-1] == "tsv":  # Narrow down to ontology map files.
#            tab = pd.

