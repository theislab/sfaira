import numpy as np
import os
import pandas as pd
import pydoc
import pytest

from sfaira.data.dataloaders.loaders import DatasetSuperGroupLoaders

from sfaira.unit_tests.directories import DIR_DATA_LOADERS_CACHE

SYMBOL_COL = "target"
ID_COL = "target_id"
UNKNOWN_IDS = ["UNKNOWN", "NOT_A_CELL"]


def test_celltype_update():
    """
    Warning: changes files in package if ontology mismatches are found!
    """
    dssg = DatasetSuperGroupLoaders(data_path=DIR_DATA_LOADERS_CACHE, meta_path=DIR_DATA_LOADERS_CACHE,
                                    cache_path=DIR_DATA_LOADERS_CACHE)
    for dsg in dssg.dataset_groups:
        collection = dsg.collection_id
        for ds in dsg.datasets.values():
            path_dsg = os.path.dirname(str(pydoc.locate(f"sfaira.data.dataloaders.loaders.{collection}.FILE_PATH")))
            for f in os.listdir(path_dsg):
                if f.split(".")[-1] == "tsv":  # Narrow down to ontology map files.
                    fn_map = os.path.join(path_dsg, f)
                    if os.path.exists(fn_map):
                        # Access reading and value protection mechanisms from first data set loaded in group.
                        tab = ds._read_ontology_class_map(fn=fn_map)
                        tab[SYMBOL_COL] = [
                            ds.celltypes_universe.onto_cl.convert_to_name(x)
                            if x not in UNKNOWN_IDS else x
                            for x in tab[ID_COL].values
                        ]
                        tab.to_csv(path_or_buf=fn_map, index=False, sep="\t")
                    else:
                        print(f"did not find map for {collection}, {ds.fn_ontology_class_map_tsv}: {fn_map}")
