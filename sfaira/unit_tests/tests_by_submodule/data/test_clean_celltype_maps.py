from sfaira.consts import OC
from sfaira.data.dataloaders.loaders import DatasetSuperGroupLoaders


# TODO export this into a maintenance module.
def maintenance_test_map_celltype_to_ontology(version: str = "3_0_0"):
    # Paths do not matter here as data sets are not loaded for these operations.
    dsgl = DatasetSuperGroupLoaders(data_path="", meta_path="", cache_path="")
    OC.set_schema_version(version=version)
    for x in dsgl.dataset_groups:
        x.update_ontology_symbols()
