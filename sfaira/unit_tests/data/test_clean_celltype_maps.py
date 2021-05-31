from sfaira.data.dataloaders.loaders import DatasetSuperGroupLoaders


def test_map_celltype_to_ontology():
    # Paths do not matter here as data sets are not loaded for these operations.
    dsgl = DatasetSuperGroupLoaders(
        data_path="~",
        meta_path="~",
        cache_path="~"
    )
    for x in dsgl.dataset_groups:
        print(x.ids)
        x.clean_ontology_class_map()
