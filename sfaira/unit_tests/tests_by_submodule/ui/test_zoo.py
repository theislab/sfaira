from sfaira.ui import ModelZoo


def test_for_fatal_embedding():
    model_id = "embedding_homosapiens-lung-linear-0.1-0.1_mylab"
    zoo = ModelZoo()
    zoo.model_id = model_id
    assert zoo.model_id == model_id
    assert zoo.model_class == "embedding"
    assert zoo.model_name == "homosapiens-lung-linear-0.1-0.1"
    assert zoo.organisation == "mylab"
    _ = zoo.topology_container
    _ = zoo.topology_container.topology
    _ = zoo.topology_container.gc


def test_for_fatal_celltype():
    model_id = "celltype_homosapiens-lung-mlp-0.0.1-0.1_mylab"
    zoo = ModelZoo()
    zoo.model_id = model_id
    assert zoo.model_id == model_id
    assert zoo.model_class == "celltype"
    assert zoo.model_name == "homosapiens-lung-mlp-0.0.1-0.1"
    assert zoo.organisation == "mylab"
    _ = zoo.topology_container
    _ = zoo.topology_container.topology
    _ = zoo.topology_container.gc
