import numpy as np

from sfaira.data import load_store
from sfaira.versions.genomes.genomes import GenomeContainer

from sfaira.unit_tests.data_for_tests.loaders import RELEASE_MOUSE, PrepareData


def _get_cart(store_format, feature_space, map_fn=None, batch_size=1, **kwargs):
    """
    TODO move into PrepareData?
    """
    store_path = PrepareData().prepare_store(store_format=store_format)
    store = load_store(cache_path=store_path, store_format=store_format)
    store.subset(attr_key="organism", values=["Mus musculus"])
    if feature_space == "single":
        store = store.stores["Mus musculus"]
        if "idx" in kwargs.keys() and isinstance(kwargs["idx"], dict):
            assert "Mus musculus" in kwargs["idx"].keys()
            kwargs["idx"] = kwargs["idx"]["Mus musculus"]
    gc = GenomeContainer(release=RELEASE_MOUSE, organism="Mus musculus")
    gc.set(**{"biotype": "protein_coding"})
    store.genome_container = gc

    if map_fn is None:

        def map_fn(x_, obs_):
            return (np.asarray(x_), ),

    g = store.checkout(map_fn=map_fn, batch_size=batch_size, **kwargs)
    return g
