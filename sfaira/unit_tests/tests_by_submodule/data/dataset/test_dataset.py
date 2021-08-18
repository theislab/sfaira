import numpy as np
import os
import pytest

from sfaira.data import DatasetSuperGroup
from sfaira.data import Universe

from sfaira.unit_tests.data_for_tests.loaders import ASSEMBLY_MOUSE, prepare_dsg
from sfaira.unit_tests.directories import DIR_TEMP, DIR_DATA_LOADERS_CACHE


def test_dsgs_instantiate():
    _ = Universe(data_path=DIR_DATA_LOADERS_CACHE, meta_path=DIR_DATA_LOADERS_CACHE, cache_path=DIR_DATA_LOADERS_CACHE)


def test_dsgs_crossref():
    """
    Tests if crossref attributes can be retrieved for all data loader entries with DOI journal defined.
    Attributes tested:
        - title
    """
    universe = Universe(data_path=DIR_DATA_LOADERS_CACHE, meta_path=DIR_DATA_LOADERS_CACHE,
                        cache_path=DIR_DATA_LOADERS_CACHE)
    for k, v in universe.datasets.items():
        title = v.title
        if title is None:
            if v.doi_journal is not None and "no_doi" not in v.doi_journal:
                raise ValueError(f"did not retrieve title for data set {k} with DOI: {v.doi_journal}.")


@pytest.mark.parametrize("organ", ["intestine", "ileum"])
def test_dsgs_subset_dataset_wise(organ: str):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    ds = prepare_dsg(load=False)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=[organ])
    ds.load()
    for x in ds.dataset_groups:
        for k, v in x.datasets.items():
            assert v.organism == "mouse", v.organism
            assert v.ontology_container_sfaira.organ.is_a(query=v.organ, reference=organ), v.organ


def test_dsgs_config_write_load():
    fn = os.path.join(DIR_TEMP, "config.csv")
    ds = prepare_dsg(load=False)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.write_config(fn=fn)
    ds2 = prepare_dsg()
    ds2.load_config(fn=fn)
    assert np.all(ds.ids == ds2.ids)


def test_dsgs_adata():
    ds = prepare_dsg(load=False)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    _ = ds.adata


def test_dsgs_load():
    ds = prepare_dsg(load=False)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()


@pytest.mark.parametrize("organ", ["lung"])
@pytest.mark.parametrize("celltype", ["T cell"])
def test_dsgs_subset_cell_wise(organ: str, celltype: str):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    ds = prepare_dsg(load=False)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=[organ])
    ds.load()
    ds.subset_cells(key="cellontology_class", values=celltype)
    for x in ds.dataset_groups:
        for k, v in x.datasets.items():
            assert v.organism == "mouse", v.id
            assert v.ontology_container_sfaira.organ.is_a(query=v.organ, reference=organ), v.organ
            for y in np.unique(v.adata.obs[v._adata_ids.cell_type].values):
                assert v.ontology_container_sfaira.cell_type.is_a(query=y, reference=celltype), y


@pytest.mark.parametrize("match_to_reference", ["Mus_musculus.GRCm38.102", {"mouse": ASSEMBLY_MOUSE}])
@pytest.mark.parametrize("remove_gene_version", [False, True])
@pytest.mark.parametrize("subset_genes_to_type", [None, "protein_coding"])
def test_dsgs_streamline_features(match_to_reference: str, remove_gene_version: bool, subset_genes_to_type: str):
    ds = prepare_dsg(load=False)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.streamline_features(remove_gene_version=remove_gene_version, match_to_reference=match_to_reference,
                           subset_genes_to_type=subset_genes_to_type)


def test_dsg_load():
    ds = prepare_dsg(load=False)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    ds.load()


def test_dsg_adata():
    ds = prepare_dsg(load=False)
    ds.subset(key="organism", values=["mouse"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    _ = ds.adata
