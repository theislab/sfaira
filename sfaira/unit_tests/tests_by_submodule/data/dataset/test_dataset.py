import numpy as np
import os
import pytest

import sfaira.versions.genomes
from sfaira.consts import AdataIdsSfaira
from sfaira.data import DatasetSuperGroup, DatasetInteractive
from sfaira.data import Universe

from sfaira.unit_tests.data_for_tests.loaders import RELEASE_HUMAN, PrepareData
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


@pytest.mark.parametrize("organ", ["lung"])
def test_dsgs_subset_dataset_wise(organ: str):
    """
    Tests if subsetting results only in datasets of the desired characteristics.
    """
    ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="organism", values=["Homo sapiens"])
    ds.subset(key="organ", values=[organ])
    ds.load()
    for x in ds.dataset_groups:
        for k, v in x.datasets.items():
            assert v.organism == "Homo sapiens", v.organism
            assert v.ontology_container_sfaira.organ.is_a(query=v.organ, reference=organ), v.organ


def test_dsgs_config_write_load():
    fn = os.path.join(DIR_TEMP, "config.csv")
    ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="organism", values=["Homo sapiens"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.write_config(fn=fn)
    ds2 = PrepareData().prepare_dsg()
    ds2.load_config(fn=fn)
    assert np.all(ds.ids == ds2.ids)


def test_dsgs_adata():
    ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="organism", values=["Homo sapiens"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    _ = ds.adata_ls


def test_dsgs_load():
    ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="organism", values=["Homo sapiens"])
    ds.subset(key="organ", values=["lung"])
    ds.load()


@pytest.mark.parametrize("celltype", ["T cell"])
def test_dsgs_subset_cell_wise(celltype: str):
    """
    Tests if sub-setting results only in datasets of the desired characteristics.
    """
    organ = "lung"
    ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="organism", values=["Homo sapiens"])
    ds.subset(key="organ", values=[organ])
    ds.load()
    ds.subset_cells(key="cell_type", values=celltype)
    for x in ds.dataset_groups:
        for k, v in x.datasets.items():
            assert v.organism == "Homo sapiens", v.id
            assert v.ontology_container_sfaira.organ.is_a(query=v.organ, reference=organ), v.organ
            for y in np.unique(v.adata.obs[v._adata_ids.cell_type].values):
                assert v.ontology_container_sfaira.cell_type.is_a(query=y, reference=celltype), y


@pytest.mark.parametrize("match_to_release", [RELEASE_HUMAN, {"Homo sapiens": RELEASE_HUMAN}])
@pytest.mark.parametrize("remove_gene_version", [False, True])
@pytest.mark.parametrize("subset_genes_to_type", [None, "protein_coding"])
def test_dsgs_streamline_features(match_to_release: str, remove_gene_version: bool, subset_genes_to_type: str):
    ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="organism", values=["Homo sapiens"])
    ds.subset(key="organ", values=["lung"])
    ds.load()
    ds.streamline_features(remove_gene_version=remove_gene_version, match_to_release=match_to_release,
                           subset_genes_to_type=subset_genes_to_type)
    gc = sfaira.versions.genomes.GenomeContainer(
        organism="Homo Sapiens",
        release=match_to_release["Homo sapiens"] if isinstance(match_to_release, dict) else match_to_release)
    gc.set(biotype=subset_genes_to_type)
    for x in ds.datasets.values():
        assert x.adata.var["gene_symbol"].tolist() == gc.symbols


def test_dsg_load():
    ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="organism", values=["Homo sapiens"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    ds.load()


def test_dsg_adata():
    ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="organism", values=["Homo sapiens"])
    ds.subset(key="organ", values=["lung"])
    ds = DatasetSuperGroup(dataset_groups=[ds])
    _ = ds.adata


def test_ds_interactive():
    adata_ids = AdataIdsSfaira()
    # Prepare object:
    ds = PrepareData().prepare_dsg(load=False)
    ds.subset(key="doi_journal", values=["no_doi_mock1"])
    ds.load()
    adata = ds.adata_ls[0]
    di = DatasetInteractive(data=adata, feature_id_col="index")
    di.organism = "Homo sapiens"
    di.organ = "lung"
    di.cell_type_obs_key = "free_annotation"
    # Test that adata is accessible in non-streamlined object:
    _ = di.adata
    # Test streamlining:
    di.streamline_features(match_to_release=RELEASE_HUMAN)
    di.streamline_metadata(schema="sfaira")
    # Test entries in streamlined object:
    adata_di = di.adata
    assert adata_ids.cell_type in adata_di.obs.columns
    assert adata_ids.cell_type + adata_ids.onto_id_suffix in adata_di.obs.columns
    assert adata_ids.organ in adata_di.uns.keys()
    assert np.all(adata_di.obs[adata_ids.cell_type].values == adata.obs["free_annotation"].values)
    assert adata_di.uns[adata_ids.organ] == "lung"
