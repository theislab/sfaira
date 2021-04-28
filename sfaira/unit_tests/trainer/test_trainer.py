import anndata
import numpy as np
import os
import pytest
from typing import Union

from sfaira.data import DistributedStore
from sfaira.interface import ModelZoo, ModelZooCelltype, ModelZooEmbedding
from sfaira.train import TrainModelCelltype, TrainModelEmbedding
from sfaira.unit_tests.utils import cached_store_writing, simulate_anndata

dir_data = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data")
dir_meta = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data/meta")

ASSEMBLY = "Mus_musculus.GRCm38.102"
TARGETS = ["T cell", "stromal cell"]


class HelperTrainerBase:

    data: Union[anndata.AnnData, DistributedStore]
    trainer: Union[TrainModelCelltype, TrainModelEmbedding]
    zoo: ModelZoo

    def __init__(self, zoo: ModelZoo):
        self.model_id = zoo.model_id
        self.tc = zoo.topology_container

    def _simulate(self) -> anndata.AnnData:
        """
        Simulate basic data example used for unit test.

        :return: Simulated data set.
        """
        return simulate_anndata(n_obs=100, genes=self.tc.gc.ensembl, targets=TARGETS)

    def load_adata(self):
        """
        Sets attribute .data with simulated data.
        """
        self.data = self._simulate()

    def load_store(self):
        store_path = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=ASSEMBLY)
        store = DistributedStore(cache_path=store_path)
        self.data = store

    def load_data(self, data_type):
        np.random.seed(1)
        if data_type == "adata":
            self.load_adata()
        else:
            self.load_store()

    def test_for_fatal(self, cls):
        self.load_data(data_type="adata")
        trainer = cls(
            data=self.data,
            model_path=dir_meta,
        )
        trainer.zoo.set_model_id(model_id=self.model_id)
        trainer.init_estim(override_hyperpar={})


def test_for_fatal_embedding():
    model_id = "embedding_human-lung_linear_mylab_0.1_0.1"
    zoo = ModelZooEmbedding()
    zoo.set_model_id(model_id=model_id)
    test_trainer = HelperTrainerBase(zoo=zoo)
    test_trainer.test_for_fatal(cls=TrainModelEmbedding)


def test_for_fatal():
    model_id = "celltype_human-lung_mlp_mylab_0.0.1_0.1"
    zoo = ModelZooCelltype()
    zoo.set_model_id(model_id=model_id)
    test_trainer = HelperTrainerBase(zoo=zoo)
    test_trainer.test_for_fatal(cls=TrainModelCelltype)
