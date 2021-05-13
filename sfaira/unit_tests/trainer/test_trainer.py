import anndata
import numpy as np
import os
from typing import Union

from sfaira.data import load_store
from sfaira.interface import ModelZoo
from sfaira.train import TrainModelCelltype, TrainModelEmbedding
from sfaira.unit_tests.utils import cached_store_writing, simulate_anndata

dir_data = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data")
dir_meta = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_data/meta")

ASSEMBLY = "Mus_musculus.GRCm38.102"
TARGETS = ["T cell", "stromal cell"]


class HelperTrainerBase:

    data: Union[anndata.AnnData, load_store]
    trainer: Union[TrainModelCelltype, TrainModelEmbedding]

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
        store_path = cached_store_writing(dir_data=dir_data, dir_meta=dir_meta, assembly=ASSEMBLY, organism="mouse")
        store = load_store(cache_path=store_path)
        self.data = store

    def load_data(self, data_type):
        np.random.seed(1)
        if data_type == "adata":
            self.load_adata()
        else:
            self.load_store()

    def test_init(self, cls):
        self.load_data(data_type="adata")
        self.trainer = cls(
            data=self.data,
            model_path=dir_meta,
        )
        self.trainer.zoo.model_id = self.model_id
        self.trainer.init_estim(override_hyperpar={})

    def test_save(self):
        self.trainer.estimator.train(epochs=1, max_steps_per_epoch=1, test_split=0.1, validation_split=0.1,
                                     optimizer="adam", lr=0.005)
        self.trainer.save(fn=os.path.join(dir_data, "trainer_test"), model=True, specific=True)


def test_save_embedding():
    model_id = "embedding_human-lung-linear-0.1-0.1_mylab"
    zoo = ModelZoo()
    zoo.model_id = model_id
    test_trainer = HelperTrainerBase(zoo=zoo)
    test_trainer.test_init(cls=TrainModelEmbedding)
    test_trainer.test_save()


def test_save_celltypes():
    model_id = "celltype_human-lung-mlp-0.0.1-0.1_mylab"
    zoo = ModelZoo()
    zoo.model_id = model_id
    test_trainer = HelperTrainerBase(zoo=zoo)
    test_trainer.test_init(cls=TrainModelCelltype)
    test_trainer.test_save()
