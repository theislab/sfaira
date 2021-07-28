import anndata
import numpy as np
import os
from typing import Union

from sfaira.data import load_store
from sfaira.ui import ModelZoo
from sfaira.unit_tests.estimators.test_estimator import TestHelperEstimatorBase, TARGETS
from sfaira.train import TrainModelCelltype, TrainModelEmbedding
from sfaira.versions.metadata import CelltypeUniverse, OntologyCl, OntologyUberon

dir_temp = os.path.join(os.path.dirname(__file__), "temp")


class HelperTrainerBase(TestHelperEstimatorBase):

    data: Union[anndata.AnnData, load_store]
    trainer: Union[TrainModelCelltype, TrainModelEmbedding]

    def __init__(self, zoo: ModelZoo):
        self.model_id = zoo.model_id
        self.tc = zoo.topology_container

    def load_data(self, data_type):
        np.random.seed(1)
        if data_type == "adata":
            self.load_adata()
        else:
            self.load_store()

    def test_init(self, cls, **kwargs):
        if not os.path.exists(dir_temp):
            os.mkdir(dir_temp)
        self.load_data(data_type="adata")
        self.trainer = cls(
            data=self.data,
            model_path=os.path.join(dir_temp, "model"),
            **kwargs
        )
        self.trainer.zoo.model_id = self.model_id
        self.trainer.init_estim(override_hyperpar={})

    def test_save(self):
        if not os.path.exists(dir_temp):
            os.mkdir(dir_temp)
        self.trainer.estimator.train(epochs=1, max_steps_per_epoch=1, test_split=0.1, validation_split=0.1,
                                     optimizer="adam", lr=0.005)
        self.trainer.save(fn=os.path.join(dir_temp, "trainer"), model=True, specific=True)


def test_save_embedding():
    model_id = "embedding_human-lung-linear-0.1-0.1_mylab"
    zoo = ModelZoo()
    zoo.model_id = model_id
    test_trainer = HelperTrainerBase(zoo=zoo)
    test_trainer.test_init(cls=TrainModelEmbedding)
    test_trainer.test_save()


def test_save_celltypes():
    # Create temporary cell type universe to give to trainer.
    tmp_fn = os.path.join(dir_temp, "universe_temp.csv")
    cl = OntologyCl(branch="v2021-02-01")
    uberon = OntologyUberon()
    cu = CelltypeUniverse(cl=cl, uberon=uberon)
    cu.write_target_universe(fn=tmp_fn, x=TARGETS)
    del cu
    model_id = "celltype_human-lung-mlp-0.0.1-0.1_mylab"
    zoo = ModelZoo()
    zoo.model_id = model_id
    test_trainer = HelperTrainerBase(zoo=zoo)
    test_trainer.test_init(cls=TrainModelCelltype, fn_target_universe=tmp_fn)
    test_trainer.test_save()
