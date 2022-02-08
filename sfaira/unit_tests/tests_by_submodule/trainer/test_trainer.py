import abc

import anndata
import numpy as np
import os
import pathlib
from typing import Union

from sfaira.consts.ontologies import DEFAULT_UBERON, DEFAULT_CL
from sfaira.data import load_store
from sfaira.train import TrainModelCelltype, TrainModelEmbedding
from sfaira.ui import ModelZoo
from sfaira.versions.metadata import CelltypeUniverse, OntologyCl, OntologyUberon

from sfaira.unit_tests.tests_by_submodule.estimators import HelperEstimatorBase, TARGETS
from sfaira.unit_tests import DIR_TEMP


def get_cu():
    """
    Get file name of a target universe for loading by trainer.
    """
    # Create temporary cell type universe to give to trainer.
    if not os.path.exists(DIR_TEMP):
        pathlib.Path(DIR_TEMP).mkdir(parents=True, exist_ok=True)
    fn = os.path.join(DIR_TEMP, "universe_temp.csv")
    cl = OntologyCl(branch=DEFAULT_CL)
    uberon = OntologyUberon(branch=DEFAULT_UBERON)
    cu = CelltypeUniverse(cl=cl, uberon=uberon)
    cu.write_target_universe(fn=fn, x=TARGETS)
    del cu
    return fn


class HelperTrainerBase:

    data: Union[anndata.AnnData, load_store]
    trainer: Union[TrainModelCelltype, TrainModelEmbedding]
    dir_temp = DIR_TEMP

    def __init__(self, zoo: ModelZoo):
        self.model_id = zoo.model_id
        self.tc = zoo.topology_container
        self.run_id = self.model_id + "_cv0"

    def load_adata(self, **kwargs):
        """
        This is inherited from estimator test helper.
        """
        pass

    def load_store(self, **kwargs):
        """
        This is inherited from estimator test helper.
        """
        pass

    def load_data(self, data_type):
        """
        Builds training data according to reference used in model definition.

        :param data_type:
        :return:
        """
        np.random.seed(1)
        if data_type == "adata":
            self.load_adata(organism="Homo sapiens", match_to_reference=self.tc.gc.release)
        else:
            self.load_store(organism="Homo sapiens", match_to_reference=self.tc.gc.release)

    def test_init(self, cls, estimator_kwargs: dict = {}, **kwargs):
        if not os.path.exists(self.dir_temp):
            pathlib.Path(self.dir_temp).mkdir(parents=True, exist_ok=True)
        self.load_data(data_type="adata")
        self.trainer = cls(
            data=self.data,
            model_path=os.path.join(self.dir_temp, "model"),
            **kwargs
        )
        self.trainer.zoo.model_id = self.model_id
        self.trainer.init_estim(override_hyperpar={}, **estimator_kwargs)

    def test_save(self):
        if not os.path.exists(self.dir_temp):
            pathlib.Path(self.dir_temp).mkdir(parents=True, exist_ok=True)
        self.trainer.estimator.train(
            epochs=3,
            max_steps_per_epoch=1,
            test_split=0.1,
            validation_split=0.1,
            optimizer="adam",
            lr=0.005,
        )
        self.trainer.save(
            fn=os.path.join(self.dir_temp, self.run_id), model=True, specific=True
        )


class HelperTrainer(HelperEstimatorBase, HelperTrainerBase):
    pass


def test_save_embedding():
    model_id = "embedding_homosapiens-lung-linear-0.1-0.1_mylab"
    zoo = ModelZoo()
    zoo.model_id = model_id
    test_trainer = HelperTrainer(zoo=zoo)
    test_trainer.test_init(cls=TrainModelEmbedding)
    test_trainer.test_save()


def test_save_celltypes():
    tmp_fn = get_cu()
    model_id = "celltype_homosapiens-lung-mlp-0.0.1-0.1_mylab"
    zoo = ModelZoo()
    zoo.model_id = model_id
    test_trainer = HelperTrainer(zoo=zoo)
    test_trainer.test_init(cls=TrainModelCelltype, fn_target_universe=tmp_fn)
    test_trainer.test_save()
