import abc
import anndata
import numpy as np
import pandas as pd
import pickle
from typing import Union

from sfaira.consts import AdataIdsSfaira
from sfaira.data.store.base import DistributedStoreBase
from sfaira.data import DistributedStoreSingleFeatureSpace, Universe
from sfaira.estimators import EstimatorKeras, EstimatorKerasCelltype, EstimatorKerasEmbedding
from sfaira.ui import ModelZoo


class TrainModel:

    data: Union[anndata.AnnData, DistributedStoreSingleFeatureSpace]
    estimator: Union[EstimatorKeras, None]

    def __init__(
            self,
            data: Union[str, anndata.AnnData, Universe, DistributedStoreSingleFeatureSpace],
    ):
        # Check if handling backed anndata or base path to directory of raw files:
        if isinstance(data, str) and data.split(".")[-1] == "h5ad":
            self.data = anndata.read(data, backed='r')
            if len(self.data.obs.columns) == 0:
                fn_backed_obs = ".".join(data.split(".")[:-1]) + "_obs.csv"
                self.data.obs = pd.read_csv(fn_backed_obs)
        elif isinstance(data, anndata.AnnData):
            self.data = data
        elif isinstance(data, list) and isinstance(data[0], anndata.AnnData):
            self.data = data
        elif isinstance(data, Universe):
            self.data = data.adata
        elif isinstance(data, DistributedStoreBase):
            self.data = data
        else:
            raise ValueError(f"did not recognize data of type {type(data)}")
        self.zoo = ModelZoo()
        self._adata_ids = AdataIdsSfaira()

    def load_into_memory(self):
        """
        Loads backed objects from DistributedStoreBase into single adata object in memory in .data slot.
        :return:
        """
        if isinstance(self.data, DistributedStoreSingleFeatureSpace):
            adata = None
            for k, v in self.data.indices.items():
                x = self.data.adata_by_key[k][v, :].to_memory()
                x.obs["dataset_id"] = k
                if adata is None:
                    adata = x
                else:
                    adata = adata.concatenate(x)
            self.data = adata

    @property
    @abc.abstractmethod
    def topology_dict(self) -> dict:
        pass

    @abc.abstractmethod
    def init_estim(self):
        pass

    @abc.abstractmethod
    def save_eval(self, fn: str, **kwargs):
        pass

    @abc.abstractmethod
    def _save_specific(self, fn: str, **kwargs):
        pass

    def save(
            self,
            fn: str,
            model: bool = True,
            specific: bool = True
    ):
        """
        Save weights and summary statistics.
        """
        assert self.estimator is not None, "initialize estimator first"
        if model:
            self.estimator.model.training_model.save_weights(fn + "_weights.h5")
        self.save_eval(fn=fn)
        with open(fn + '_history.pickle', 'wb') as f:
            pickle.dump(obj=self.estimator.history, file=f)
        with open(fn + '_hyperparam.pickle', 'wb') as f:
            pickle.dump(obj=self.estimator.train_hyperparam, file=f)
        with open(fn + '_model_hyperparam.pickle', 'wb') as f:
            pickle.dump(obj=self.estimator.model.hyperparam, file=f)
        if specific:
            self._save_specific(fn=fn)

    def n_counts(self, idx):
        return np.asarray(
            self.estimator.data.X[np.sort(idx), :].sum(axis=1)[np.argsort(idx)]
        ).flatten()


class TrainModelEmbedding(TrainModel):

    estimator: EstimatorKerasEmbedding

    def __init__(
            self,
            model_path: str,
            data: Union[str, anndata.AnnData, Universe, DistributedStoreSingleFeatureSpace],
    ):
        super(TrainModelEmbedding, self).__init__(data=data)
        self.estimator = None
        self.model_dir = model_path

    @property
    def topology_dict(self) -> dict:
        topology_dict = self.zoo.topology_container.topology
        return topology_dict

    def init_estim(
            self,
            override_hyperpar: Union[dict, None] = None
    ):
        assert self.zoo.model_id is not None, "choose model in zoo first"
        self.estimator = EstimatorKerasEmbedding(
            data=self.data,
            model_dir=self.model_dir,
            model_id=self.zoo.model_id,
            model_topology=self.zoo.topology_container
        )
        self.estimator.init_model(override_hyperpar=override_hyperpar)
        print(f"TRAINER: initialised model with {self.estimator.topology_container.n_var} features.")

    def save_eval(self, fn: str, **kwargs):
        evaluation_train = self.estimator.evaluate_any(idx=self.estimator.idx_train)
        evaluation_val = self.estimator.evaluate_any(idx=self.estimator.idx_eval)
        evaluation_test = self.estimator.evaluate_any(idx=self.estimator.idx_test)
        evaluation_all = self.estimator.evaluate_any(idx=None)
        evaluation = {
            'test': evaluation_test,
            'val': evaluation_val,
            'train': evaluation_train,
            'all': evaluation_all
        }
        with open(fn + '_evaluation.pickle', 'wb') as f:
            pickle.dump(obj=evaluation, file=f)

    def _save_specific(self, fn: str, **kwargs):
        """
        Save embedding prediction:

        :param fn:
        :return:
        """
        embedding = self.estimator.predict_embedding()
        df_summary = self.estimator.obs_test
        df_summary = df_summary[[k for k in df_summary.columns if k in self._adata_ids.obs_keys]]
        df_summary["ncounts"] = self.n_counts(idx=self.estimator.idx_test)
        np.save(file=fn + "_embedding", arr=embedding)
        df_summary.to_csv(fn + "_covar.csv")
        with open(fn + "_topology.pickle", "wb") as f:
            pickle.dump(obj=self.topology_dict, file=f)


class TrainModelCelltype(TrainModel):

    estimator: EstimatorKerasCelltype

    def __init__(
            self,
            model_path: str,
            data: Union[str, anndata.AnnData, Universe, DistributedStoreSingleFeatureSpace],
            fn_target_universe: str,
    ):
        super(TrainModelCelltype, self).__init__(data=data)
        self.estimator = None
        self.model_dir = model_path
        self.fn_target_universe = fn_target_universe

    @property
    def topology_dict(self) -> dict:
        topology_dict = self.zoo.topology_container.topology
        # Load target universe leaves into topology dict:
        topology_dict["output"]["targets"] = self.estimator.celltype_universe.onto_cl.leaves
        return topology_dict

    def init_estim(
            self,
            override_hyperpar: Union[dict, None] = None
    ):
        assert self.zoo.model_id is not None, "choose model in zoo first"
        self.estimator = EstimatorKerasCelltype(
            data=self.data,
            model_dir=self.model_dir,
            model_id=self.zoo.model_id,
            model_topology=self.zoo.topology_container
        )
        self.estimator.celltype_universe.load_target_universe(self.fn_target_universe)
        self.estimator.init_model(override_hyperpar=override_hyperpar)
        print(f"TRAINER: initialised model with {self.estimator.topology_container.n_var} features and "
              f"{self.estimator.ntypes} labels: \n{self.estimator.ontology_names}.")

    def save_eval(self, fn: str, eval_weighted: bool = False, **kwargs):
        evaluation = {
            'train': self.estimator.evaluate_any(idx=self.estimator.idx_train, weighted=False),
            'val': self.estimator.evaluate_any(idx=self.estimator.idx_eval, weighted=False),
            'test': self.estimator.evaluate_any(idx=self.estimator.idx_test, weighted=False),
            'all': self.estimator.evaluate_any(idx=None, weighted=False)
        }
        with open(fn + '_evaluation.pickle', 'wb') as f:
            pickle.dump(obj=evaluation, file=f)
        if eval_weighted:
            evaluation_weighted = {
                'train': self.estimator.evaluate_any(idx=self.estimator.idx_train, weighted=True),
                'val': self.estimator.evaluate_any(idx=self.estimator.idx_eval, weighted=True),
                'test': self.estimator.evaluate_any(idx=self.estimator.idx_test, weighted=True),
                'all': self.estimator.evaluate_any(idx=None, weighted=True)
            }
            with open(fn + '_evaluation_weighted.pickle', 'wb') as f:
                pickle.dump(obj=evaluation_weighted, file=f)

    def _save_specific(self, fn: str, **kwargs):
        """
        Save true and predicted labels on test set:

        :param fn:
        :return:
        """
        obs = self.estimator.data.obs
        ytrue = self.estimator.ytrue()
        yhat = self.estimator.predict()
        df_summary = self.estimator.obs_test
        df_summary = df_summary[[k for k in df_summary.columns if k in self._adata_ids.obs_keys]]
        df_summary["ncounts"] = self.n_counts(idx=self.estimator.idx_test)
        np.save(file=fn + "_ytrue", arr=ytrue)
        np.save(file=fn + "_yhat", arr=yhat)
        df_summary.to_csv(fn + "_covar.csv")
        with open(fn + '_ontology_names.pickle', 'wb') as f:
            pickle.dump(obj=self.estimator.ontology_names, file=f)
        with open(fn + '_ontology_ids.pickle', 'wb') as f:
            pickle.dump(obj=self.estimator.ontology_ids, file=f)
        with open(fn + "_topology.pickle", "wb") as f:
            pickle.dump(obj=self.topology_dict, file=f)

        cell_counts = obs['cell_type'].value_counts().to_dict()
        with open(fn + '_celltypes_valuecounts_wholedata.pickle', 'wb') as f:
            pickle.dump(obj=[cell_counts], file=f)
