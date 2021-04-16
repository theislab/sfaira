import abc
import anndata
import numpy as np
import pandas as pd
import pickle
from typing import Union

from sfaira.data import Universe
from sfaira.estimators import EstimatorKerasCelltype, EstimatorKerasEmbedding
from sfaira.interface import ModelZooEmbedding, ModelZooCelltype


class TrainModel:

    def __init__(
            self,
            config_path: str,
            data_path: str,
            meta_path: str,
            cache_path: str,
    ):
        # Check if handling backed anndata or base path to directory of raw files:
        if data_path.split(".")[-1] == "h5ad":
            self.data = anndata.read(data_path, backed='r')
            if len(self.data.obs.columns) == 0:
                fn_backed_obs = ".".join(data_path.split(".")[:-1]) + "_obs.csv"
                self.data.obs = pd.read_csv(fn_backed_obs)
        else:
            dataset = Universe(data_path=data_path, meta_path=meta_path, cache_path=cache_path)
            dataset.load_config(config_path)
            self.set_data(dataset)

    @abc.abstractmethod
    def set_data(self, dataset):
        pass

    @abc.abstractmethod
    def init_estim(self):
        pass

    @property
    def adata(self):
        """
        Get adata object depending on whether backed or a property of a container class.

        :return:
        """
        if self.data is None:
            raise ValueError("self.data not set yet")
        elif isinstance(self.data, anndata.AnnData):
            return self.data
        else:
            raise ValueError(f"self.data type not recognized: {type(self.data)}")

    @abc.abstractmethod
    def _save_specific(
            self,
            fn: str
    ):
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


class TrainModelEmbedding(TrainModel):

    def __init__(
            self,
            config_path: str,
            data_path: str,
            meta_path: str,
            cache_path: str,
            model_path: str,
    ):
        super(TrainModelEmbedding, self).__init__(config_path=config_path, data_path=data_path, meta_path=meta_path, cache_path=cache_path)
        self.zoo = ModelZooEmbedding(model_lookuptable=None)
        self.estimator = None
        self.model_dir = model_path

    def set_data(self, dataset):
        dataset.load(match_to_reference=True)
        self.data = dataset.adata

    def init_estim(
            self,
            override_hyperpar: Union[dict, None] = None
    ):
        assert self.zoo.model_id is not None, "choose model in zoo first"
        self.estimator = EstimatorKerasEmbedding(
            data=self.adata,
            model_dir=self.model_dir,
            model_id=self.zoo.model_id,
            organism=self.zoo.organism,
            organ=self.zoo.organ,
            model_type=self.zoo.model_type,
            model_topology=self.zoo.model_topology
        )
        self.estimator.init_model(override_hyperpar=override_hyperpar)

    def save_eval(
            self,
            fn: str
    ):
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

    def _save_specific(
            self,
            fn: str
    ):
        """
        Save embedding prediction:

        :param fn:
        :return:
        """
        embedding = self.estimator.predict_embedding()
        df_summary = self.estimator.obs_test[
            ["dataset", "cell_ontology_class", "state_exact", "author", "year", "assay_sc",
             "assay_differentiation", "assay_type_differentiation", "cell_line", "sample_source"]
        ]
        df_summary["ncounts"] = np.asarray(
            self.estimator.data.X[np.sort(self.estimator.idx_test), :].sum(axis=1)[np.argsort(self.estimator.idx_test)]
        ).flatten()
        np.save(file=fn + "_embedding", arr=embedding)
        df_summary.to_csv(fn + "_covar.csv")


class TrainModelCelltype(TrainModel):

    def __init__(
            self,
            config_path: str,
            data_path: str,
            meta_path: str,
            cache_path: str,
            model_path: str,
    ):
        super(TrainModelCelltype, self).__init__(config_path=config_path, data_path=data_path, meta_path=meta_path, cache_path=cache_path)
        self.zoo = ModelZooCelltype(model_lookuptable=None)
        self.estimator = None
        self.model_dir = model_path

    def set_data(self, dataset):
        dataset.subset("annotated", True)
        dataset.load(match_to_reference=True)
        self.data = dataset.adata

    def init_estim(
            self,
            override_hyperpar: Union[dict, None] = None
    ):
        assert self.zoo.model_id is not None, "choose model in zoo first"
        self.estimator = EstimatorKerasCelltype(
            data=self.adata,
            model_dir=self.model_dir,
            model_id=self.zoo.model_id,
            organism=self.zoo.organism,
            organ=self.zoo.organ,
            model_type=self.zoo.model_type,
            model_topology=self.zoo.model_topology
        )
        self.estimator.init_model(override_hyperpar=override_hyperpar)

    def save_eval(
            self,
            fn: str
    ):
        evaluation = {
            'train': self.estimator.evaluate_any(idx=self.estimator.idx_train, weighted=False),
            'val': self.estimator.evaluate_any(idx=self.estimator.idx_eval, weighted=False),
            'test': self.estimator.evaluate_any(idx=self.estimator.idx_test, weighted=False),
            'all': self.estimator.evaluate_any(idx=None, weighted=False)
        }
        with open(fn + '_evaluation.pickle', 'wb') as f:
            pickle.dump(obj=evaluation, file=f)
        evaluation_weighted = {
            'train': self.estimator.evaluate_any(idx=self.estimator.idx_train, weighted=True),
            'val': self.estimator.evaluate_any(idx=self.estimator.idx_eval, weighted=True),
            'test': self.estimator.evaluate_any(idx=self.estimator.idx_test, weighted=True),
            'all': self.estimator.evaluate_any(idx=None, weighted=True)
        }
        with open(fn + '_evaluation_weighted.pickle', 'wb') as f:
            pickle.dump(obj=evaluation_weighted, file=f)

    def _save_specific(
            self,
            fn: str
    ):
        """
        Save true and predicted labels on test set:

        :param fn:
        :return:
        """
        ytrue = self.estimator.ytrue()
        yhat = self.estimator.predict()
        df_summary = self.estimator.obs_test[
            ["dataset", "cell_ontology_class", "state_exact", "author", "year", "assay_sc",
             "assay_differentiation", "assay_type_differentiation", "cell_line", "sample_source"]
        ]
        df_summary["ncounts"] = np.asarray(self.estimator.data.X[self.estimator.idx_test, :].sum(axis=1)).flatten()
        np.save(file=fn + "_ytrue", arr=ytrue)
        np.save(file=fn + "_yhat", arr=yhat)
        df_summary.to_csv(fn + "_covar.csv")
        with open(fn + '_ontology_names.pickle', 'wb') as f:
            pickle.dump(obj=self.estimator.ids, file=f)

        cell_counts = self.data.obs_concat(keys=['cell_ontology_class'])['cell_ontology_class'].value_counts().to_dict()
        cell_counts_leaf = cell_counts.copy()
        for k in cell_counts.keys():
            if k not in self.estimator.ids:
                if k not in self.estimator.celltypes_version.ontology.node_ids:
                    raise(ValueError(f"Celltype '{k}' not found in celltype universe"))
                for leaf in self.estimator.celltypes_version.ontology.node_ids:
                    if leaf not in cell_counts_leaf.keys():
                        cell_counts_leaf[leaf] = 0
                    cell_counts_leaf[leaf] += 1 / len(self.estimator.celltypes_version.ontology.node_ids)
                del cell_counts_leaf[k]
        with open(fn + '_celltypes_valuecounts_wholedata.pickle', 'wb') as f:
            pickle.dump(obj=[cell_counts, cell_counts_leaf], file=f)
