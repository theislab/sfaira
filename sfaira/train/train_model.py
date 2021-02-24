import abc
import anndata
import numpy as np
import pandas as pd
import pickle
from typing import Union

from sfaira.data import DatasetGroup, DatasetSuperGroup
from sfaira.estimators import EstimatorKeras, EstimatorKerasCelltype, EstimatorKerasEmbedding
from sfaira.interface import ModelZoo, ModelZooEmbedding, ModelZooCelltype


class TargetZoos:
    """
    Class that provides access to all available dataset human in sfaira.

    Parameters
    ----------
    data_path : str
        Path to the files for this dataset on disk
    meta_path : str
        Path to the meta files for this dataset on disk
    """

    def __init__(self, data_path: Union[str, None], meta_path: Union[str, None] = None, cache_path: Union[str, None] = None):
        if data_path is not None:
            from sfaira.data.dataloaders.anatomical_groups import mouse, human
            self.data_mouse = {
                "bladder": mouse.DatasetGroupBladder(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "brain": mouse.DatasetGroupBrain(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "diaphragm": mouse.DatasetGroupDiaphragm(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "adipose": mouse.DatasetGroupAdipose(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "heart": mouse.DatasetGroupHeart(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "kidney": mouse.DatasetGroupKidney(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "colon": mouse.DatasetGroupColon(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "muscle": mouse.DatasetGroupMuscle(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "liver": mouse.DatasetGroupLiver(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "lung": mouse.DatasetGroupLung(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "mammarygland": mouse.DatasetGroupMammaryGland(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "bone": mouse.DatasetGroupBone(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "femalegonad": mouse.DatasetGroupFemalegonad(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "pancreas": mouse.DatasetGroupPancreas(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "blood": mouse.DatasetGroupBlood(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "placenta": mouse.DatasetGroupPlacenta(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "prostate": mouse.DatasetGroupProstate(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "rib": mouse.DatasetGroupRib(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "skin": mouse.DatasetGroupSkin(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "ileum": mouse.DatasetGroupIleum(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "spleen": mouse.DatasetGroupSpleen(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "stomach": mouse.DatasetGroupStomach(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "malegonad": mouse.DatasetGroupMalegonad(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "thymus": mouse.DatasetGroupThymus(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "tongue": mouse.DatasetGroupTongue(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "trachea": mouse.DatasetGroupTrachea(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                "uterus": mouse.DatasetGroupUterus(data_path=data_path, cache_path=cache_path),
            }
            self.data_human = {
                'adipose': human.DatasetGroupAdipose(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'adrenalgland': human.DatasetGroupAdrenalgland(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'artery': human.DatasetGroupArtery(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'bladder': human.DatasetGroupBladder(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'blood': human.DatasetGroupBlood(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'bone': human.DatasetGroupBone(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'brain': human.DatasetGroupBrain(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'calvaria': human.DatasetGroupCalvaria(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'cervix': human.DatasetGroupCervix(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'chorionicvillus': human.DatasetGroupChorionicvillus(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'colon': human.DatasetGroupColon(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'duodenum': human.DatasetGroupDuodenum(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'epityphlon': human.DatasetGroupEpityphlon(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'esophagus': human.DatasetGroupEsophagus(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'eye': human.DatasetGroupEye(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'fallopiantube': human.DatasetGroupFallopiantube(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'femalegonad': human.DatasetGroupFemalegonad(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'gallbladder': human.DatasetGroupGallbladder(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'heart': human.DatasetGroupHeart(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'hesc': human.DatasetGroupHesc(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'ileum': human.DatasetGroupIleum(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'jejunum': human.DatasetGroupJejunum(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'kidney': human.DatasetGroupKidney(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'liver': human.DatasetGroupLiver(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'lung': human.DatasetGroupLung(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'malegonad': human.DatasetGroupMalegonad(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'muscle': human.DatasetGroupMuscle(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'omentum': human.DatasetGroupOmentum(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'pancreas': human.DatasetGroupPancreas(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'placenta': human.DatasetGroupPlacenta(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'pleura': human.DatasetGroupPleura(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'prostate': human.DatasetGroupProstate(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'rectum': human.DatasetGroupRectum(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'rib': human.DatasetGroupRib(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'skin': human.DatasetGroupSkin(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'spinalcord': human.DatasetGroupSpinalcord(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'spleen': human.DatasetGroupSpleen(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'stomach': human.DatasetGroupStomach(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'thymus': human.DatasetGroupThymus(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'thyroid': human.DatasetGroupThyroid(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'trachea': human.DatasetGroupTrachea(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'ureter': human.DatasetGroupUreter(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
                'uterus': human.DatasetGroupUterus(data_path=data_path, meta_path=meta_path, cache_path=cache_path),
            }

        else:
            self.data_human = None
            self.data_mouse = None

    def write_celltypes_tocsv_mouse(self, fn: str):
        for x in self.data_mouse.keys():
            ds = self.data_mouse[x]
            self._write_celltypes_tocsv(fn, x, ds)

    def write_celltypes_tocsv_human(self, fn: str):
        for x in self.data_human.keys():
            ds = self.data_human[x]
            self._write_celltypes_tocsv(fn, x, ds)

    def _write_celltypes_tocsv(self, fn: str, x: str, ds: DatasetGroup):
        ds.load(annotated_only=True, remove_gene_version=False, match_to_reference=None)
        if len(ds.adata_ls) > 0:
            obs = ds.obs_concat(keys=["cell_ontology_class", "cell_ontology_id"])
            obs.index = range(0, obs.shape[0])
            strids = []
            listids = []
            for i in obs.index:
                if type(obs.loc[i]['cell_ontology_class']) != list:
                    strids.append(i)
                else:
                    listids.append(i)
            remaining = []
            for _, l in obs.iloc[listids].iterrows():
                if type(l['cell_ontology_id']) == list:
                    if not len(l['cell_ontology_class']) == len(l['cell_ontology_id']):
                        raise ValueError(
                            "Number of cell type labels and cell type ontologies for this cell do not match")
                    for i in range(len(l['cell_ontology_class'])):
                        remaining.append({
                            "cell_ontology_class": l['cell_ontology_class'][i],
                            "cell_ontology_id": l['cell_ontology_id'][i]
                        })
                else:
                    for i in range(len(l['cell_ontology_class'])):
                        remaining.append({
                            "cell_ontology_class": l['cell_ontology_class'][i],
                            "cell_ontology_id": None
                        })
            obs = obs.loc[strids]
            for i in remaining:
                obs = obs.append(i, ignore_index=True)
            obs = obs.drop_duplicates()
            obs = obs.sort_values(by="cell_ontology_class")
            obs.index = range(0, obs.shape[0])
            obs.to_csv(fn + x + ".csv")


class TrainModel(TargetZoos):

    estimator: Union[None, EstimatorKeras]
    zoo: Union[None, ModelZoo]
    model_dir: str
    data: Union[DatasetGroup, DatasetSuperGroup, anndata.AnnData, str, None]

    def __init__(self, data_path: str, meta_path: str):
        # Check if handling backed anndata or base path to directory of raw files:
        if data_path.split(".")[-1] == "h5ad":
            self.data = anndata.read(data_path, backed='r')
            if len(self.data.obs.columns) == 0:
                fn_backed_obs = ".".join(data_path.split(".")[:-1]) + "_obs.csv"
                self.data.obs = pd.read_csv(fn_backed_obs)
        else:
            super(TrainModel, self).__init__(data_path=data_path, meta_path=meta_path)
            self.data = None

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
        elif isinstance(self.data, DatasetGroup) or isinstance(self.data, DatasetSuperGroup):
            return self.data.adata
        else:
            raise ValueError("self.data type not recognized: %s " % type(self.data))

    def mouse_target(self, organ: str):
        self.set_data(data_group=self.data_mouse[organ])

    def human_target(self, organ: str):
        self.set_data(data_group=self.data_human[organ])

    def set_data(
            self,
            data_group: Union[DatasetGroup, DatasetSuperGroup]
    ):
        """
        Set input data group.
        :return:
        """
        self.data = data_group

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
            data_path: str,
            meta_path: str,
            model_path: str
    ):
        super(TrainModelEmbedding, self).__init__(data_path=data_path, meta_path=meta_path)
        self.zoo = ModelZooEmbedding(model_lookuptable=None)
        self.estimator = None
        self.model_dir = model_path

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
            ["dataset", "cell_ontology_class", "state_exact", "lab", "year", "subtissue", "protocol"]
        ]
        df_summary["ncounts"] = np.asarray(
            self.estimator.data.X[np.sort(self.estimator.idx_test), :].sum(axis=1)[np.argsort(self.estimator.idx_test)]
        ).flatten()
        np.save(file=fn + "_embedding", arr=embedding)
        df_summary.to_csv(fn + "_covar.csv")


class TrainModelCelltype(TrainModel):

    def __init__(
            self,
            data_path: str,
            meta_path: str,
            model_path: str
    ):
        super(TrainModelCelltype, self).__init__(data_path=data_path, meta_path=meta_path)
        self.zoo = ModelZooCelltype(model_lookuptable=None)
        self.estimator = None
        self.model_dir = model_path

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
            ["dataset", "cell_ontology_class", "state_exact", "lab", "year", "subtissue", "protocol"]
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
