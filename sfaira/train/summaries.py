import numpy as np
import pandas
import pickle
import json
import shutil
import warnings
from typing import Union, List
import os
from .train_model import TargetZoos
from .external import SPECIES_DICT

from .external import EstimatorKerasEmbedding


def _tp(yhat, ytrue):
    """
    Class wise true positive count.

    :param yhat:
    :param ytrue:
    :return:
    """
    yhat_true = np.asarray(yhat == np.max(yhat, axis=1, keepdims=True), dtype="float32")
    return np.sum(yhat_true * ytrue, axis=0)


def _fp(yhat, ytrue):
    """
    Class wise false positive count.

    :param yhat:
    :param ytrue:
    :return:
    """
    yhat_true = np.asarray(yhat == np.max(yhat, axis=1, keepdims=True), dtype="float32")
    return np.sum(yhat_true * (1. - ytrue), axis=0)


def _tn(yhat, ytrue):
    """
    Class wise true negative count.

    :param yhat:
    :param ytrue:
    :return:
    """
    yhat_true = np.asarray(yhat < np.max(yhat, axis=1, keepdims=True), dtype="float32")
    return np.sum(yhat_true * (1. - ytrue), axis=0)


def _fn(yhat, ytrue):
    """
    Class wise false negative count.

    :param yhat:
    :param ytrue:
    :return:
    """
    yhat_true = np.asarray(yhat < np.max(yhat, axis=1, keepdims=True), dtype="float32")
    return np.sum(yhat_true * ytrue, axis=0)


def accuracy(yhat, ytrue):
    """
    Class wise accuracy.

    :param yhat:
    :param ytrue:
    :return:
    """
    return (_tp(yhat, ytrue) + _tn(yhat, ytrue)) / yhat.shape[0]


def f1(yhat, ytrue):
    """
    Class wise F1.

    :param yhat:
    :param ytrue:
    :return:
    """
    precision = _tp(yhat, ytrue) / (_tp(yhat, ytrue) + _fp(yhat, ytrue))
    recall = _tp(yhat, ytrue) / (_tp(yhat, ytrue) + _fn(yhat, ytrue))
    return 2 * 1 / (1 / precision + 1 / recall)


def auc_roc(yhat, ytrue):
    """
    Class wise AUC ROC.

    :param yhat:
    :param ytrue:
    :return:
    """
    import sklearn

    auc_roc = np.array([
        sklearn.metrics.roc_auc_score(y_true=ytrue[:, i], y_score=yhat[:, i])
        for i in range(ytrue.shape[0])
    ])
    return auc_roc


class GridsearchContainer:

    histories: Union[None, dict]
    evals: Union[None, dict]
    run_ids: Union[None, list]
    gs_keys: Union[None, dict]
    summary_tab: Union[None, pandas.DataFrame]
    cv: bool
    source_path: str
    model_id_len: Union[None, int]

    def __init__(
            self,
            source_path: str,
            cv: bool
    ):
        self.histories = None
        self.evals = None
        self.run_ids = None
        self.gs_keys = None
        self.cv = cv
        self.source_path = source_path
        self.summary_tab = None

    def load_gs(
            self,
            gs_ids: List[str]
    ):
        """
        Loads all relevant data of a grid search.

        :param gs_ids:
        :return:
        """
        res_dirs = [self.source_path + x + "/results/" for x in gs_ids]
        run_ids = [
            np.sort(np.unique([
                x.split("_history.pickle")[0]
                for x in os.listdir(indir)
                if "_history.pickle" in x
            ]))
            for i, indir in enumerate(res_dirs)
        ]
        histories = {}
        evals = {}
        hyperpars = {}
        model_hyperpars = {}
        run_ids_proc = []
        gs_keys = []
        for i, indir in enumerate(res_dirs):
            for x in run_ids[i]:
                fn_history = indir + x + "_history.pickle"
                if os.path.isfile(fn_history):
                    with open(fn_history, 'rb') as f:
                        histories[x] = pickle.load(f)
                else:
                    print("file %s not found" % (x + "_history.pickle"))
                fn_eval = indir + x + "_evaluation.pickle"
                if os.path.isfile(fn_eval):
                    with open(fn_eval, 'rb') as f:
                        evals[x] = pickle.load(f)
                else:
                    print("file %s not found" % (x + "_evaluation.pickle"))
                fn_hp = indir + x + "_hyperparam.pickle"
                if os.path.isfile(fn_hp):
                    with open(fn_hp, 'rb') as f:
                        hyperpars[x] = pickle.load(f)
                else:
                    print("file %s not found" % (x + "_hyperparam.pickle"))
                fn_mhp = indir + x + "_model_hyperparam.pickle"
                if os.path.isfile(fn_mhp):
                    with open(fn_mhp, 'rb') as f:
                        model_hyperpars[x] = pickle.load(f)
                else:
                    pass
                    #TODO add: print("file %s not found" % (x + "_model_hyperparam.pickle"))
                run_ids_proc.append(x)
                gs_keys.append(indir.split("/")[-3])

        self.run_ids = run_ids_proc
        self.gs_keys = dict(zip(run_ids_proc, gs_keys))
        self.evals = evals
        self.hyperpars = hyperpars
        self.model_hyperpars = model_hyperpars
        self.histories = histories

    def load_y(
            self,
            hat_or_true: str,
            run_id: str
    ):
        fn = self.source_path + self.gs_keys[run_id] + "/results/" + run_id + f"_y{hat_or_true}.npy"
        return np.load(fn)

    def best_model_by_partition(
            self,
            partition_select: str,
            metric_select: str,
            cv_mode: str = "mean",
            subset: dict = {},
            return_run_only: bool = False,
            grouping: list = ["organ", "model_type"]
    ):
        """

        :param partition_select:
        :param metric_select:
        :param cv_mode:
        :param subset:
        :param return_run_only:
        :param grouping:
        :return:
        """
        model_ids = []
        run_ids = []
        for id, df in self.summary_tab.groupby(grouping):
            if df.shape[0] > 0:
                model_id_temp, run_id_temp, cv_id_temp = self.get_best_model_ids(
                    tab=df,
                    partition_select=partition_select,
                    metric_select=metric_select,
                    subset=subset,
                    cv_mode=cv_mode
                )
                model_ids.append(model_id_temp)
                run_ids.append(run_id_temp)
        if return_run_only:
            return self.summary_tab.loc[[x in run_ids for x in self.summary_tab["run"].values], :]
        else:
            return self.summary_tab.loc[[x in model_ids for x in self.summary_tab["model_gs_id"].values], :]

    def get_best_model_ids(
            self,
            tab,
            metric_select: str,
            partition_select: str,
            subset: dict = {},
            cv_mode: str = "mean"
    ):
        """

        :param tab:
        :param metric_select:
        :param partition_select:
        :param subset:
        :param cv_mode:
        :return:
        """
        for k, v in subset.items():
            tab = tab.loc[tab[k].values == v, :]

        if metric_select.endswith('accuracy') \
                or metric_select.endswith('acc_agg') \
                or metric_select.endswith('f1') \
                or metric_select.endswith('tpr'):
            ascending = False
            if cv_mode == "min":
                raise Warning("selected cv_mode min with metric_id %s, likely not intended" % metric_select)
        elif metric_select.endswith('loss') \
                or metric_select.endswith('mse') \
                or metric_select.endswith('negll') \
                or metric_select.endswith('custom_cce_agg') \
                or metric_select.endswith('fpr'):
            ascending = True
            if cv_mode == "max":
                raise Warning("selected cv_mode max with metric_id %s, likely not intended" % metric_select)
        else:
            raise ValueError("measure %s not recognized" % metric_select)

        if partition_select not in ["test", "val", "train"]:
            raise ValueError("partition %s not recognised" % partition_select)

        metric_select = partition_select + "_" + metric_select

        if cv_mode.lower() == "mean":
            best_model = tab.groupby("run", as_index=False)[metric_select].mean().\
                sort_values([metric_select], ascending=ascending)
        elif cv_mode.lower() == "median":
            best_model = tab.groupby("run", as_index=False)[metric_select].median().\
                sort_values([metric_select], ascending=ascending)
        elif cv_mode.lower() == "max":
            best_model = tab.groupby("run", as_index=False)[metric_select].max().\
                sort_values([metric_select], ascending=ascending)
        elif cv_mode.lower() == "min":
            best_model = tab.groupby("run", as_index=False)[metric_select].min().\
                sort_values([metric_select], ascending=ascending)
        else:
            raise ValueError("cv_mode %s not recognized" % cv_mode)

        best_run_id = best_model['run'].values[0] if best_model.shape[0] > 0 else None

        best_cv = tab[tab["run"] == best_run_id].\
            sort_values([metric_select], ascending=ascending)['cv'].values[0] if best_run_id is not None \
            else None

        best_model_id = tab[tab["run"] == best_run_id]. \
            sort_values([metric_select], ascending=ascending)['model_gs_id'].values[0] if best_run_id is not None \
            else None

        return best_model_id, best_run_id, best_cv

    @property
    def cv_keys(self) -> List[str]:
        """
        Returns keys of cross-validation used in dictionaries in this class.

        :return: list of string keys
        """
        return np.unique(self.summary_tab["cv"].values).tolist()

    def save_best_weight(
            self,
            path: str,
            partition: str = "val",
            metric: str = "loss",
            subset: dict = {}
    ):
        """
        Copies weight file from best hyperparameter setting from grid search directory to zoo directory with cleaned
        file name.

        :param path: Target file to save to. This is intended to be the zoo directory ready for upload.
        :param partition:
        :param metric:
        :param subset:
        :return:
        """
        assert not self.cv, "not implemented for CV yet"
        model_id, _, _ = self.get_best_model_ids(
            tab=self.summary_tab,
            partition_select=partition,
            metric_select=metric,
            cv_mode="mean",
            subset=subset,
        )
        shutil.copyfile(
            self.source_path + self.gs_keys[model_id] + "/results/" + model_id + "_weights.h5",
            path + model_id + "_weights.h5"
        )

    def plot_completions(
            self,
            groupby=["depth", "width", "lr", "dropout", "l1", "l2"],
            height_fig=7,
            width_fig=7
    ):
        """
        Plot number of completed grid search points by category.

        :param groupby:
        :param height_fig:
        :param width_fig:
        :return:
        """
        import matplotlib.pyplot as plt
        import seaborn as sns

        if self.summary_tab is None:
            self.create_summary_tab()
        sns_tab = self.summary_tab.copy()

        # Build figure.
        organs = np.unique(sns_tab["organ"].values)
        model_types = np.unique(sns_tab["model_type"].values)
        hm = np.zeros((len(organs), len(model_types)))
        for i, m in enumerate(model_types):
            for j, o in enumerate(organs):
                n_by_gridpoint = sns_tab.loc[
                    np.logical_and(
                        sns_tab["model_type"].values == m,
                        sns_tab["organ"].values == o
                    ), :
                ].groupby(groupby).size().values
                # Assume that largest number of successful completions is maximum (all completed:
                hm[j, i] = np.sum(n_by_gridpoint == np.max(n_by_gridpoint)) if len(n_by_gridpoint) > 0 else 0
        sns_data_heatmap = pandas.DataFrame(
            hm, index=organs, columns=model_types
        )
        fig, axs = plt.subplots(1, 1, figsize=(height_fig, width_fig))
        with sns.axes_style("dark"):
            axs = sns.heatmap(
                sns_data_heatmap,
                annot=True, fmt=".2f",
                ax=axs,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'n'}
            )
        return fig, axs

    def plot_best_model_by_hyperparam(
            self,
            metric_select: str,
            param_hue='lr',
            partition_select: str = "val",
            partition_show: str = "test",
            subset: dict = {},
            param_x=['lr', 'depth', 'width', 'dropout', 'l1', 'l2'],
            show_swarm: bool = False,
            panel_width: float = 4.,
            panel_height: float = 2.
    ):
        """
        Produces boxplots for all hyperparameters choices by organ.

        :param partition: "train" or "eval" or "test" partition of data.
        :param metric_select: Metric to plot.
        :param param_x: Hyper-parameter for x-axis partition.
        :param param_hue: Hyper-parameter for hue-axis partition.
        :param panel_width:
        :param panel_height:
        :return:
        """
        import seaborn as sns
        import matplotlib.pyplot as plt

        params = [param for param in param_x if len(np.unique(self.summary_tab[param])) > 1 and param != param_hue]

        organs = np.unique(self.summary_tab["organ"].values)
        fig, ax = plt.subplots(
            nrows=len(organs), ncols=len(params),
            figsize=(panel_width * len(params), panel_height * len(organs))
        )
        if len(organs) == 1:
            ax = np.expand_dims(ax, axis=0)
        for j, param in enumerate(params):
            summary_table_param = self.best_model_by_partition(
                partition_select=partition_select,
                metric_select=metric_select,
                cv_mode="mean",
                subset=subset,
                return_run_only=False,
                grouping=["organ", param, param_hue]
            )
            summary_table_param.sort_values([param, param_hue])
            for i, organ in enumerate(organs):
                summary_table = summary_table_param.loc[summary_table_param["organ"].values == organ, :]
                # Plot each metric:
                ycol = partition_show + "_" + metric_select
                if len(organs) == 1 and len(params) == 1:
                    ax = np.array([ax])
                sns.boxplot(
                    x=param, hue=param_hue, y=ycol,
                    data=summary_table, ax=ax[i, j] if len(ax.shape) == 2 else ax[i] if len(ax.shape) == 1 else ax
                )
                if show_swarm:
                    sns.swarmplot(
                        x=param, hue=param_hue, y=ycol,
                        data=summary_table, ax=ax[i, j] if len(ax.shape) == 2 else ax[i] if len(ax.shape) == 1 else ax
                    )
                if j == 0:
                    if len(ax.shape) == 2:
                        ax[i, j].set_ylabel(organ + "\n" + ycol)
                    else:
                        ax[i].set_ylabel(organ + "\n" + ycol)
        return fig, ax

    def plot_training_history(
            self,
            metric_select: str,
            metric_show: str,
            partition_select: str = "val",
            subset: dict = {},
            cv_key: Union[str, None] = None,
            log_loss: bool = False
    ):
        """
        Plot train and validation loss during training and learning rate reduction for each organ

        The partition that is shown in train+val by default because these are the only ones recorded during training.

        :param metric_select: metric to select best model by
        :param metric_show: metric to show as function of training progress, together with loss and learing rate.
        :param partition_select: "train" or "eval" or "test" partition of data to select fit by.
        :param metric_select: Metric to select fit by.
        :param cv_key: Index of cross-validation to plot training history for.
        :param log_loss:
        :return:
        """
        import seaborn as sns
        import matplotlib.pyplot as plt
        panel_width = 5
        panel_height = 3

        organs = np.unique(self.summary_tab["organ"].values)
        fig, ax = plt.subplots(
            nrows=len(organs), ncols=3,
            figsize=(panel_width * 3, panel_height * len(organs))
        )
        if len(organs) == 1:
            ax = np.expand_dims(ax, axis=0)
        for i, organ in enumerate(organs):
            model_gs_id, _, _ = self.get_best_model_ids(
                tab=self.summary_tab,
                partition_select=partition_select,
                metric_select=metric_select,
                cv_mode="mean",
                subset=dict(list(subset.items()) + [("organ", organ)]),
            )
            if cv_key is None:
                sns_data = []
                for run in np.unique(
                        self.summary_tab.loc[self.summary_tab["model_gs_id"].values == model_gs_id, "run"].values
                ).tolist():
                    sns_data_temp = pandas.DataFrame(self.histories[run])
                    sns_data_temp["epoch"] = np.arange(0, sns_data_temp.shape[0])
                    sns_data_temp["cv"] = run.split("_")[-1]
                    sns_data.append(sns_data_temp)
                sns_data = pandas.concat(sns_data, axis=0)
            else:
                cv = cv_key
                sns_data = pandas.DataFrame(self.histories[model_gs_id + "_" + cv])
                sns_data["epoch"] = np.arange(0, sns_data.shape[0])
                sns_data["cv"] = cv

            # loss
            sns_data_loss = pandas.concat([pandas.DataFrame({
                "epoch": sns_data["epoch"].values,
                "cv": sns_data["cv"].values,
                "loss": np.log(sns_data[x].values) if log_loss else sns_data[x].values,
                "partition": x
            }) for i, x in enumerate(["loss", "val_loss"])])
            sns.lineplot(
                x="epoch", y="loss", style="partition", hue="cv",
                data=sns_data_loss, ax=ax[i, 0]
            )
            ax[i, 0].set_ylabel(organ + "\nloss")
            ax[i, 0].legend_.remove()

            # metric
            if metric_show not in sns_data.columns:
                raise ValueError("metric %s not found in %s" % (metric_show, str(sns_data.columns)))
            sns_data_metric = pandas.concat([pandas.DataFrame({
                "epoch": sns_data["epoch"].values,
                "cv": sns_data["cv"].values,
                metric_show: sns_data[metric_show].values,
                "partition": x
            }) for i, x in enumerate([metric_show, "val_" + metric_show])])
            sns.lineplot(
                x="epoch", y=metric_show, style="partition", hue="cv",
                data=sns_data_metric, ax=ax[i, 1]
            )
            ax[i, 1].set_ylabel(organ + "\n" + metric_show)
            ax[i, 1].legend_.remove()

            # lr
            sns_data_lr = pandas.DataFrame({
                "epoch": sns_data["epoch"].values,
                "cv": sns_data["cv"].values,
                "lr": np.log(sns_data["lr"].values) / np.log(10)
            })
            sns.lineplot(
                x="epoch", y="lr", hue="cv",
                data=sns_data_lr, ax=ax[i, 2]
            )
            ax[i, 2].set_ylabel("log10 learning rate")
            ax[i, 2].legend_.remove()
        return fig, ax

    def write_best_hyparam(
            self,
            write_path,
            subset: dict = {},
            partition: str = "test",
            metric: str = "custom_negll",
            cvs: Union[None, List[int]] = None
    ):
        best_model_id = self.get_best_model_ids(
            tab=self.summary_tab,
            subset=subset,
            partition_select=partition,
            metric_select=metric,
        )[0]

        if best_model_id is not None:
            if cvs is None:
                file_path_base = os.path.join(
                    self.source_path,
                    self.gs_keys[best_model_id],
                    'results',
                    best_model_id
                )
            else:
                file_path_base = os.path.join(
                    self.source_path,
                    self.gs_keys[best_model_id + "_cv" + str(cvs[0])],
                    'results',
                    best_model_id + "_cv" + str(cvs[0])
                )

            # Read model hyperparameter
            with open(file_path_base + "_model_hyperparam.pickle", 'rb') as file:
                hyparam_model = pickle.load(file)

            # Read optimizer hyperparameter
            with open(file_path_base + "_hyperparam.pickle", 'rb') as file:
                hyparam_optim = pickle.load(file)

            # Write both hyperparameter dicts
            with open(os.path.join(write_path, best_model_id[:-12] + "_best_hyperparam.txt"), 'w') as file:
                file.write(json.dumps({"model": hyparam_model, "optimizer": hyparam_optim}))
        return


class SummarizeGridsearchCelltype(GridsearchContainer):
    loss_idx: int
    acc_idx: int

    def __init__(
            self,
            source_path: str,
            cv: bool,
            model_id_len: int = 7
    ):
        super(SummarizeGridsearchCelltype, self).__init__(
            source_path=source_path,
            cv=cv
        )
        self.model_id_len = model_id_len

    def load_ontology_names(
            self,
            run_id: str
    ):
        """
        Loads ontology ids from a specific model of a previously loaded grid search.

        :param run_id:
        :return:
        """
        fn = self.source_path + self.gs_keys[run_id] + "/results/" + run_id + "_ontology_names.pickle"
        if not os.path.isfile(fn):
            raise FileNotFoundError(f"file {run_id}_ontology_names.pickle not found")
        with open(fn, 'rb') as f:
            ids = pickle.load(f)
        return ids

    def create_summary_tab(self):
        """
        metrics = list(self.evals.values())[0]['val'].keys()
        hyperpar = list(self.hyperpars.keys())
        model_hyperpar = list(self.hyperpars.keys())
        self.summary_tab = pandas.DataFrame(dict(
            list({
                "cv": [id_i.split("_")[-1] if self.cv else "1" for id_i in self.run_ids],
                "model": ["_".join(id_i.split("_")[:self.model_id_len]) for id_i in self.run_ids],
                "model_type": [id_i.split("_")[3] for id_i in self.run_ids],
                "run": self.run_ids,
            }.items()) +
            list(dict([(hp, [self.hyperpars[id_i][hp] for id_i in self.run_ids]) for hp in hyperpar]).items()) +
            list(dict([(hp, [self.model_hyperpar[id_i][hp] for id_i in self.run_ids]) for hp in model_hyperpar]).items()) +
            list(dict(
                [("train_" + m, [self.evals[id_i]["train"][m] for id_i in self.run_ids]) for m in metrics]).items()) +
            list(dict(
                [("test_" + m, [self.evals[id_i]["test"][m] for id_i in self.run_ids]) for m in metrics]).items()) +
            list(dict([("val_" + m, [self.evals[id_i]["val"][m] for id_i in self.run_ids]) for m in metrics]).items()) +
            list(dict([("all_" + m, [self.evals[id_i]["all"][m] for id_i in self.run_ids]) for m in metrics]).items())
        ))
        :return:
        """
        metrics = list(self.evals.values())[0]['val'].keys()
        self.summary_tab = pandas.DataFrame(dict(
            list({
                "depth": [id_i.split("_")[self.model_id_len + 0] for id_i in self.run_ids],
                "width": [id_i.split("_")[self.model_id_len + 1] for id_i in self.run_ids],
                "lr": [id_i.split("_")[self.model_id_len + 2] for id_i in self.run_ids],
                "dropout": [id_i.split("_")[self.model_id_len + 3] for id_i in self.run_ids],
                "l1": [id_i.split("_")[self.model_id_len + 4] for id_i in self.run_ids],
                "l2": [id_i.split("_")[self.model_id_len + 5] for id_i in self.run_ids],
                "cv": [id_i.split("_")[-1] if self.cv else "cv0" for id_i in self.run_ids],
                "model": ["_".join(id_i.split("_")[:self.model_id_len]) for id_i in self.run_ids],
                "organ": [id_i.split("_")[2] for id_i in self.run_ids],
                "model_type": [
                    "linear" if (id_i.split("_")[3] == "mlp" and id_i.split("_")[5].split(".")[1] == "0")
                    else id_i.split("_")[3]
                    for id_i in self.run_ids
                ],
                "model_gs_id": ["_".join(id_i.split("_")[:(self.model_id_len + 6)]) for id_i in self.run_ids],
                "run": self.run_ids
            }.items()) +
            list(dict([("train_" + m, [self.evals[x]["train"][m] for x in self.run_ids]) for m in metrics]).items()) +
            list(dict([("test_" + m, [self.evals[x]["test"][m] for x in self.run_ids]) for m in metrics]).items()) +
            list(dict([("val_" + m, [self.evals[x]["val"][m] for x in self.run_ids]) for m in metrics]).items()) +
            list(dict([("all_" + m, [self.evals[x]["all"][m] for x in self.run_ids]) for m in metrics]).items())
        ))
        if self.summary_tab.shape[0] == 0:
            raise ValueError("summary_tab was empty")

    def best_model_celltype(
            self,
            subset: dict = {},
            partition: str = "val",
            metric: str = "loss",
            cvs: Union[None, List[int]] = None
    ):
        model_id, _, _ = self.get_best_model_ids(
            tab=self.summary_tab,
            partition_select=partition,
            metric_select=metric,
            cv_mode="mean",
            subset=subset,
        )
        if model_id is not None:
            if cvs is not None:
                fns = [
                    self.source_path + self.gs_keys[model_id + "_cv" + str(x)] + "/results/" + model_id + "_cv" + str(x)
                    for x in cvs
                ]
            else:
                fns = [self.source_path + self.gs_keys[model_id] + "/results/" + model_id]
            covar = [pandas.read_csv(x + "_covar.csv") for x in fns]
            return model_id, covar
        else:
            return None, [None]

    def plot_best(
            self,
            rename_levels=[],
            partition_select: str = "val",
            partition_show: str = "test",
            metric_select: str = "acc",
            metric_show: str = "acc",
            collapse_cv: str = "max",
            vmin=None,
            vmax=None,
            height_fig=7,
            width_fig=7
    ):
        """
        Plot accuracy or other metric heatmap by organ and model type.

        :param rename_levels:
        :param metric: Metric to plot in heatmap.

             - acc
             - f1
        :param collapse_cv: How to collapse values from cross validation into single scalar:

            - mean
            - median
            - max
        :param ylim:
        :param xrot:
        :param height_fig:
        :param width_fig:
        :return:
        """
        import matplotlib.pyplot as plt
        import seaborn as sns

        if self.summary_tab is None:
            self.create_summary_tab()

        # Choose the best over categories based on mean loss in CV.
        # Keep variation across CV.
        sns_tab = self.best_model_by_partition(
            partition_select=partition_select,
            metric_select=metric_select,
            return_run_only=False,
            grouping=["organ", "model_type"]
        )
        for rename_level in rename_levels:
            levels_new = sns_tab[rename_level[0]].values
            levels_new[levels_new == rename_level[1]] = rename_level[2]
            sns_tab[rename_level[0]] = levels_new

        # Build figure.
        organs = np.unique(sns_tab["organ"].values)
        model_types = np.unique(sns_tab["model_type"].values)
        hm = np.zeros((len(organs), len(model_types))) + np.nan
        mask = np.isnan(hm)
        for i, m in enumerate(model_types):
            for j, o in enumerate(organs):
                data_temp = sns_tab.loc[
                    np.logical_and(
                        sns_tab["model_type"].values == m,
                        sns_tab["organ"].values == o
                    ), partition_show + "_" + metric_show
                ]
                if data_temp.shape[0] > 0:
                    if self.cv:
                        if collapse_cv == "mean":
                            hm[j, i] = np.mean(data_temp.values)
                        elif collapse_cv == "median":
                            hm[j, i] = np.median(data_temp.values)
                        elif collapse_cv == "max":
                            hm[j, i] = np.max(data_temp.values)
                        elif collapse_cv == "min":
                            hm[j, i] = np.min(data_temp.values)
                        else:
                            raise ValueError(f"collapse_cv {collapse_cv} not recognized")
                        mask[j, i] = False
                    else:
                        hm[j, i] = data_temp.values[0]
                        mask[j, i] = False
        if vmin is not None:
            hm = np.maximum(hm, np.asarray(vmin))
        if vmax is not None:
            hm = np.minimum(hm, np.asarray(vmin))
        sns_data_heatmap = pandas.DataFrame(
            hm, index=organs, columns=model_types
        )
        fig, axs = plt.subplots(1, 1, figsize=(height_fig, width_fig))
        with sns.axes_style("dark"):
            axs = sns.heatmap(
                sns_data_heatmap, #mask=mask,
                annot=True, fmt=".2f",
                ax=axs, vmin=0, vmax=1,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': partition_show + "_" + metric_show},
                cmap=None
            )
        return fig, axs, sns_data_heatmap

    def plot_best_classwise_heatmap(
            self,
            organ: str,
            organism: str,
            datapath: str,
            celltype_version: str = "0",
            partition_select: str = "val",
            metric_select: str = "custom_cce_agg",
            metric_show: str = "f1",
            collapse_cv: str = "mean",
            min_cells: int = 10,
            height_fig: int = 7,
            width_fig: int = 7
    ):
        """
        Plot evaluation metric heatmap for specified organ by cell classes and model types.

        :param organ: Organ to plot in heatmap.
        :param organism: Species that the gridsearch was run on
        :param datapath: Path to the local sfaira data repository
        :param celltype_version: Version in sfaira celltype universe
        :param partition_select: Based on which partition to select the best model
            - train
            - val
            - test
            - all
        :param metric_select: Based on which metric to select the best model
            - loss
            - accuracy
            - custom_cce_agg
            - acc_agg
            - f1
            - tpr
            - fpr
        :param metric_show: Which classwise metric to plot.
            - accuracy
            - f1
        :param collapse_cv: How to collapse over the single cv runs.
        :param min_cells: Minimum number of cells of a type must be present in the whole dataset for that class to be included in the plot.
        :param height_fig: Figure height.
        :param width_fig: Figure width.
        :return: fig, axs, sns_data_heatmap
        """

        import matplotlib.pyplot as plt
        import seaborn as sns

        if self.summary_tab is None:
            self.create_summary_tab()

        # Choose the best over categories based on mean loss in CV.
        # Keep variation across CV.
        sns_tab = self.best_model_by_partition(
            partition_select=partition_select,
            metric_select=metric_select,
            return_run_only=False,
            grouping=["organ", "model_type"]
        )
        sns_tab = sns_tab[sns_tab['organ'] == organ]

        tz = TargetZoos(path=datapath)
        if organism == "human":
            dataset = tz.data_human[organ]
        elif organism == "mouse":
            dataset = tz.data_mouse[organ]
        else:
            raise(ValueError(f"Supplied organism {organism} not recognised. Should be one of ('mouse', 'human')"))
        dataset.load_all()
        cell_counts = dataset.obs_concat(keys=['cell_ontology_class'])['cell_ontology_class'].value_counts().to_dict()

        celltype_versions = SPECIES_DICT.copy()
        celltype_versions[organism][organ].set_version(celltype_version)
        leafnodes = celltype_versions[organism][organ].ids
        ontology = celltype_versions[organism][organ].ontology[celltype_version]["names"]

        celltypelist = list(cell_counts.keys()).copy()
        for k in celltypelist:
            if k not in leafnodes:
                if k not in ontology.keys():
                    raise(ValueError(f"Celltype '{k}' not found in celltype universe"))
                for leaf in ontology[k]:
                    if leaf not in cell_counts.keys():
                        cell_counts[leaf] = 0
                    cell_counts[leaf] += 1/len(ontology[k])
                del cell_counts[k]

        # Compute class-wise metrics
        vals = []
        for i, run_id in enumerate(sns_tab["run"].values):
            yhat = self.load_y(hat_or_true='hat', run_id=run_id)
            ytrue = self.load_y(hat_or_true='true', run_id=run_id)
            if metric_show == "acc":
                m = accuracy(yhat, ytrue)
            elif metric_show == "f1":
                m = f1(yhat, ytrue)
            else:
                raise ValueError("did not recognize metric_show %s" % metric_show)
            vals.append(m)
        sns_tab[metric_show + "_classwise"] = vals

        # Build figure.
        model_types = sns_tab["model_type"].unique()
        classes = self.load_ontology_names(run_id=sns_tab["run"].values[0])
        if 'unknown' not in classes and 'Unknown' not in classes:
            classes = classes + ['Unknown']
            cell_counts['Unknown'] = 0
        hm = np.zeros((len(classes), len(model_types))) + np.nan
        # mask = np.isnan(hm)
        for i, m in enumerate(model_types):
            data_temp = np.vstack(sns_tab.loc[sns_tab["model_type"].values == m, metric_show + "_classwise"].values)
            if data_temp.shape[0] > 0:
                if self.cv:
                    if collapse_cv == "mean":
                        hm[:, i] = np.nanmean(data_temp, axis=0)
                    elif collapse_cv == "median":
                        hm[:, i] = np.nanmedian(data_temp, axis=0)
                    elif collapse_cv == "max":
                        hm[:, i] = np.nanmax(data_temp, axis=0)
                    elif collapse_cv == "min":
                        hm[:, i] = np.nanmin(data_temp, axis=0)
                    else:
                        raise ValueError(f"collapse_cv {collapse_cv} not recognized")
                else:
                    hm[:, i] = data_temp.values[0]
        n_cells = []
        for c in classes:
            if c in cell_counts.keys():
                n_cells.append(np.round(cell_counts[c]))
            else:
                warnings.warn(f"Celltype {c} from cell ontology now found in {organism} {organ} dataset")
                n_cells.append(np.nan)
        n_cells = np.array(n_cells)[:, None]
        sns_data_heatmap = pandas.DataFrame(
            np.hstack((n_cells, hm)),
            index=classes,
            columns=['Number of cells in whole dataset'] + list(model_types)
        )
        sns_data_heatmap = sns_data_heatmap[sns_data_heatmap['Number of cells in whole dataset'] >= min_cells]
        mask = np.zeros(sns_data_heatmap.shape).astype(bool)
        mask[:, 0] = True
        with sns.axes_style("dark"):
            fig, axs = plt.subplots(1, 1, figsize=(width_fig, height_fig))
            axs = sns.heatmap(
                sns_data_heatmap, mask=mask,
                annot=True, fmt=".2f",
                ax=axs, vmin=0, vmax=1,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': "test_" + metric_show},
                cmap=None
            )
            axs = sns.heatmap(
                data=sns_data_heatmap, mask=~mask,
                annot=True, fmt=".0f",
                ax=axs, alpha=0,
                xticklabels=True, yticklabels=True,
                annot_kws={"color": "black"},
                cbar=False
            )
        return fig, axs, sns_data_heatmap

    def plot_best_classwise_scatter(
            self,
            organ: str,
            organism: str,
            datapath: str,
            celltype_version: str = "0",
            partition_select: str = "val",
            metric_select: str = "custom_cce_agg",
            metric_show: str = "f1",
            collapse_cv: str = "mean",
            min_cells: int = 10,
            height_fig: int = 7,
            width_fig: int = 7,
            annotate_thres_ncells: int = 1000,
            annotate_thres_f1: float = 0.5
    ):
        """
        Plot evaluation metric scatterplot for specified organ by cell classes and model types.

        :param organ: Organ to plot in heatmap.
        :param organism: Species that the gridsearch was run on
        :param datapath: Path to the local sfaira data repository
        :param celltype_version: Version in sfaira celltype universe
        :param partition_select: Based on which partition to select the best model
            - train
            - val
            - test
            - all
        :param metric_select: Based on which metric to select the best model
            - loss
            - accuracy
            - custom_cce_agg
            - acc_agg
            - f1
            - tpr
            - fpr
        :param metric_show: Which classwise metric to plot.
            - accuracy
            - f1
        :param collapse_cv: How to collapse over the single cv runs.
        :param min_cells: Minimum number of cells of a type must be present in the whole dataset for that class to be included in the plot.
        :param height_fig: Figure height.
        :param width_fig: Figure width.
        :param annotate_thres_ncells:
        :param annotate_thres_f1:
        :return: fig, axs, sns_data_scatter
        """

        import matplotlib.pyplot as plt
        import seaborn as sns

        if self.summary_tab is None:
            self.create_summary_tab()

        # Choose the best over categories based on mean loss in CV.
        # Keep variation across CV.
        sns_tab = self.best_model_by_partition(
            partition_select=partition_select,
            metric_select=metric_select,
            return_run_only=False,
            grouping=["organ", "model_type"]
        )
        sns_tab = sns_tab[sns_tab['organ'] == organ]

        tz = TargetZoos(path=datapath)
        if organism == "human":
            dataset = tz.data_human[organ]
        elif organism == "mouse":
            dataset = tz.data_mouse[organ]
        else:
            raise(ValueError(f"Supplied organism {organism} not recognised. Should be one of ('mouse', 'human')"))
        dataset.load_all()
        cell_counts = dataset.obs_concat(keys=['cell_ontology_class'])['cell_ontology_class'].value_counts().to_dict()

        celltype_versions = SPECIES_DICT.copy()
        celltype_versions[organism][organ].set_version(celltype_version)
        leafnodes = celltype_versions[organism][organ].ids
        ontology = celltype_versions[organism][organ].ontology[celltype_version]["names"]

        celltypelist = list(cell_counts.keys()).copy()
        for k in celltypelist:
            if k not in leafnodes:
                if k not in ontology.keys():
                    raise(ValueError(f"Celltype '{k}' not found in celltype universe"))
                for leaf in ontology[k]:
                    if leaf not in cell_counts.keys():
                        cell_counts[leaf] = 0
                    cell_counts[leaf] += 1/len(ontology[k])
                del cell_counts[k]

        # Compute class-wise metrics
        vals = []
        for i, run_id in enumerate(sns_tab["run"].values):
            yhat = self.load_y(hat_or_true='hat', run_id=run_id)
            ytrue = self.load_y(hat_or_true='true', run_id=run_id)
            if metric_show == "acc":
                m = accuracy(yhat, ytrue)
            elif metric_show == "f1":
                m = f1(yhat, ytrue)
            else:
                raise ValueError("did not recognize metric_show %s" % metric_show)
            vals.append(m)
        sns_tab[metric_show + "_classwise"] = vals

        # Build figure.
        model_types = sns_tab["model_type"].unique()
        classes = self.load_ontology_names(run_id=sns_tab["run"].values[0])
        if 'unknown' not in classes and 'Unknown' not in classes:
            classes = classes + ['Unknown']
            cell_counts['Unknown'] = 0
        hm = np.zeros((len(classes), len(model_types))) + np.nan
        # mask = np.isnan(hm)
        for i, m in enumerate(model_types):
            data_temp = np.vstack(sns_tab.loc[sns_tab["model_type"].values == m, metric_show + "_classwise"].values)
            if data_temp.shape[0] > 0:
                if self.cv:
                    if collapse_cv == "mean":
                        hm[:, i] = np.nanmean(data_temp, axis=0)
                    elif collapse_cv == "median":
                        hm[:, i] = np.nanmedian(data_temp, axis=0)
                    elif collapse_cv == "max":
                        hm[:, i] = np.nanmax(data_temp, axis=0)
                    elif collapse_cv == "min":
                        hm[:, i] = np.nanmin(data_temp, axis=0)
                    else:
                        raise ValueError(f"collapse_cv {collapse_cv} not recognized")
                else:
                    hm[:, i] = data_temp.values[0]
        n_cells = []
        for c in classes:
            if c in cell_counts.keys():
                n_cells.append(np.round(cell_counts[c]))
            else:
                warnings.warn(f"Celltype {c} from cell ontology now found in {organism} {organ} dataset")
                n_cells.append(np.nan)
        n_cells = np.array(n_cells)[:, None]
        sns_data_scatter = pandas.DataFrame(
            np.hstack((n_cells, hm)),
            index=classes,
            columns=['Number of cells in whole dataset'] + list(model_types)
        )
        sns_data_scatter = sns_data_scatter[sns_data_scatter['Number of cells in whole dataset'] >= min_cells]
        sns_data_scatter = pandas.melt(sns_data_scatter,
                                       id_vars=['Number of cells in whole dataset'],
                                       value_vars=list(model_types),
                                       var_name='Model type',
                                       value_name='Classwise f1 score',
                                       ignore_index=False
                                       )

        with sns.axes_style("dark"):
            fig, axs = plt.subplots(1, 1, figsize=(width_fig, height_fig))
            axs = sns.scatterplot(x='Number of cells in whole dataset',
                                  y='Classwise f1 score',
                                  style='Model type',
                                  data=sns_data_scatter,
                                  ax=axs
                                  )
            for line in range(0, sns_data_scatter.shape[0]):
                if (sns_data_scatter['Number of cells in whole dataset'][line] > annotate_thres_ncells) \
                        and (sns_data_scatter['Classwise f1 score'][line] > annotate_thres_f1):
                    axs.text(sns_data_scatter['Number of cells in whole dataset'][line] + 100,
                             sns_data_scatter['Classwise f1 score'][line],
                             sns_data_scatter.index[line],
                             horizontalalignment='left',
                             size='medium',
                             color='black',
                             weight='semibold'
                             )

        return fig, axs, sns_data_scatter


class SummarizeGridsearchEmbedding(GridsearchContainer):
    loss_idx: int
    mse_idx: int

    def __init__(
            self,
            source_path: str,
            cv: bool,
            loss_idx: int = 0,
            mse_idx: int = 1,
            model_id_len: int = 7
    ):
        super(SummarizeGridsearchEmbedding, self).__init__(
            source_path=source_path,
            cv=cv
        )
        self.loss_idx = loss_idx
        self.mse_idx = mse_idx
        self.model_id_len = model_id_len

    def create_summary_tab(self):
        metrics = list(self.evals.values())[0]['val'].keys()
        self.summary_tab = pandas.DataFrame(dict(
            list({
                "depth": [id_i.split("_")[self.model_id_len + 0] for id_i in self.run_ids],
                "width": [id_i.split("_")[self.model_id_len + 1] for id_i in self.run_ids],
                "lr": [id_i.split("_")[self.model_id_len + 2] for id_i in self.run_ids],
                "dropout": [id_i.split("_")[self.model_id_len + 3] for id_i in self.run_ids],
                "l1": [id_i.split("_")[self.model_id_len + 4] for id_i in self.run_ids],
                "l2": [id_i.split("_")[self.model_id_len + 5] for id_i in self.run_ids],
                "cv": [id_i.split("_")[-1] if self.cv else "1" for id_i in self.run_ids],
                "model": ["_".join(id_i.split("_")[:self.model_id_len]) for id_i in self.run_ids],
                "organ": [id_i.split("_")[2] for id_i in self.run_ids],
                "model_type": [id_i.split("_")[3] for id_i in self.run_ids],
                "model_gs_id": ["_".join(id_i.split("_")[:(self.model_id_len + 6)]) for id_i in self.run_ids],
                "run": self.run_ids,
        }.items()) +
            list(dict([("train_" + m, [self.evals[x]["train"][m] for x in self.run_ids]) for m in metrics]).items()) +
            list(dict([("test_" + m, [self.evals[x]["test"][m] for x in self.run_ids]) for m in metrics]).items()) +
            list(dict([("val_" + m, [self.evals[x]["val"][m] for x in self.run_ids]) for m in metrics]).items()) +
            list(dict([("all_" + m, [self.evals[x]["all"][m] for x in self.run_ids]) for m in metrics]).items())
        ))
        if self.summary_tab.shape[0] == 0:
            raise ValueError("summary_tab was empty")

    def best_model_embedding(
            self,
            subset: dict = {},
            partition: str = "val",
            metric: str = "loss",
            cvs: Union[None, List[int]] = None
    ):
        model_id, _, _ = self.get_best_model_ids(
            tab=self.summary_tab,
            partition_select=partition,
            metric_select=metric,
            cv_mode="mean",
            subset=subset,
        )
        if model_id is not None:
            if cvs is not None:
                fns = [
                    self.source_path + self.gs_keys[model_id + "_cv" + str(x)] + "/results/" + model_id + "_cv" + str(x)
                    for x in cvs
                ]
            else:
                fns = [self.source_path + self.gs_keys[model_id] + "/results/" + model_id]
            embedding = [np.load(x + "_embedding.npy") for x in fns]
            covar = [pandas.read_csv(x + "_covar.csv") for x in fns]
            return model_id, embedding, covar
        else:
            return None, [None], [None]

    def plot_best(
            self,
            rename_levels=[],
            partition_select: str = "val",
            partition_show: str = "test",
            metric_select: str = "ll",
            metric_show: str = "ll",
            collapse_cv: str = "min",
            vmin=None,
            vmax=None,
            height_fig=7,
            width_fig=7
    ):
        """

        :param rename_levels:
        :param collapse_cv:
        :param metric:
        :param ylim:
        :param xrot:
        :param height_fig:
        :param width_fig:
        :return:
        """
        import matplotlib.pyplot as plt
        import seaborn as sns

        if self.summary_tab is None:
            self.create_summary_tab()

        # Choose the best over categories based on mean loss in CV.
        # Keep variation across CV.
        sns_tab = self.best_model_by_partition(
            partition_select=partition_select,
            metric_select=metric_select,
            return_run_only=False,
            grouping=["organ", "model_type"]
        )
        for rename_level in rename_levels:
            levels_new = sns_tab[rename_level[0]].values
            levels_new[levels_new == rename_level[1]] = rename_level[2]
            sns_tab[rename_level[0]] = levels_new

        # Build figure.
        organs = np.unique(sns_tab["organ"].values)
        model_types = np.unique(sns_tab["model_type"].values)
        hm = np.zeros((len(organs), len(model_types))) + np.nan
        mask = np.isnan(hm)
        for i, m in enumerate(model_types):
            for j, o in enumerate(organs):
                data_temp = sns_tab.loc[
                    np.logical_and(
                        sns_tab["model_type"].values == m,
                        sns_tab["organ"].values == o
                    ), partition_show + "_" + metric_show
                ]
                if data_temp.shape[0] > 0:
                    if self.cv:
                        if collapse_cv == "mean":
                            hm[j, i] = np.mean(data_temp.values)
                        elif collapse_cv == "median":
                            hm[j, i] = np.median(data_temp.values)
                        elif collapse_cv == "max":
                            hm[j, i] = np.max(data_temp.values)
                        elif collapse_cv == "min":
                            hm[j, i] = np.min(data_temp.values)
                        else:
                            raise ValueError("collapse_cv % s not recognized" % collapse_cv)
                    else:
                        hm[j, i] = data_temp.values[0]
                        mask[j, i] = False
        if vmin is not None:
            hm = np.maximum(hm, np.asarray(vmin))
        if vmax is not None:
            hm = np.minimum(hm, np.asarray(vmin))
        sns_data_heatmap = pandas.DataFrame(
            hm, index=organs, columns=model_types
        )
        fig, axs = plt.subplots(1, 1, figsize=(height_fig, width_fig))
        with sns.axes_style("dark"):
            axs = sns.heatmap(
                sns_data_heatmap, #mask=mask,
                annot=True, fmt=".2f",
                ax=axs,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': partition_show + "_" + metric_show}
            )
        return fig, axs, sns_data_heatmap

    from typing import Union, List

    def get_gradients_by_celltype(
            self,
            organ: str,
            organism: str,
            model_type: Union[str, List[str]],
            metric_select: str,
            datapath,
            test_data=True,
            partition_select: str = "val",
            ignore_cache=False,
            min_cells=10,
    ):
        """
        Compute gradients across latent units with respect to input features for each cell type.

        :param organ:
        :param organism:
        :param model_type:
        :param metric_select:
        :param datapath:
        :param test_data:
        :param partition_select:
        :param ignore_cache:
        :param min_cells:
        :return: (cell types, input features) cumulative gradients
        """
        model_id, _, _ = self.get_best_model_ids(
            tab=self.summary_tab,
            metric_select=metric_select,
            partition_select=partition_select,
            subset={
                "organ": organ,
                "model_type": model_type,
            }
        )
        # check cached file

        resultspath = os.path.join(self.source_path, self.gs_keys[model_id], 'results')

        if os.path.isfile(os.path.join(resultspath, model_id + '_grads.pickle')) and not ignore_cache:
            print('Load gradients from cached file...')
            with open(os.path.join(resultspath, model_id + '_grads.pickle'), 'rb') as f:
                gradients_raw = pickle.load(f)
        else:
            print('Compute gradients (1/3): load data')
            # load data
            tz = TargetZoos(path=datapath)
            if organism == "human":
                dataset = tz.data_human[organ]
            elif organism == "mouse":
                dataset = tz.data_mouse[organ]
            else:
                raise (ValueError(f"Supplied organism {organism} not recognised. Should be one of ('mouse', 'human')"))
            dataset.load_all(annotated_only=True)

            print('Compute gradients (2/3): load embedding')
            # load embedding
            adata = dataset.adata
            topology = model_id
            embedding = EstimatorKerasEmbedding(
                data=adata,
                model_dir="",
                model_id="",
                species=organism,
                organ=organ,
                model_type=model_type,
                model_topology=model_id.split('_')[5]
            )
            embedding.init_model()
            embedding.model.training_model.load_weights(os.path.join(resultspath, model_id + '_weights.h5'))

            # compute gradients
            print('Compute gradients (3/3): cumulate gradients')
            gradients_raw = embedding.compute_gradients_input(test_data=test_data, batch_size=256, per_celltype=True)
            with open(os.path.join(resultspath, model_id + '_grads.pickle'), 'wb') as f:
                pickle.dump(gradients_raw, f, pickle.HIGHEST_PROTOCOL)
            print('Gradients saved to cache file!')

        # filter by minimum number cells min_cells
        filtered_grads = {}
        celltypes = []
        for celltype in gradients_raw['gradients'].keys():
            if gradients_raw['counts'][celltype] > min_cells:
                filtered_grads.update({celltype: gradients_raw['gradients'][celltype]})
                celltypes.append(celltype)

        return np.concatenate([
            np.mean(a, axis=0, keepdims=True)
            for a in list(filtered_grads.values())
        ], axis=0), celltypes

    def plot_gradient_distr(
            self,
            organ: str,
            organism: str,
            model_type: Union[str, List[str]],
            metric_select: str,
            datapath,
            test_data=True,
            partition_select: str = "val",
            normalize=True,
            remove_inactive=True,
            min_cells=10,
            bw=0.02,
            xlim=None,
            by_type=True,
            height_fig=7,
            width_fig=7,
            hist=False,
            ignore_cache=False,
            save=None,
    ):
        import seaborn as sns
        import matplotlib.pyplot as plt

        if by_type and isinstance(model_type, list):
            raise ValueError("cannot plot by type and by model")
        if isinstance(model_type, str):
            model_type = [model_type]

        if self.summary_tab is None:
            self.create_summary_tab()

        # best model for each organ and model_type
        avg_grads = {}
        celltypes = {}
        for modelt in model_type:
            avg_grads[modelt], celltypes[modelt] = self.get_gradients_by_celltype(
                organ=organ,
                organism=organism,
                model_type=modelt,
                metric_select=metric_select,
                datapath=datapath,
                test_data=test_data,
                partition_select=partition_select,
                ignore_cache=ignore_cache,
                min_cells=min_cells,
            )

            if normalize:
                avg_grads[modelt] = np.abs(avg_grads[modelt])
                avg_grads[modelt] = (avg_grads[modelt] - np.min(avg_grads[modelt], axis=1, keepdims=True)) / \
                                    np.maximum(
                                        np.max(avg_grads[modelt], axis=1, keepdims=True) - np.min(avg_grads[modelt],
                                                                                                  axis=1,
                                                                                                  keepdims=True), 1e-8)

        fig, axs = plt.subplots(1, 1, figsize=(width_fig, height_fig))

        if len(avg_grads.values()) == 1:
            threshold = np.mean(list(avg_grads.values())[0]) * 0.05
            avg_grads_mask = np.mean(list(avg_grads.values())[0], axis=0) > threshold
            active_grads = list(avg_grads.values())[0][:, avg_grads_mask]
            plt.axvline(threshold, color='k', linestyle='dashed', linewidth=1,
                        label="active gene threshold"
                        )
            plt.axvline(np.mean(active_grads), color='k', linestyle='solid', linewidth=1,
                        label="average gradient\nof active genes")
            print('number of active inputs: ', active_grads.shape[1])

        for k, v in avg_grads.items():
            if by_type:
                v_mask = np.mean(v, axis=0) > threshold
                for i, x in enumerate(v):
                    if remove_inactive:
                        x = x[v_mask]
                    if not hist:
                        sns.kdeplot(x, bw_method=bw, ax=axs)
            else:
                if remove_inactive:
                    threshold = np.mean(v) * 0.05
                    v_mask = np.mean(v, axis=0) > threshold
                    v = v[:, v_mask]
                if not hist:
                    sns.kdeplot(np.asarray(v).flatten(), bw_method=bw, label=k, ax=axs)

        if xlim is not None:
            axs.set_xlim(xlim)
        plt.legend(loc="best")
        plt.xlabel(r'$\rm{mean}_{i=1,...,D} \frac{\partial z_i}{\partial x}$')
        if hist:
            plt.ylabel('# genes')
        plt.tight_layout()
        if save is not None:
            plt.savefig(save)
        plt.show()

    def plot_gradient_cor(
            self,
            organ: str,
            organism: str,
            model_type: Union[str, List[str]],
            metric_select: str,
            datapath,
            test_data=True,
            partition_select: str = "val",
            height_fig=7,
            width_fig=7,
            ignore_cache=False,
            min_cells=10,
            by_type=True,
            vmin=0.,
            vmax=1.,
            save=None,
    ):
        """
        Plot correlation heatmap of gradient vectors accumulated on input features between cell types or models.

        :param organ:
        :param organism:
        :param model_type:
        :param metric_select:
        :param datapath:
        :param test_data:
        :param partition_select:
        :param height_fig:
        :param width_fig:
        :param ignore_cache:
        :param min_cells:
        :param by_type:
        :param vmin:
        :param vmax:
        :param save:
        :return:
        """
        import seaborn as sns
        import matplotlib.pyplot as plt

        if by_type and isinstance(model_type, list):
            raise ValueError("cannot plot by type and by model")
        if isinstance(model_type, str):
            model_type = [model_type]

        if self.summary_tab is None:
            self.create_summary_tab()

        # best model for each organ and model_type
        avg_grads = {}
        celltypes = {}
        for modelt in model_type:
            avg_grads[modelt], celltypes[modelt] = self.get_gradients_by_celltype(
                organ=organ,
                organism=organism,
                model_type=modelt,
                metric_select=metric_select,
                datapath=datapath,
                test_data=test_data,
                partition_select=partition_select,
                ignore_cache=ignore_cache,
                min_cells=min_cells,
            )

        fig, axs = plt.subplots(1, 1, figsize=(width_fig, height_fig))
        if by_type:
            v = avg_grads[model_type[0]]
            celltypes_coord = celltypes[model_type[0]]
            cell_names = [str(i) for i in range(v.shape[0])]
            cormat = pandas.DataFrame(
                np.corrcoef(v),
                index=celltypes_coord,
                columns=celltypes_coord
            )
            sns.heatmap(cormat, vmin=vmin, vmax=vmax, ax=axs)
        else:
            pass

        plt.tight_layout()
        if save is not None:
            plt.savefig(save)
        plt.show()