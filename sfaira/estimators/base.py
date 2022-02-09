import os
import warnings
from typing import List, Union

import anndata
import numpy as np
import pandas as pd

from sfaira.consts import AdataIdsSfaira, OCS, AdataIds
from sfaira.data.store.stores.base import StoreBase
from sfaira.data.store.stores.multi import StoresAnndata
from sfaira.data.store.stores.single import StoreSingleFeatureSpace
from sfaira.models.base import BasicModel
from sfaira.versions.metadata import CelltypeUniverse, OntologyCl, OntologyObo
from sfaira.versions.topologies import TopologyContainer


def split_idx(data: StoreSingleFeatureSpace, test_split, val_split):
    """
    Split training and evaluation data.
    """
    np.random.seed(1)
    all_idx = np.arange(0, data.n_obs)  # n_obs is both a property of AnnData and DistributedStoreBase
    if isinstance(test_split, float) or isinstance(test_split, int):
        idx_test = np.sort(np.random.choice(
            a=all_idx,
            size=round(data.n_obs * test_split),
            replace=False,
        ))
    elif isinstance(test_split, dict):
        in_test = np.ones((data.n_obs,), dtype=int) == 1
        for k, v in test_split.items():
            if isinstance(v, bool) or isinstance(v, int) or isinstance(v, str):
                v = [v]
            elif isinstance(v, tuple):
                v = list(v)
            elif isinstance(v, np.ndarray):
                v = v.tolist()
            if not isinstance(v, list):
                raise ValueError(f"conversion of v {v} to list failed")
            if np.any([isinstance(vi, list) or isinstance(vi, tuple) or isinstance(vi, np.ndarray) for vi in v]):
                raise ValueError(f"found nested list v {v}, only use bool, scalar or string elements in v.")
            idx = data.get_subset_idx(attr_key=k, values=v, excluded_values=None)
            # Build continuous vector across all sliced data sets and establish which observations are kept
            # in subset.
            in_test_k = np.ones((data.n_obs,), dtype=int) == 0
            counter = 0
            for kk, vv in data.indices.items():
                if kk in idx.keys() and len(idx[kk]) > 0:
                    in_test_k[np.where([x in idx[kk] for x in vv])[0] + counter] = True
                counter += len(vv)
            in_test = np.logical_and(in_test, in_test_k)
        idx_test = np.sort(np.where(in_test)[0])
    else:
        raise ValueError("type of test_split %s not recognized" % type(test_split))
    print(f"Found {len(idx_test)} out of {data.n_obs} cells that correspond to test data set")
    assert len(idx_test) < data.n_obs, f"test set covers full data set, apply a more restrictive test " \
                                       f"data definiton ({len(idx_test)}, {data.n_obs})"
    idx_train_eval = all_idx[~np.isin(all_idx, idx_test)]
    np.random.seed(1)
    idx_eval = np.sort(np.random.choice(
        a=idx_train_eval,
        size=round(len(idx_train_eval) * val_split),
        replace=False
    ))
    idx_train = np.sort(idx_train_eval[~np.isin(idx_train_eval, idx_eval)])

    # Check that none of the train, test, eval partitions are empty
    if not len(idx_test):
        warnings.warn("Test partition is empty!")
    if not len(idx_eval):
        raise ValueError("The evaluation dataset is empty.")
    if not len(idx_train):
        raise ValueError("The train dataset is empty.")
    return idx_train, idx_eval, idx_test


class EstimatorBase:

    data: StoreSingleFeatureSpace
    model: Union[BasicModel, None]
    topology_container: TopologyContainer
    model_id: Union[str, None]
    weights: Union[np.ndarray, None]
    model_dir: Union[str, None]
    history: Union[dict, None]
    train_hyperparam: Union[dict, None]
    idx_train: Union[np.ndarray, None]
    idx_eval: Union[np.ndarray, None]
    idx_test: Union[np.ndarray, None]
    adata_ids: AdataIds

    def __init__(
            self,
            data: Union[anndata.AnnData, List[anndata.AnnData], StoreSingleFeatureSpace],
            model_dir: Union[str, None],
            model_class: str,
            model_id: Union[str, None],
            model_topology: TopologyContainer,
            weights_md5: Union[str, None] = None,
            cache_path: str = os.path.join('cache', ''),
            adata_ids: AdataIds = AdataIdsSfaira()
    ):
        self.model = None
        self.model_dir = model_dir
        self.model_id = model_id
        self.model_class = model_class
        self.topology_container = model_topology
        if isinstance(data, anndata.AnnData):
            data = StoresAnndata(adatas=data).stores[self.organism]
        if isinstance(data, list) or isinstance(data, tuple):
            for x in data:
                assert isinstance(x, anndata.AnnData), f"found element in list that was not anndata but {type(x)}"
            data = StoresAnndata(adatas=data).stores[self.organism]
        self.data = data
        # Prepare store with genome container sub-setting:
        # This class is tailored for DistributedStoreSingleFeatureSpace but we test for the base class here in the
        # constructor so that genome_container can also be set in inheriting classes that may be centred around
        # different child classes of DistributedStoreBase.
        if isinstance(self.data, StoreBase):
            self.data.genome_container = self.topology_container.gc

        self.history = None
        self.train_hyperparam = None
        self.idx_train = None
        self.idx_eval = None
        self.idx_test = None
        self.md5 = weights_md5
        self.cache_path = cache_path
        self._adata_ids = adata_ids

    @property
    def model_type(self):
        return self.topology_container.model_type

    @property
    def organism(self):
        return self.topology_container.organism

    def _get_class_dict(
            self,
            obs_key: str
    ):
        y = self.obs[obs_key]
        for i, val in enumerate(y):
            if type(val) == list:
                y[i] = " / ".join(val)
        labels = np.unique(y)
        label_dict = {}
        for i, label in enumerate(labels):
            label_dict.update({label: float(i)})
        return label_dict

    def split_train_val_test(self, val_split: float, test_split: Union[float, dict]):
        """
        Split indices in store into train, valiation and test split.
        """
        idx_train, idx_eval, idx_test = split_idx(data=self.data, test_split=test_split, val_split=val_split)
        self.idx_train = idx_train
        self.idx_eval = idx_eval
        self.idx_test = idx_test

    @property
    def using_store(self) -> bool:
        return isinstance(self.data, StoreSingleFeatureSpace)

    @property
    def obs_train(self) -> pd.DataFrame:
        return self.data.checkout(idx=self.idx_train).obs

    @property
    def obs_eval(self) -> pd.DataFrame:
        return self.data.checkout(idx=self.idx_eval).obs

    @property
    def obs_test(self) -> pd.DataFrame:
        return self.data.checkout(idx=self.idx_test).obs

    @property
    def obs(self) -> pd.DataFrame:
        return self.data.checkout(idx=self._process_idx_for_eval(None)).obs

    def _process_idx_for_eval(self, idx):
        """
        Defaults to all observations if no indices are defined.
        """
        if idx is None:
            idx = np.arange(0, self.data.n_obs)
        return idx


class EstimatorBaseEmbedding(EstimatorBase):

    pass


class EstimatorBaseCelltype(EstimatorBase):

    def __init__(
            self,
            data: Union[anndata.AnnData, StoreSingleFeatureSpace],
            model_dir: Union[str, None],
            model_id: Union[str, None],
            model_topology: TopologyContainer,
            weights_md5: Union[str, None] = None,
            cache_path: str = os.path.join('cache', ''),
            celltype_ontology: Union[OntologyObo, None] = None,
            max_class_weight: float = 1e3,
            remove_unlabeled_cells: bool = True,
            adata_ids: AdataIds = AdataIdsSfaira()
    ):
        super(EstimatorBaseCelltype, self).__init__(
            data=data,
            model_dir=model_dir,
            model_class="celltype",
            model_id=model_id,
            model_topology=model_topology,
            weights_md5=weights_md5,
            cache_path=cache_path,
            adata_ids=adata_ids
        )

        if remove_unlabeled_cells:
            # Remove cells without type label from store:
            self.data.subset(attr_key="cell_type", excluded_values=[
                self._adata_ids.unknown_metadata_identifier,
                self._adata_ids.not_a_cell_celltype_identifier,
            ])
        assert "cl" in self.topology_container.output.keys(), self.topology_container.output.keys()
        assert "targets" in self.topology_container.output.keys(), self.topology_container.output.keys()
        self.max_class_weight = max_class_weight
        if celltype_ontology is None:
            celltype_ontology = OntologyCl(branch=self.topology_container.output["cl"])
        self.celltype_universe = CelltypeUniverse(
            cl=celltype_ontology,
            uberon=OCS.organ,
        )
        # Set leaves if they are defined in topology:
        if self.topology_container.output["targets"] is not None:
            self.celltype_universe.onto_cl.leaves = self.topology_container.output["targets"]

    @property
    def ntypes(self):
        return self.celltype_universe.onto_cl.n_leaves

    @property
    def ontology_ids(self):
        return self.celltype_universe.onto_cl.convert_to_id(self.celltype_universe.onto_cl.leaves)

    @property
    def ontology_names(self):
        return self.celltype_universe.onto_cl.convert_to_name(self.celltype_universe.onto_cl.leaves)

    def _one_hot_encoder(self):
        leave_maps = self.celltype_universe.onto_cl.prepare_maps_to_leaves(include_self=True)

        def encoder(x) -> np.ndarray:
            if isinstance(x, str):
                x = [x]
            # Encodes unknowns to empty rows.
            idx = [
                leave_maps[y] if y not in [
                    self._adata_ids.unknown_metadata_identifier,
                    self._adata_ids.not_a_cell_celltype_identifier,
                ] else np.array([])
                for y in x
            ]
            oh = np.zeros((len(x), self.ntypes,), dtype="float32")
            for i, y in enumerate(idx):
                scale = len(y)
                if scale > 0:
                    oh[i, y] = 1. / scale
            return oh

        return encoder

    def _get_celltype_out(self, idx: Union[np.ndarray, None]):
        """
        TODO depreceate, carry over weight code to _get_generator
        Build one hot encoded cell type output tensor and observation-wise weight matrix.
        """
        if idx is None:
            idx = np.arange(0, self.data.n_obs)
        # One whether "unknown" is already included, otherwise add one extra column.
        onehot_encoder = self._one_hot_encoder()
        y = np.concatenate([
            onehot_encoder(z)
            for z in self.obs[self._adata_ids.cell_type + self._adata_ids.onto_id_suffix].values[idx].tolist()
        ], axis=0)
        # Distribute aggregated class weight for computation of weights:
        freq = np.mean(y / np.sum(y, axis=1, keepdims=True), axis=0, keepdims=True)
        weights = 1. / np.matmul(y, freq.T)  # observation wise weight matrix
        # Threshold weights:
        weights = np.asarray(
            np.minimum(weights, np.zeros_like(weights) + self.max_class_weight),
            dtype="float32"
        ).flatten()
        return weights, y
