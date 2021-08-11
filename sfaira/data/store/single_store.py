import abc
import anndata
import dask.array
import dask.dataframe
import numpy as np
import os
import pandas as pd
import pickle
import scipy.sparse
import sys
import time
from typing import Dict, List, Tuple, Union

from sfaira.consts import AdataIdsSfaira, OCS
from sfaira.data.dataloaders.base.utils import is_child, UNS_STRING_META_IN_OBS
from sfaira.data.store.base import DistributedStoreBase
from sfaira.versions.genomes.genomes import GenomeContainer

"""
Distributed stores are array-like classes that sit on groups of on disk representations of anndata instances files.
Depending on the file format of the count matrix on disk, different in memory representations are sensible.
In particular, if .X is saved as zarr array, one can use lazy dask arrays to operate across sets of count matrices,
heavily reducing the complexity of the code required here and often increasing access speed.

This instances sit on groups of data objects that have the same feature space, e.g. are from the same organism.
You can operate on multiple organisms at the same time by using an umbrella class over a set of such instances by using
DistributedStoreMultipleFeatureSpaceBase.

DistributedStoreBase is base class for any file format on disk.
DistributedStoreDao wraps an on-disk representation of anndata instance in the sfaira "dao" format.
DistributedStoreH5ad wraps an on-disk representation of anndata instances as a h5ad file.
DistributedStoreAnndata wraps in-memory anndata instance.

Note that in all cases, you can use standard anndata reading functions to load a single object into memory.
"""


def _process_batch_size(x: int, idx: np.ndarray) -> int:
    if x > len(idx):
        batch_size_new = len(idx)
        print(f"WARNING: reducing retrieval batch size according to data availability in store "
              f"from {x} to {batch_size_new}")
        x = batch_size_new
    return x


class DistributedStoreSingleFeatureSpace(DistributedStoreBase):

    """
    Data set group class tailored to data access requirements common in high-performance computing (HPC).

    This class does not inherit from DatasetGroup because it entirely relies on the cached objects.
    This class is centred around .adata_by_key and .indices.

    .adata_by_key is a dictionary (by id) of backed anndata instances that point to individual h5ads.
    This dictionary is intialised with all h5ads in the store.
    As the store is sub-setted, key-value pairs are deleted from this dictionary.

    .indices have keys that correspond to keys in .adata_by_key and contain index vectors of observations in the anndata
    instances in .adata_by_key which are still kept.
    These index vectors are a form of lazy slicing that does not require data set loading or re-writing.
    As the store is sub-setted, key-value pairs are deleted from this dictionary if no observations from a given key
    match the sub-setting.
    If a subset of observations from a key matches the subsetting operation, the index set in the corresponding value is
    reduced.
    All data retrieval operations work on .indices: Generators run over these indices when retrieving observations for
    example.
    """

    _adata_by_key: Dict[str, anndata.AnnData]
    _indices: Dict[str, np.ndarray]
    _obs_by_key: Union[None, Dict[str, dask.dataframe.DataFrame]]
    data_source: str

    def __init__(self, adata_by_key: Dict[str, anndata.AnnData], indices: Dict[str, np.ndarray],
                 obs_by_key: Union[None, Dict[str, dask.dataframe.DataFrame]] = None, data_source: str = "X"):
        self.adata_by_key = adata_by_key
        self.indices = indices
        self.obs_by_key = obs_by_key
        self.ontology_container = OCS
        self._genome_container = None
        self._adata_ids_sfaira = AdataIdsSfaira()
        self.data_source = data_source
        self._celltype_universe = None

    @property
    def idx(self) -> np.ndarray:
        """
        Global indices.
        """
        idx_global = np.arange(0, np.sum([len(v) for v in self.indices.values()]))
        return idx_global

    def _validate_idx(self, idx: Union[np.ndarray, list]) -> np.ndarray:
        """
        Validate global index vector.
        """
        assert np.max(idx) < self.n_obs, f"maximum of supplied index vector {np.max(idx)} exceeds number of modelled " \
                                         f"observations {self.n_obs}"
        assert len(idx) == len(np.unique(idx)), f"there were {len(idx) - len(np.unique(idx))} repeated indices in idx"
        if isinstance(idx, np.ndarray):
            assert len(idx.shape) == 1, idx.shape
            assert idx.dtype == np.int
        else:
            assert isinstance(idx, list)
            assert isinstance(idx[0], int) or isinstance(idx[0], np.int)
            idx = np.asarray(idx)
        return idx

    @property
    def organisms_by_key(self) -> Dict[str, str]:
        """
        Data set-wise organism label as dictionary of data set keys.
        """
        ks = self.indices.keys()
        organisms = [self._adata_by_key[k].uns[self._adata_ids_sfaira.organism] for k in ks]
        # Flatten list, assumes that each data set maps to one organism:
        organisms = [x[0] if (isinstance(x, list) or isinstance(x, tuple)) else x for x in organisms]
        return dict(list(zip(ks, organisms)))

    @property
    def organism(self):
        """
        Organism of store.
        """
        organisms = np.sort(np.unique(list(self.organisms_by_key.values())))
        assert len(organisms) == 1, organisms
        return organisms[0]

    def _validate_feature_space_homogeneity(self) -> List[str]:
        """
        Assert that the data sets which were kept have the same feature names.

        :return: List of feature names in shared feature space or dictionary of list of features.
        """
        reference_k = list(self._adata_by_key.keys())[0]
        var_names = self._adata_by_key[reference_k].var_names.tolist()
        for k in list(self._adata_by_key.keys()):
            assert len(var_names) == len(self._adata_by_key[k].var_names), \
                f"number of features in store differed in object {k} compared to {reference_k}"
            assert np.all(var_names == self._adata_by_key[k].var_names), \
                f"var_names in store were not matched in object {k} compared to {reference_k}"
        return var_names

    @property
    def adata_by_key(self) -> Dict[str, anndata.AnnData]:
        """
        Anndata instance for each selected data set in store, sub-setted by selected cells.
        """
        return self._adata_by_key

    @adata_by_key.setter
    def adata_by_key(self, x: Dict[str, anndata.AnnData]):
        self._adata_by_key = x

    @property
    def data_by_key(self):
        """
        Data matrix for each selected data set in store, sub-setted by selected cells.
        """
        return dict([(k, v.X) for k, v in self.adata_by_key.items()])

    @property
    def adata_memory_footprint(self) -> Dict[str, float]:
        """
        Memory foot-print of data set k in MB.
        """
        return dict([(k, sys.getsizeof(v) / np.power(1024, 2)) for k, v in self.adata_by_key.items()])

    @property
    def indices(self) -> Dict[str, np.ndarray]:
        """
        Indices of observations that are currently exposed in adata of this instance.

        This depends on previous subsetting.
        """
        return self._indices

    @indices.setter
    def indices(self, x: Dict[str, np.ndarray]):
        """
        Setter imposes a few constraints on indices:

            1) checks that keys are contained ._adata_by_key.keys()
            2) checks that indices are contained in size of values of ._adata_by_key
            3) checks that indces are not duplicated
            4) checks that indices are sorted
        """
        for k, v in x.items():
            assert k in self._adata_by_key.keys(), f"did not find key {k}"
            assert np.max(v) < self._adata_by_key[k].n_obs, f"found index for key {k} that exceeded data set size"
            assert len(v) == len(np.unique(v)), f"found duplicated indices for key {k}"
            assert np.all(np.diff(v) >= 0), f"indices not sorted for key {k}"
        self._indices = x

    @property
    def obs_by_key(self) -> Dict[str, Union[pd.DataFrame, dask.dataframe.DataFrame]]:
        if self._obs_by_key is not None:
            # Run sanity checks to validate that this external obs can be used in the context of .adata_by_key:
            assert np.all(list(self._adata_by_key.keys()) == list(self._obs_by_key.keys()))
            assert np.all([self._obs_by_key[k].shape[0] == self._adata_by_key[k].shape[0]
                           for k in self._obs_by_key.keys()])
            return self._obs_by_key
        else:
            return dict([(k, v.obs) for k, v in self.adata_by_key.items()])

    @obs_by_key.setter
    def obs_by_key(self, x: Union[None, Dict[str, dask.dataframe.DataFrame]]):
        if x is not None:
            for k, v in x.items():
                if not (isinstance(v, dask.dataframe.DataFrame) or isinstance(v, pd.DataFrame)):
                    raise ValueError(f"value of entry {k} was not a dask.dataframe.DataFrame but {type(v)}")
        self._obs_by_key = x

    @property
    def genome_container(self) -> Union[GenomeContainer, None]:
        return self._genome_container

    @genome_container.setter
    def genome_container(self, x: Union[GenomeContainer]):
        var_names = self._validate_feature_space_homogeneity()
        # Validate genome container choice:
        # Make sure that all var names defined in genome container are also contained in loaded data sets.
        assert np.all([y in var_names for y in x.ensembl]), \
            "did not find variable names from genome container in store"
        self._genome_container = x

    @property
    def dataset_weights(self):
        return self._dataset_weights

    @dataset_weights.setter
    def dataset_weights(self, x: Dict[str, float]):
        assert np.all([k in self.adata_by_key.keys() for k in x.keys()]), "did not recognize some keys"
        assert np.all([k in x.keys() for k in self.indices.keys()]), "some data sets in index were omitted"
        self._dataset_weights = x

    def get_subset_idx(self, attr_key, values: Union[str, List[str], None],
                       excluded_values: Union[str, List[str], None]) -> dict:
        """
        Get indices of subset list of adata objects based on cell-wise properties.

        :param attr_key: Property to subset by. Options:

            - "assay_differentiation" points to self.assay_differentiation_obs_key
            - "assay_sc" points to self.assay_sc_obs_key
            - "assay_type_differentiation" points to self.assay_type_differentiation_obs_key
            - "cell_line" points to self.cell_line
            - "cellontology_class" points to self.cellontology_class_obs_key
            - "developmental_stage" points to self.developmental_stage_obs_key
            - "ethnicity" points to self.ethnicity_obs_key
            - "organ" points to self.organ_obs_key
            - "organism" points to self.organism_obs_key
            - "sample_source" points to self.sample_source_obs_key
            - "sex" points to self.sex_obs_key
            - "state_exact" points to self.state_exact_obs_key
        :param values: Classes to overlap to. Supply either values or excluded_values.
        :param excluded_values: Classes to exclude from match list. Supply either values or excluded_values.
        :return dictionary of files and observation indices by file.
        """
        if not isinstance(values, list):
            values = [values]
        assert (values is None or excluded_values is not None) or (values is not None or excluded_values is None), \
            "supply either values or excluded_values"

        def get_idx(adata, obs, k, v, xv, dataset):
            # Use cell-wise annotation if data set-wide maps are ambiguous:
            # This can happen if the different cell-wise annotations are summarised as a union in .uns.
            if getattr(self._adata_ids_sfaira, k) in adata.uns.keys() and \
                    adata.uns[getattr(self._adata_ids_sfaira, k)] != UNS_STRING_META_IN_OBS and \
                    getattr(self._adata_ids_sfaira, k) not in obs.columns:
                values_found = adata.uns[getattr(self._adata_ids_sfaira, k)]
                if isinstance(values_found, np.ndarray):
                    values_found = values_found.tolist()
                elif not isinstance(values_found, list):
                    values_found = [values_found]
                if len(values_found) > 1:
                    values_found = None  # Go to cell-wise annotation.
                else:
                    # Replicate unique property along cell dimension.
                    values_found = [values_found[0] for _ in range(adata.n_obs)]
            else:
                values_found = None
            if values_found is None:
                if getattr(self._adata_ids_sfaira, k) in obs.columns:
                    values_found = obs[getattr(self._adata_ids_sfaira, k)].values
                else:
                    values_found = []
                    print(f"WARNING: did not find attribute {k} in data set {dataset}")
            values_found_unique = np.unique(values_found)
            try:
                ontology = getattr(self.ontology_container, k)
            except AttributeError:
                raise ValueError(f"{k} not a valid property of ontology_container object")
            # Test only unique elements found in ontology to save time.
            if v is not None:
                values_found_unique_matched = [
                    x for x in values_found_unique if np.any([
                        is_child(query=x, ontology=ontology, ontology_parent=y)
                        for y in v
                    ])
                ]
            else:
                values_found_unique_matched = [
                    x for x in values_found_unique if np.all([
                        not is_child(query=x, ontology=ontology, ontology_parent=y)
                        for y in xv
                    ])
                ]
            idx = np.where([x in values_found_unique_matched for x in values_found])[0]
            return idx

        indices = {}
        for key in self.indices.keys():
            if key not in self.adata_by_key.keys():
                raise ValueError(f"data set {key} queried by indices does not exist in store (.adata_by_key)")
            # Get indices of idx_old to keep:
            adata_k = self.adata_by_key[key]
            obs_k = self.obs_by_key[key]
            idx_old = self.indices[key]
            # Cannot index on view here as indexing on view of views of backed anndata objects is not yet supported.
            idx_subset = get_idx(adata=adata_k, obs=obs_k, k=attr_key, v=values, xv=excluded_values, dataset=key)
            # Keep intersection of old and new hits.
            idx_new = np.sort(list(set(np.asarray(idx_old).tolist()).intersection(
                set(np.asarray(idx_subset).tolist()))))
            if len(idx_new) > 0:
                indices[key] = np.asarray(idx_new, dtype="int32")
        return indices

    def subset(self, attr_key, values: Union[str, List[str], None] = None,
               excluded_values: Union[str, List[str], None] = None, verbose: int = 1):
        """
        Subset list of adata objects based on cell-wise properties.

        Subsetting is done based on index vectors, the objects remain untouched.

        :param attr_key: Property to subset by. Options:

            - "assay_differentiation" points to self.assay_differentiation_obs_key
            - "assay_sc" points to self.assay_sc_obs_key
            - "assay_type_differentiation" points to self.assay_type_differentiation_obs_key
            - "cell_line" points to self.cell_line
            - "cellontology_class" points to self.cellontology_class_obs_key
            - "developmental_stage" points to self.developmental_stage_obs_key
            - "ethnicity" points to self.ethnicity_obs_key
            - "organ" points to self.organ_obs_key
            - "organism" points to self.organism_obs_key
            - "sample_source" points to self.sample_source_obs_key
            - "sex" points to self.sex_obs_key
            - "state_exact" points to self.state_exact_obs_key
        :param values: Classes to overlap to. Supply either values or excluded_values.
        :param excluded_values: Classes to exclude from match list. Supply either values or excluded_values.
        """
        self.indices = self.get_subset_idx(attr_key=attr_key, values=values, excluded_values=excluded_values)
        if self.n_obs == 0 and verbose > 0:
            print("WARNING: store is now empty.")

    def write_config(self, fn: Union[str, os.PathLike]):
        """
        Writes a config file that describes the current data sub-setting.

        This config file can be loaded later to recreate a sub-setting.
        This config file contains observation-wise subsetting information.

        :param fn: Output file without file type extension.
        """
        with open(fn + '.pickle', 'wb') as f:
            pickle.dump(self.indices, f)

    def load_config(self, fn: Union[str, os.PathLike]):
        """
        Load a config file and recreates a data sub-setting.
        This config file contains observation-wise subsetting information.

        :param fn: Output file without file type extension.
        """
        with open(fn, 'rb') as f:
            self.indices = pickle.load(f)
        # Make sure all declared data sets are found in store:
        for x in self.indices.keys():
            if x not in self.adata_by_key.keys():
                raise ValueError(f"did not find object with name {x} in currently loaded universe")

    @property
    def var_names(self) -> List[str]:
        """
        Feature names of selected genes by organism in store.
        """
        var_names = self._validate_feature_space_homogeneity()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is None:
            return var_names
        else:
            return self.genome_container.ensembl

    @property
    def n_vars(self) -> int:
        """
        Number of selected features per organism in store
        """
        var_names = self._validate_feature_space_homogeneity()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is None:
            return len(var_names)
        else:
            return self.genome_container.n_var

    @property
    def n_obs(self) -> int:
        """
        Number of observations selected in store.
        """
        return np.sum([len(v) for v in self.indices.values()])

    @property
    def shape(self) -> Tuple[int, int]:
        return self.n_obs, self.n_vars

    @abc.abstractmethod
    def _generator(
            self,
            idx_gen: iter,
            var_idx: Union[np.ndarray, None],
            obs_keys: List[str],
    ) -> iter:
        pass

    def _index_curation_helper(
            self,
            idx: Union[np.ndarray, None],
            batch_size: int,
    ) -> Tuple[Union[np.ndarray, None], Union[np.ndarray, None], int]:
        """
        Process indices and batch size input for generator production.

        Feature indicces are formatted based on previously loaded genome container.

        :param idx: Global idx to query from store. These is an array with indices corresponding to a contiuous index
            along all observations in self.adata_by_key, ordered along a hypothetical concatenation along the keys of
            self.adata_by_key. If None, all observations are seleccted.
        :param batch_size: Number of observations read from disk in each batched access (generator invocation).
        :return: Tuple:
            - idx: Processed observation index vector for generator to access.
            - var_idx Processed feature index vector for generator to access.
            - batch_size Processed retrieval batch size for generator to access.
        """
        # Make sure that features are ordered in the same way in each object so that generator yields consistent cell
        # vectors.
        var_names = self._validate_feature_space_homogeneity()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is not None:
            var_names_target = self.genome_container.ensembl
            var_idx = np.sort([var_names.index(x) for x in var_names_target])
            # Check if index vector is just full ordered list of indices, in this case, sub-setting is unnecessary.
            if len(var_idx) == len(var_names) and np.all(var_idx == np.arange(0, len(var_names))):
                var_idx = None
        else:
            var_idx = None
        # Select all cells if idx was None:
        if idx is None:
            idx = np.arange(0, self.n_obs)
        idx = self._validate_idx(idx)
        batch_size = _process_batch_size(x=batch_size, idx=idx)
        return idx, var_idx, batch_size

    def generator(
            self,
            idx: Union[np.ndarray, None] = None,
            batch_size: int = 1,
            obs_keys: List[str] = [],
            return_dense: bool = True,
            randomized_batch_access: bool = False,
            random_access: bool = False,
            **kwargs
    ) -> iter:
        """
        Yields an unbiased generator over observations in the contained data sets.

        :param idx: Global idx to query from store. These is an array with indices corresponding to a contiuous index
            along all observations in self.adata_by_key, ordered along a hypothetical concatenation along the keys of
            self.adata_by_key. If None, all observations are seleccted.
        :param batch_size: Number of observations read from disk in each batched access (generator invocation).
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adata_by_key.
        :param return_dense: Whether to force return count data .X as dense batches. This allows more efficient feature
            indexing if the store is sparse (column indexing on csr matrices is slow).
        :param randomized_batch_access: Whether to randomize batches during reading (in generator). Lifts necessity of
            using a shuffle buffer on generator, however, batch composition stays unchanged over epochs unless there
            is overhangs in retrieval_batch_size in the raw data files, which often happens and results in modest
            changes in batch composition.
            Do not use randomized_batch_access and random_access.
        :param random_access: Whether to fully shuffle observations before batched access takes place. May
            slow down access compared randomized_batch_access and to no randomization.
            Do not use randomized_batch_access and random_access.
        :return: Generator function which yields batch_size at every invocation.
            The generator returns a tuple of (.X, .obs).
        """
        idx, var_idx, batch_size = self._index_curation_helper(idx=idx, batch_size=batch_size)
        if randomized_batch_access and random_access:
            raise ValueError("Do not use randomized_batch_access and random_access.")
        n_obs = len(idx)
        remainder = n_obs % batch_size
        n_batches = int(n_obs // batch_size + int(remainder > 0))

        def idx_gen():
            """
            Yields index objects for one epoch of all data.

            These index objects are used by generators that have access to the data objects to build data batches.

            :returns: Tuple of:
                - Ordering of observations in epoch.
                - Batch start and end indices for batch based on ordering defined in first output.
            """
            batch_starts_ends = [
                (int(x * batch_size), int(np.minimum((x * batch_size) + batch_size, n_obs)))
                for x in np.arange(0, n_batches)
            ]
            batch_range = np.arange(0, len(batch_starts_ends))
            if randomized_batch_access:
                np.random.shuffle(batch_range)
            batch_starts_ends = [batch_starts_ends[i] for i in batch_range]
            obs_idx = idx.copy()
            if random_access:
                np.random.shuffle(obs_idx)
            yield obs_idx, batch_starts_ends

        return self._generator(idx_gen=idx_gen(), var_idx=var_idx, obs_keys=obs_keys), n_batches

    def generator_balanced(
            self,
            idx: Union[np.ndarray, None] = None,
            balance_obs: Union[str, None] = None,
            balance_damping: float = 0.,
            batch_size: int = 1,
            obs_keys: List[str] = [],
            **kwargs
    ) -> iter:
        """
        Yields a data set balanced generator.

        Yields one random batch per dataset. Assumes that data sets are annotated in .obs.
        Uses self.dataset_weights if this are given to sample data sets with different frequencies.
        Can additionally also balance across one meta data annotation within each data set.

        Assume you have a data set with two classes (A=80, B=20 cells) in a column named "cellontology_class".
        The single batch for this data set produced by this generator in each epoch contains N cells.
        If balance_obs is False, these N cells are the result of a draw without replacement from all 100 cells in this
        dataset in which each cell receives the same weight / success probability of 1.0.
        If balance_obs is True, these N cells are the result of a draw without replacement from all 100 cells in this
        data set with individual success probabilities such that classes are balanced: 0.2 for A and 0.8 for B.

        :param idx: Global idx to query from store. These is an array with indicies corresponding to a contiuous index
            along all observations in self.adata_by_key, ordered along a hypothetical concatenation along the keys of
            self.adata_by_key.
        :param balance_obs: .obs column key to balance samples from each data set over.
            Note that each data set must contain this column in its .obs table.
        :param balance_damping: Damping to apply to class weighting induced by balance_obs. The class-wise
            wise sampling probabilities become `max(balance_damping, (1. - frequency))`
        :param batch_size: Number of observations read from disk in each batched access (generator invocation).
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adata_by_key.
        :return: Generator function which yields batch_size at every invocation.
            The generator returns a tuple of (.X, .obs).
        """
        idx, var_idx, batch_size = self._index_curation_helper(idx=idx, batch_size=batch_size)

        def idx_gen():
            batch_starts_ends = []
            idx_proc = []
            val_dataset = self.obs[self._adata_ids_sfaira.dataset].values
            datasets = np.unique(val_dataset)
            if self.dataset_weights is not None:
                weights = np.array([self.dataset_weights[x] for x in datasets])
                p = weights / np.sum(weights)
                datasets = np.random.choice(a=datasets, replace=True, size=len(datasets), p=p)
            if balance_obs is not None:
                val_meta = self.obs[balance_obs].values
            for x in datasets:
                idx_x = np.where(val_dataset == x)[0]
                n_obs = len(idx_x)
                batch_size_o = int(np.minimum(batch_size, n_obs))
                batch_starts_ends.append(np.array([(0, batch_size_o), ]))
                if balance_obs is None:
                    p = np.ones_like(idx_x) / len(idx_x)
                else:
                    if balance_obs not in self.obs.columns:
                        raise ValueError(f"did not find column {balance_obs} in {self.organism}")
                    val_meta_x = val_meta[idx_x]
                    class_freq = dict([(y, np.mean(val_meta_x == y)) for y in np.unique(val_meta_x)])
                    class_freq_x_by_obs = np.array([class_freq[y] for y in val_meta_x])
                    damped_freq_coefficient = np.maximum(balance_damping, (1. - class_freq_x_by_obs))
                    p = np.ones_like(idx_x) / len(idx_x) * damped_freq_coefficient
                idx_x_sample = np.random.choice(a=idx_x, replace=False, size=batch_size_o, p=p)
                idx_proc.append(idx_x_sample)
            idx_proc = np.asarray(idx_proc)
            yield idx_proc, batch_starts_ends

        return self._generator(idx_gen=idx_gen(), var_idx=var_idx, obs_keys=obs_keys)

    @property
    @abc.abstractmethod
    def X(self):
        pass

    @property
    @abc.abstractmethod
    def obs(self) -> Union[pd.DataFrame]:
        pass


class DistributedStoreH5ad(DistributedStoreSingleFeatureSpace):

    in_memory: bool

    def __init__(self, in_memory: bool, **kwargs):
        super(DistributedStoreH5ad, self).__init__(**kwargs)
        self._x_as_dask = False
        self.in_memory = in_memory

    @property
    def adata_sliced(self) -> Dict[str, anndata.AnnData]:
        """
        Only exposes the subset and slices of the adata instances contained in ._adata_by_key defined in .indices.
        """
        return dict([(k, self._adata_by_key[k][v, :]) for k, v in self.indices.items()])

    @property
    def indices_global(self) -> dict:
        """
        Increasing indices across data sets which can be concatenated into a single index vector with unique entries
        for cells.

        E.g.: For two data sets of 10 cells each, the return value would be {A:[0..9], B:[10..19]}.
        Note that this operates over pre-selected indices, if this store was subsetted before resulting in only the
        second half B to be kept, the return value would be {A:[0..9], B:[10..14]}, where .indices would be
        {A:[0..9], B:[15..19]}.
        """
        counter = 0
        indices = {}
        for k, v in self.indices.items():
            indices[k] = np.arange(counter, counter + len(v))
            counter += len(v)
        return indices

    @property
    def X(self):
        if self.in_memory:
            assert np.all([isinstance(v.X, scipy.sparse.spmatrix) for v in self.adata_by_key.values()])
            return scipy.sparse.vstack([v.X for v in self.adata_by_key.values()])
        else:
            raise NotImplementedError()

    @property
    def obs(self) -> Union[pd.DataFrame]:
        """
        Assemble .obs table of subset of selected data.

        :return: .obs data frame.
        """
        return pd.concat([
            self._adata_by_key[k].obs.iloc[v, :]
            for k, v in self.indices.items()
        ], axis=0, join="inner", ignore_index=False, copy=False)

    def _generator(
            self,
            idx_gen: iter,
            var_idx: np.ndarray,
            obs_keys: List[str] = [],
            return_dense: bool = False,
    ) -> iter:
        """
        Yields data batches as defined by index sets emitted from index generator.

        :param idx_gen: Generator that yield two elements in each draw:
            - np.ndarray: The cells to emit.
            - List[Tuple[int, int]: Batch start and end indices.
        :param var_idx: The features to emit.
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adata_by_key.
        :return: Generator function which yields batch_size at every invocation.
            The generator returns a tuple of (.X, .obs).
        """
        adata_sliced = self.adata_sliced
        # Speed up access to single object by skipping index overlap operations:
        single_object = len(adata_sliced.keys()) == 1
        if not single_object:
            idx_dict_global = dict([(k, set(v)) for k, v in self.indices_global.items()])

        def generator():
            for idx, batch_starts_ends in idx_gen:
                for s, e in batch_starts_ends:
                    idx_i = idx[s:e]
                    # Match adata objects that overlap to batch:
                    if single_object:
                        idx_i_dict = dict([(k, np.sort(idx_i)) for k in adata_sliced.keys()])
                    else:
                        idx_i_set = set(idx_i)
                        idx_i_dict = dict([
                            (k, np.sort(list(idx_i_set.intersection(v))))
                            for k, v in idx_dict_global.items()
                        ])
                    # Only retain non-empty.
                    idx_i_dict = dict([(k, v) for k, v in idx_i_dict.items() if len(v) > 0])
                    # I) Prepare data matrix.
                    x = [
                        adata_sliced[k].X[v, :]
                        for k, v in idx_i_dict.items()
                    ]
                    # Move from ArrayView to numpy if backed and dense:
                    x = [
                        xx.toarray() if isinstance(xx, anndata._core.views.ArrayView) else xx
                        for xx in x
                    ]
                    # Do dense conversion now so that col-wise indexing is not slow, often, dense conversion
                    # would be done later anyway.
                    if return_dense:
                        x = [np.asarray(xx.todense()) if isinstance(xx, scipy.sparse.spmatrix) else xx for xx in x]
                        is_dense = True
                    else:
                        is_dense = isinstance(x[0], np.ndarray)
                    # Concatenate blocks in observation dimension:
                    if len(x) > 1:
                        if is_dense:
                            x = np.concatenate(x, axis=0)
                        else:
                            x = scipy.sparse.vstack(x)
                    else:
                        x = x[0]
                    if var_idx is not None:
                        x = x[:, var_idx]
                    # Prepare .obs.
                    obs = pd.concat([
                        adata_sliced[k].obs[obs_keys].iloc[v, :]
                        for k, v in idx_i_dict.items()
                    ], axis=0, join="inner", ignore_index=True, copy=False)
                    yield x, obs

        return generator


class DistributedStoreDao(DistributedStoreSingleFeatureSpace):

    _dataset_weights: Union[None, Dict[str, float]]
    _x: Union[None, dask.array.Array]
    _x_by_key: Union[None, dask.array.Array]

    def __init__(self, x_by_key, **kwargs):
        super(DistributedStoreDao, self).__init__(**kwargs)
        self._x = None
        self._x_as_dask = True
        self._x_by_key = x_by_key

    @property
    def indices(self) -> Dict[str, np.ndarray]:
        return super(DistributedStoreDao, self).indices

    @indices.setter
    def indices(self, x: Dict[str, np.ndarray]):
        """
        Extends setter in super class by wiping .X cache.

        Setter imposes a few constraints on indices:

            1) checks that keys are contained ._adata_by_key.keys()
            2) checks that indices are contained in size of values of ._adata_by_key
            3) checks that indces are not duplicated
            4) checks that indices are sorted
        """
        self._x = None
        for k, v in x.items():
            assert k in self._adata_by_key.keys(), f"did not find key {k}"
            assert np.max(v) < self._adata_by_key[k].n_obs, f"found index for key {k} that exceeded data set size"
            assert len(v) == len(np.unique(v)), f"found duplicated indices for key {k}"
            assert np.all(np.diff(v) >= 0), f"indices not sorted for key {k}"
        self._indices = x

    @property
    def data_by_key(self):
        """
        Data matrix for each selected data set in store, sub-setted by selected cells.
        """
        # Accesses _x_by_key rather than _adata_by_key as long as the dask arrays are stored there.
        return dict([(k, self._x_by_key[k][v, :]) for k, v in self.indices.items()])

    @property
    def X(self) -> dask.array.Array:
        """
        One dask array of all cells.

        Requires feature dimension to be shared.
        """
        if self._x is None:
            if self.data_source == "X":
                # TODO avoiding anndata .X here
                # assert np.all([isinstance(self._adata_by_key[k].X, dask.array.Array) for k in self.indices.keys()])
                assert np.all([isinstance(self._x_by_key[k], dask.array.Array) for k in self.indices.keys()])
                self._x = dask.optimize(dask.array.vstack([
                    self._x_by_key[k][v, :]
                    for k, v in self.indices.items()
                ]))[0]
            else:
                raise ValueError(f"Did not recognise data_source={self.data_source}.")
        return self._x

    @property
    def obs(self) -> pd.DataFrame:
        """
        Assemble .obs table of subset of selected data.

        Resulting index is increasing integers starting with zero.

        :return: .obs data frame.
        """
        # TODO Using loc indexing here instead of iloc, this might be faster on larger tables?
        return pd.concat([
            self.adata_by_key[k].obs.loc[self.adata_by_key[k].obs.index[v], :]
            for k, v in self.indices.items()
        ], axis=0, join="inner", ignore_index=True, copy=False)

    def _generator(
            self,
            idx_gen: iter,
            var_idx: np.ndarray,
            obs_keys: List[str] = [],
    ) -> iter:
        """
        Yields data batches as defined by index sets emitted from index generator.

        :param idx_gen: Generator that yield two elements in each draw:
            - np.ndarray: The cells to emit.
            - List[Tuple[int, int]: Batch start and end indices.
        :param var_idx: The features to emit.
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adata_by_key.
        :return: Generator function which yields batch_size at every invocation.
            The generator returns a tuple of (.X, .obs).
        """
        # Normalise cell indices such that each organism is indexed starting at zero:
        # This is required below because each organism is represented as its own dask array.
        # TODO this might take a lot of time as the dask array is built.
        t0 = time.time()
        x = self.X
        print(f"init X: {time.time() - t0}")
        t0 = time.time()
        obs = self.obs[obs_keys]
        # Redefine index so that .loc indexing can be used instead of .iloc indexing:
        obs.index = np.arange(0, obs.shape[0])
        print(f"init obs: {time.time() - t0}")

        def generator():
            # Can all data sets corresponding to one organism as a single array because they share the second dimension
            # and dask keeps expression data and obs out of memory.
            for idx, batch_starts_ends in idx_gen:
                idx = np.sort(idx)
                x_temp = x[idx, :]
                obs_temp = obs.loc[obs.index[idx], :]  # TODO better than iloc?
                for s, e in batch_starts_ends:
                    x_i = x_temp[s:e, :]
                    if var_idx is not None:
                        x_i = x_i[:, var_idx]
                    # Exploit fact that index of obs is just increasing list of integers, so we can use the .loc[]
                    # indexing instead of .iloc[]:
                    obs_i = obs_temp.loc[obs_temp.index[s:e], :]
                    yield x_i, obs_i

        return generator
