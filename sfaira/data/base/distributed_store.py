import abc
import anndata
import dask.array
import numpy as np
import os
import pandas as pd
import pickle
import scipy.sparse
import sys
from typing import Dict, List, Tuple, Union

from sfaira.consts import AdataIdsSfaira, OCS
from sfaira.data.base.dataset import is_child, UNS_STRING_META_IN_OBS
from sfaira.data.base.zarr_andata import read_zarr
from sfaira.versions.genomes import GenomeContainer

"""
Distributed stores are array-like classes that sit on groups of on disk representations of anndata instances files.
Depending on the file format of the count matrix on disk, different in memory representations are sensible.
In particular, if .X is saved as zarr array, one can use lazy dask arrays to operate across sets of count matrices,
heavily reducing the complexity of the code required here and often increasing access speed.

DistributedStoreBase is base class for any file format on disk.
DistributedStoreZarr is adapted to classes that store an anndata instance as a zarr group.
DistributedStoreH5ad is adapted to classes that store an anndata instance as a h5ad file.

Note that in all cases, you can use standard anndata reading functions to load a single object into memory.
"""


def access_helper(adata, s, e, j, return_dense, obs_keys) -> tuple:
    x = adata.X[s:e, :]
    # Do dense conversion now so that col-wise indexing is not slow, often, dense conversion
    # would be done later anyway.
    if return_dense and isinstance(x, scipy.sparse.spmatrix):
        x = x.todense()
    if j is not None:
        x = x[:, j]
    obs = adata.obs[obs_keys].iloc[s:e, :]
    return x, obs


class DistributedStoreBase(abc.ABC):
    """
    Data set group class tailored to data access requirements common in high-performance computing (HPC).

    This class does not inherit from DatasetGroup because it entirely relies on the cached objects.
    This class is centred around .adatas and .indices.

    .adatas is a dictionary (by id) of backed anndata instances that point to individual h5ads.
    This dictionary is intialised with all h5ads in the store.
    As the store is subsetted, key-value pairs are deleted from this dictionary.

    .indices have keys that correspond to keys in .adatas and contain index vectors of observations in the anndata
    instances in .adatas which are still kept.
    These index vectors are a form of lazy slicing that does not require data set loading or re-writing.
    As the store is subsetted, key-value pairs are deleted from this dictionary if no observations from a given key
    match the subsetting.
    If a subset of observations from a key matches the subsetting operation, the index set in the corresponding value is
    reduced.
    All data retrievel operations work on .indices: Generators run over these indices when retrieving observations for
    example.
    """

    _adatas: Dict[str, anndata.AnnData]
    _indices: Dict[str, np.ndarray]

    def __init__(self, adatas: Dict[str, anndata.AnnData], indices: Dict[str, np.ndarray]):
        self.adatas = adatas
        self.indices = indices
        self.ontology_container = OCS
        self._genome_container = None
        self._adata_ids_sfaira = AdataIdsSfaira()
        self._celltype_universe = None

    def _validate_idx(self, idx: Union[np.ndarray, list]) -> np.ndarray:
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

    def _validate_feature_space_homogeneity(self) -> List[str]:
        """
        Assert that the data sets which were kept have the same feature names.
        """
        var_names = self._adatas[list(self.indices.keys())[0]].var_names.tolist()
        for k, v in self.indices.items():
            assert len(var_names) == len(self._adatas[k].var_names), \
                f"number of features in store differed in object {k} compared to {list(self._adatas.keys())[0]}"
            assert np.all(var_names == self._adatas[k].var_names), \
                f"var_names in store were not matched in object {k} compared to {list(self._adatas.keys())[0]}"
        return var_names

    def _generator_helper(
            self,
            idx: Union[np.ndarray, None] = None,
    ) -> Tuple[Union[np.ndarray, None], Union[np.ndarray, None]]:
        # Make sure that features are ordered in the same way in each object so that generator yields consistent cell
        # vectors.
        _ = self._validate_feature_space_homogeneity()
        var_names_store = self.adata_dict[list(self.indices.keys())[0]].var_names.tolist()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is not None:
            var_names_target = self.genome_container.ensembl
            var_idx = np.sort([var_names_store.index(x) for x in var_names_target])
            # Check if index vector is just full ordered list of indices, in this case, sub-setting is unnecessary.
            if len(var_idx) == len(var_names_store) and np.all(var_idx == np.arange(0, len(var_names_store))):
                var_idx = None
        else:
            var_idx = None
        if idx is not None:
            idx = self._validate_idx(idx)
        return idx, var_idx

    @property
    def adatas(self) -> Dict[str, anndata.AnnData]:
        return self._adatas

    @adatas.setter
    def adatas(self, x: Dict[str, anndata.AnnData]):
        self._adatas = x

    @property
    def indices(self) -> Dict[str, np.ndarray]:
        return self._indices

    @indices.setter
    def indices(self, x: Dict[str, np.ndarray]):
        """
        Setter imposes a few constraints on indices:

            1) checks that keys are contained ._adatas.keys()
            2) checks that indices are contained in size of values of ._adatas
            3) checks that indces are not duplicated
            4) checks that indices are sorted
        """
        for k, v in x.items():
            assert k in self._adatas.keys(), f"did not find key {k}"
            assert np.max(v) < self._adatas[k].n_obs, f"found index for key {k} that exceeded data set size"
            assert len(v) == len(np.unique(v)), f"found duplicated indices for key {k}"
            assert np.all(np.diff(v) >= 0), f"indices not sorted for key {k}"
        self._indices = x

    @property
    def genome_container(self) -> Union[GenomeContainer, None]:
        return self._genome_container

    @genome_container.setter
    def genome_container(self, x: GenomeContainer):
        var_names = self._validate_feature_space_homogeneity()
        # Validate genome container choice:
        # Make sure that all var names defined in genome container are also contained in loaded data sets.
        assert np.all([y in var_names for y in x.ensembl]), \
            "did not find variable names from genome container in store"
        self._genome_container = x

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

        def get_idx(adata, k, v, xv, dataset):
            # Use cell-wise annotation if data set-wide maps are ambiguous:
            # This can happen if the different cell-wise annotations are summarised as a union in .uns.
            if getattr(self._adata_ids_sfaira, k) in adata.uns.keys() and \
                    adata.uns[getattr(self._adata_ids_sfaira, k)] != UNS_STRING_META_IN_OBS:
                values_found = adata.uns[getattr(self._adata_ids_sfaira, k)]
                if isinstance(values_found, np.ndarray):
                    values_found = values_found.tolist()
                elif not isinstance(values_found, list):
                    values_found = [values_found]
                if len(values_found) > 1:
                    values_found = None  # Go to cell-wise annotation.
                else:
                    # Replicate unique property along cell dimension.
                    values_found = [values_found[0] for i in range(adata.n_obs)]
            else:
                values_found = None
            if values_found is None:
                if getattr(self._adata_ids_sfaira, k) in adata.obs.keys():
                    values_found = adata.obs[getattr(self._adata_ids_sfaira, k)].values
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
            # TODO keep this logging for now to catch undesired behaviour resulting from loaded edges in ontologies.
            print(f"matched cell-wise keys {str(values_found_unique_matched)} in data set {dataset}")
            idx = np.where([x in values_found_unique_matched for x in values_found])[0]
            return idx

        indices = {}
        for k, adata_k in self.adata_dict.items():
            if k not in self.adatas.keys():
                raise ValueError(f"data set {k} queried by indices does not exist in store (.adatas)")
            # Get indices of idx_old to keep:
            idx_subset = get_idx(adata=adata_k, k=attr_key, v=values, xv=excluded_values, dataset=k)
            # Keep intersection of old and new hits.
            idx_old = self.indices[k]
            idx_new = np.asarray(idx_old)[idx_subset]
            if len(idx_new) > 0:
                indices[k] = np.asarray(idx_new, dtype="int32")
        return indices

    def subset(self, attr_key, values: Union[str, List[str], None] = None,
               excluded_values: Union[str, List[str], None] = None):
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
        if self.n_obs == 0:
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
        # Subset to described data sets:
        for x in self.indices.keys():
            if x not in self.adatas.keys():
                raise ValueError(f"did not find object with name {x} in currently loaded universe")
        # Only retain data sets with which are mentioned in config file.
        self.subset(attr_key="id", values=list(self.indices.keys()))

    @property
    def var_names(self):
        var_names = self._validate_feature_space_homogeneity()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is None:
            return var_names
        else:
            return self.genome_container.ensembl

    @property
    def n_vars(self):
        var_names = self._validate_feature_space_homogeneity()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is None:
            return len(var_names)
        else:
            return self.genome_container.n_var

    @property
    def n_obs(self):
        return np.sum([len(v) for v in self.indices.values()])

    @property
    def shape(self):
        return [self.n_obs, self.n_vars]

    @property
    def obs(self) -> Union[pd.DataFrame]:
        """
        Assemble .obs table of subset of selected data.

        :return: .obs data frame.
        """
        return pd.concat([
            self._adatas[k].obs.iloc[v, :]
            for k, v in self.indices.items()
        ], axis=0)

    @abc.abstractmethod
    def generator(
            self,
            idx: Union[np.ndarray, None] = None,
            batch_size: int = 1,
            obs_keys: List[str] = [],
            return_dense: bool = True,
            randomized_batch_access: bool = False,
    ) -> iter:
        pass

    @property
    @abc.abstractmethod
    def adata_dict(self) -> Dict[str, anndata.AnnData]:
        pass

    @abc.abstractmethod
    def X(self) -> Union[dask.array.Array, scipy.sparse.csr_matrix]:
        pass

    @abc.abstractmethod
    def n_counts(self, idx: Union[np.ndarray, list, None] = None) -> np.ndarray:
        pass


class DistributedStoreH5ad(DistributedStoreBase):

    def __init__(self, cache_path: Union[str, os.PathLike]):
        # Collect all data loaders from files in directory:
        adatas = {}
        indices = {}
        for f in os.listdir(cache_path):
            adata = None
            trial_path = os.path.join(cache_path, f)
            if os.path.isfile(trial_path):
                # Narrow down to supported file types:
                if f.split(".")[-1] == "h5ad":
                    print(f"Discovered {f} as .h5ad file.")
                    try:
                        adata = anndata.read_h5ad(
                            filename=trial_path,
                            backed="r",
                        )
                    except OSError as e:
                        adata = None
                        print(f"WARNING: for data set {f}: {e}")
            if adata is not None:
                adatas[adata.uns["id"]] = adata
                indices[adata.uns["id"]] = np.arange(0, adata.n_obs)
        self._x_as_dask = False
        super(DistributedStoreH5ad, self).__init__(adatas=adatas, indices=indices)

    @property
    def adata_dict(self) -> Dict[str, anndata.AnnData]:
        """
        Only exposes the subset and slices of the adata instances contained in ._adatas defined in .indices.
        """
        return dict([(k, self._adatas[k][v, :]) for k, v in self.indices.items()])

    @property
    def X(self):
        assert False

    def n_counts(self, idx: Union[np.ndarray, list, None] = None) -> np.ndarray:
        """
        Compute sum over features for each observation in index.

        :param idx: Index vector over observations in object.
        :return: Array with sum per observations: (number of observations in index,)
        """
        return np.concatenate([
            np.asarray(v.X.sum(axis=1)).flatten()
            for v in self.adatas_subset(idx=idx).values()
        ], axis=0)

    def generator(
            self,
            idx: Union[np.ndarray, None] = None,
            batch_size: int = 1,
            obs_keys: List[str] = [],
            return_dense: bool = True,
            randomized_batch_access: bool = False,
    ) -> iter:
        """
        Yields an unbiased generator over observations in the contained data sets.

        :param idx: Global idx to query from store. These is an array with indicies corresponding to a contiuous index
            along all observations in self.adatas, ordered along a hypothetical concatenation along the keys of
            self.adatas.
        :param batch_size: Number of observations read from disk in each batched access (generator invocation).
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adatas.
        :param return_dense: Whether to force return count data .X as dense batches. This allows more efficient feature
            indexing if the store is sparse (column indexing on csr matrices is slow).
        :param randomized_batch_access: Whether to randomize batches during reading (in generator). Lifts necessity of
            using a shuffle buffer on generator, however, batch composition stays unchanged over epochs unless there
            is overhangs in retrieval_batch_size in the raw data files, which often happens and results in modest
            changes in batch composition.
        :return: Generator function which yields batch_size at every invocation.
            The generator returns a tuple of (.X, .obs) with types:

                - if store format is h5ad: (Union[scipy.sparse.csr_matrix, np.ndarray], pandas.DataFrame)
        """
        idx, var_idx = self._generator_helper(idx=idx)

        def generator():
            adatas_sliced_subset = self.adatas_subset(idx=idx)
            key_batch_starts_ends = []  # List of tuples of data set key and (start, end) index set of batches.
            for k, adata in adatas_sliced_subset.items():
                n_obs = adata.shape[0]
                if n_obs > 0:  # Skip data objects without matched cells.
                    # Cells left over after batching to batch size, accounting for overhang:
                    remainder = n_obs % batch_size
                    key_batch_starts_ends_k = [
                        (k, (int(x * batch_size), int(np.minimum((x * batch_size) + batch_size, n_obs))))
                        for x in np.arange(0, n_obs // batch_size + int(remainder > 0))
                    ]
                    assert np.sum([v2 - v1 for k, (v1, v2) in key_batch_starts_ends_k]) == n_obs
                    key_batch_starts_ends.extend(key_batch_starts_ends_k)
            batch_range = np.arange(0, len(key_batch_starts_ends))
            if randomized_batch_access:
                np.random.shuffle(batch_range)
            for i in batch_range:
                k, (s, e) = key_batch_starts_ends[i]
                x, obs = access_helper(adata=adatas_sliced_subset[k], s=s, e=e, j=var_idx, return_dense=return_dense,
                                       obs_keys=obs_keys)
                yield x, obs

        return generator

    def adatas_subset(self, idx: Union[np.ndarray, list]) -> Dict[str, anndata.AnnData]:
        """
        Subsets adatas_sliced as if it was one object, ie behaves the same way as self.adata[idx]  without explicitly
        concatenating.
        """
        if idx is not None:
            idx = self._validate_idx(idx)
            indices_subsetted = {}
            counter = 0
            for k, v in self.indices.items():
                n_obs_k = len(v)
                indices_global = np.arange(counter, counter + n_obs_k)
                indices_subset_k = [x for x, y in zip(v, indices_global) if y in idx]
                if len(indices_subset_k) > 0:
                    indices_subsetted[k] = indices_subset_k
                counter += n_obs_k
            assert counter == self.n_obs
            return dict([(k, self._adatas[k][v, :]) for k, v in indices_subsetted.items()])
        else:
            return self.adata_dict

    def get_subset_idx_global(self, attr_key, values: Union[str, List[str], None] = None,
                              excluded_values: Union[str, List[str], None] = None) -> np.ndarray:
        """
        Get indices of subset list of adata objects based on cell-wise properties treating instance as single array.

        The indices are continuous across all data sets as if they were one array.

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
        :param values: Classes to overlap to.
        :return Index vector
        """
        # Get indices of of cells in target set by file.
        idx_by_dataset = self.get_subset_idx(attr_key=attr_key, values=values, excluded_values=excluded_values)
        # Translate file-wise indices into global index list across all data sets.
        idx = []
        counter = 0
        for k, v in self.indices.items():
            for x in v:
                if k in idx_by_dataset.keys() and x in idx_by_dataset[k]:
                    idx.append(counter)
                counter += 1
        return np.asarray(idx)

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


class DistributedStoreZarr(DistributedStoreBase):

    def __init__(self, cache_path: Union[str, os.PathLike]):
        # Collect all data loaders from files in directory:
        adatas = {}
        indices = {}
        for f in os.listdir(cache_path):
            adata = None
            trial_path = os.path.join(cache_path, f)
            if os.path.isdir(trial_path):
                # zarr-backed anndata are saved as directories with the elements of the array group as further sub
                # directories, e.g. a directory called "X", and a file ".zgroup" which identifies the zarr group.
                if [".zgroup" in os.listdir(trial_path)]:
                    adata = read_zarr(trial_path, use_dask=True)
                    print(f"Discovered {f} as zarr group, "
                          f"sized {round(sys.getsizeof(adata) / np.power(1024, 2), 1)}MB")
            if adata is not None:
                adatas[adata.uns["id"]] = adata
                indices[adata.uns["id"]] = np.arange(0, adata.n_obs)
        self._x_as_dask = True
        super(DistributedStoreZarr, self).__init__(adatas=adatas, indices=indices)

    @property
    def adata_dict(self) -> Dict[str, anndata.AnnData]:
        """
        Only exposes the subset and slices of the adata instances contained in ._adatas defined in .indices.
        """
        return dict([(k, self._adatas[k][v, :]) for k, v in self.indices.items()])

    @property
    def X(self) -> Union[dask.array.Array]:
        assert np.all([isinstance(self._adatas[k].X, dask.array.Array) for k in self.indices.keys()])
        return dask.array.vstack([
            self._adatas[k].X[v, :]
            for k, v in self.indices.items()
        ])

    def n_counts(self, idx: Union[np.ndarray, list, None] = None) -> np.ndarray:
        """
        Compute sum over features for each observation in index.

        :param idx: Index vector over observations in object.
        :return: Array with sum per observations: (number of observations in index,)
        """
        return np.asarray(self.X.sum(axis=1)).flatten()

    def generator(
            self,
            idx: Union[np.ndarray, None] = None,
            batch_size: int = 1,
            obs_keys: List[str] = [],
            return_dense: bool = True,
            randomized_batch_access: bool = False,
    ) -> iter:
        """
        Yields an unbiased generator over observations in the contained data sets.

        :param idx: Global idx to query from store. These is an array with indicies corresponding to a contiuous index
            along all observations in self.adatas, ordered along a hypothetical concatenation along the keys of
            self.adatas.
        :param batch_size: Number of observations read from disk in each batched access (generator invocation).
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adatas.
        :param return_dense: Whether to force return count data .X as dense batches. This allows more efficient feature
            indexing if the store is sparse (column indexing on csr matrices is slow).
        :param randomized_batch_access: Whether to randomize batches during reading (in generator). Lifts necessity of
            using a shuffle buffer on generator, however, batch composition stays unchanged over epochs unless there
            is overhangs in retrieval_batch_size in the raw data files, which often happens and results in modest
            changes in batch composition.
        :return: Generator function which yields batch_size at every invocation.
            The generator returns a tuple of (.X, .obs) with types:

                - if store format is h5ad: (Union[scipy.sparse.csr_matrix, np.ndarray], pandas.DataFrame)
        """
        idx, var_idx = self._generator_helper(idx=idx)

        def generator():
            # Can treat full data set as a single array because dask keeps expression data out of memory.
            x = self.X[idx, :]
            obs = self.obs.iloc[idx, :]
            n_obs = x.shape[0]
            remainder = n_obs % batch_size
            assert n_obs == obs.shape[0]
            batch_starts_ends = [
                (int(x * batch_size), int(np.minimum((x * batch_size) + batch_size, n_obs)))
                for x in np.arange(0, n_obs // batch_size + int(remainder > 0))
            ]
            batch_range = np.arange(0, len(batch_starts_ends))
            if randomized_batch_access:
                np.random.shuffle(batch_range)
            for i in batch_range:
                s, e = batch_starts_ends[i]
                x_i = x[s:e, :]
                obs_i = obs[obs_keys].iloc[s:e, :]
                if var_idx is not None:
                    x_i = x_i[:, var_idx]
                yield x_i, obs_i

        return generator


def DistributedStore(cache_path: Union[str, os.PathLike], store_format: str = "zarr") -> \
        Union[DistributedStoreH5ad, DistributedStoreZarr]:
    """
    Instantiates a distributed store class.

    Note: this function mimics a class constructor, therefore the upper-case usage in the name.
    This function effectively serves as a conditional constructor.

    :param cache_path: Store directory.
    :param store_format: Format of store {"h5ad", "zarr"}.

        - "h5ad": Returns instance of DistributedStoreH5ad.
        - "zarr": Returns instance of DistributedStoreZarr.
    :return: Instances of a distributed store class.
    """
    if store_format == "h5ad":
        return DistributedStoreH5ad(cache_path=cache_path)
    elif store_format == "zarr":
        return DistributedStoreZarr(cache_path=cache_path)
    else:
        raise ValueError(f"Did not recognize store_format {store_format}.")
