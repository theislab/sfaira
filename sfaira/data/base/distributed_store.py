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
from typing import Dict, List, Tuple, Union

from sfaira.consts import AdataIdsSfaira, OCS
from sfaira.data.base.dataset import is_child, UNS_STRING_META_IN_OBS
from sfaira.data.base.io_dao import read_dao
from sfaira.versions.genomes.genomes import GenomeContainer

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
    This class is centred around .adata_by_key and .indices.

    .adata_by_key is a dictionary (by id) of backed anndata instances that point to individual h5ads.
    This dictionary is intialised with all h5ads in the store.
    As the store is subsetted, key-value pairs are deleted from this dictionary.

    .indices have keys that correspond to keys in .adata_by_key and contain index vectors of observations in the anndata
    instances in .adata_by_key which are still kept.
    These index vectors are a form of lazy slicing that does not require data set loading or re-writing.
    As the store is subsetted, key-value pairs are deleted from this dictionary if no observations from a given key
    match the subsetting.
    If a subset of observations from a key matches the subsetting operation, the index set in the corresponding value is
    reduced.
    All data retrievel operations work on .indices: Generators run over these indices when retrieving observations for
    example.
    """

    _adata_by_key: Dict[str, anndata.AnnData]
    _indices: Dict[str, np.ndarray]
    _obs_by_key: Union[None, Dict[str, dask.dataframe.DataFrame]]

    def __init__(self, adata_by_key: Dict[str, anndata.AnnData], indices: Dict[str, np.ndarray],
                 obs_by_key: Union[None, Dict[str, dask.dataframe.DataFrame]] = None):
        self.adata_by_key = adata_by_key
        self.indices = indices
        self.obs_by_key = obs_by_key
        self.ontology_container = OCS
        self._genome_container = None
        self._adata_ids_sfaira = AdataIdsSfaira()
        self._celltype_universe = None

    @property
    def idx_by_organism(self) -> Dict[str, np.ndarray]:
        idx_dict = {}
        organisms_by_key = self.organisms_by_key
        for x in self.organisms:
            idx_global_in_organism = np.where(np.concatenate([
                np.ones_like(v) == 1 if organisms_by_key[k] == x else np.ones_like(v) == 0
                for k, v in self.indices.items()
            ]))[0]
            idx_dict[x] = idx_global_in_organism
        return idx_dict

    def _validate_idx(self, idx: Union[np.ndarray, list]) -> Dict[str, np.ndarray]:
        """
        Validate global index vector and transform into a dictionary over organisms.
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
        idx_dict = self.idx_by_organism
        for k, v in idx_dict.items():
            idx_dict[k] = np.asarray(list(set(v).intersection(set(idx))))
        # Sanity check that all indices were assigned uniquely to an organism:
        # TODO this check can be removed in the future or refactored into a unit test:
        assert len(idx) == np.sum([len(x) for x in idx_dict.values()])
        for i, x in enumerate(idx_dict.values()):
            for j, y in enumerate(list(idx_dict.values())[(i + 1):]):
                assert len(set(x).intersection(set(y))) == 0, f"overlap between {self.organisms[i]} and " \
                                                              f"{self.organisms[i + j]}:\n{x}\n{y}"
        return idx_dict

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
    def organisms(self):
        """
        Organisms in store.
        """
        return np.sort(np.unique(list(self.organisms_by_key.values())))

    def _validate_feature_space_homogeneity(self) -> Dict[str, List[str]]:
        """
        Assert that the data sets which were kept have the same feature names within each organism.

        :return: List of feature names in shared feature space or dictionary of list of features by organism.
        """
        organisms = self.organisms_by_key
        organisms_unique = self.organisms
        var_names = {}
        for x in organisms_unique:
            ks_x = [y for y, z in organisms.items() if z == x]
            var_names_x = self._adata_by_key[ks_x[0]].var_names.tolist()
            for k in ks_x:
                assert len(var_names_x) == len(self._adata_by_key[k].var_names), \
                    f"number of features in store differed in object {k} compared to {ks_x[0]}"
                assert np.all(var_names_x == self._adata_by_key[k].var_names), \
                    f"var_names in store were not matched in object {k} compared to {ks_x[0]}"
            var_names[x] = var_names_x
        # Return as an organism-wise dictionary of feature names.
        return var_names

    def _generator_helper(
            self,
            idx: Union[np.ndarray, None] = None,
    ) -> Tuple[Union[Dict[str, np.ndarray], None], Union[Dict[str, np.ndarray], None]]:
        # Make sure that features are ordered in the same way in each object so that generator yields consistent cell
        # vectors.
        var_names = self._validate_feature_space_homogeneity()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is not None:
            for k in self.genome_container.keys():
                assert k in var_names.keys(), f"did not match {k} from var_names in GenomeContainer " \
                                              f"{self.genome_container.keys()}"
            var_idx = {}
            for k, v in self.genome_container.items():
                var_names_target = v.ensembl
                var_idx[k] = np.sort([var_names[k].index(x) for x in var_names_target])
                # Check if index vector is just full ordered list of indices, in this case, sub-setting is unnecessary.
                if len(var_idx[k]) == len(var_names[k]) and np.all(var_idx[k] == np.arange(0, len(var_names[k]))):
                    var_idx[k] = None
        else:
            var_idx = None
        if idx is not None:
            idx = self._validate_idx(idx)
        return idx, var_idx

    @property
    def adata_by_key(self) -> Dict[str, anndata.AnnData]:
        return self._adata_by_key

    @adata_by_key.setter
    def adata_by_key(self, x: Dict[str, anndata.AnnData]):
        self._adata_by_key = x

    def adata_memory_footprint(self, k):
        """
        Memory foot-print of data set k in MB.
        """
        return sys.getsizeof(self.adata_by_key[k]) / np.power(1024, 2)

    @property
    def indices(self) -> Dict[str, np.ndarray]:
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
    def idx_dataset_start(self):
        """
        Global indices corresponding to first cell of each data set.
        """
        idx = []
        counter = 0
        for k, v in self.indices.items():
            idx.append(counter)
            counter += len(v)
        return idx

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
    def genome_container(self) -> Union[GenomeContainer, Dict[str, GenomeContainer], None]:
        return self._genome_container

    @genome_container.setter
    def genome_container(self, x: Union[GenomeContainer, Dict[str, GenomeContainer]]):
        if isinstance(x, GenomeContainer):
            # Transform into dictionary first.
            organisms = self.organisms
            if len(organisms) > 1:
                raise ValueError(f"Gave a single GenomeContainer for a store instance that has mulitiple organism: "
                                 f"{organisms}, either further subset the store or give a dictionary of "
                                 f"GenomeContainers")
            else:
                x = {organisms[0]: x}
        var_names = self._validate_feature_space_homogeneity()
        # Validate genome container choice:
        # Make sure that all var names defined in genome container are also contained in loaded data sets.
        for k, v in var_names.items():
            assert k in x.keys(), f"did not find organism {k} which occurs in store in x"
            assert np.all([y in v for y in x[k].ensembl]), \
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
        # Make sure all declared data sets are found in store:
        for x in self.indices.keys():
            if x not in self.adata_by_key.keys():
                raise ValueError(f"did not find object with name {x} in currently loaded universe")

    @property
    def var_names(self) -> Dict[str, List[str]]:
        """
        Feature names of selected genes by organism in store.
        """
        var_names = self._validate_feature_space_homogeneity()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is None:
            return var_names
        else:
            return dict([(k, v.ensembl) for k, v in self.genome_container.items()])

    @property
    def n_vars(self) -> Dict[str, int]:
        """
        Number of selected features per organism in store
        """
        var_names = self._validate_feature_space_homogeneity()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is None:
            return dict([(k, len(v)) for k, v in var_names.items()])
        else:
            return dict([(k, v.n_var) for k, v in self.genome_container.items()])

    @property
    def n_obs(self) -> int:
        """
        Number of observations selected in store.
        """
        return np.sum([len(v) for v in self.indices.values()])

    @property
    def n_obs_organism(self) -> Dict[str, int]:
        """
        Number of observations selected in store per organism as dictionary.
        """
        organisms_by_key = self.organisms_by_key
        return dict([
            (x, np.sum([len(v) for k,v in self.indices.items() if organisms_by_key[k] == x]))
            for x in self.organisms
        ])

    @property
    def shape(self) -> Dict[str, Tuple[int, int]]:
        return dict([(k, (self.n_obs_organism[k], v)) for k, v in self.n_vars.items()])

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
    def X(self) -> Union[dask.array.Array, scipy.sparse.csr_matrix]:
        pass

    @property
    @abc.abstractmethod
    def obs(self) -> Union[pd.DataFrame]:
        pass

    @abc.abstractmethod
    def n_counts(self, idx: Union[np.ndarray, list, None] = None) -> np.ndarray:
        pass


class DistributedStoreH5ad(DistributedStoreBase):

    def __init__(self, cache_path: Union[str, os.PathLike]):
        # Collect all data loaders from files in directory:
        adata_by_key = {}
        indices = {}
        for f in np.sort(os.listdir(cache_path)):
            adata = None
            trial_path = os.path.join(cache_path, f)
            if os.path.isfile(trial_path):
                # Narrow down to supported file types:
                if f.split(".")[-1] == "h5ad":
                    try:
                        adata = anndata.read_h5ad(
                            filename=trial_path,
                            backed="r",
                        )
                    except OSError as e:
                        adata = None
                        print(f"WARNING: for data set {f}: {e}")
            if adata is not None:
                adata_by_key[adata.uns["id"]] = adata
                indices[adata.uns["id"]] = np.arange(0, adata.n_obs)
        self._x_as_dask = False
        super(DistributedStoreH5ad, self).__init__(adata_by_key=adata_by_key, indices=indices)

    @property
    def adata_sliced(self) -> Dict[str, anndata.AnnData]:
        """
        Only exposes the subset and slices of the adata instances contained in ._adata_by_key defined in .indices.
        """
        return dict([(k, self._adata_by_key[k][v, :]) for k, v in self.indices.items()])

    @property
    def X(self):
        assert False

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

    def n_counts(self, idx: Union[np.ndarray, list, None] = None) -> np.ndarray:
        """
        Compute sum over features for each observation in index.

        :param idx: Index vector over observations in object.
        :return: Array with sum per observations: (number of observations in index,)
        """
        return np.concatenate([
            np.asarray(v.X.sum(axis=1)).flatten()
            for v in self.adata_by_key_subset(idx=idx).values()
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
            along all observations in self.adata_by_key, ordered along a hypothetical concatenation along the keys of
            self.adata_by_key.
        :param batch_size: Number of observations read from disk in each batched access (generator invocation).
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adata_by_key.
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
            adatas_sliced_subset = self.adata_by_key_subset(idx=idx)
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
                x, obs = access_helper(adata=adatas_sliced_subset[k], s=s, e=e, j=var_idx[self.organisms_by_key[k]],
                                       return_dense=return_dense, obs_keys=obs_keys)
                yield x, obs

        return generator

    def adata_by_key_subset(self, idx: Union[np.ndarray, list]) -> Dict[str, anndata.AnnData]:
        """
        Subsets adata_by_key as if it was one object, ie behaves the same way as self.adata[idx]  without explicitly
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
            return dict([(k, self._adata_by_key[k][v, :]) for k, v in indices_subsetted.items()])
        else:
            return self.adata_sliced

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
        :param excluded_values: Classes to exclude from match list. Supply either values or excluded_values.
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


class DistributedStoreDao(DistributedStoreBase):

    def __init__(self, cache_path: Union[str, os.PathLike], columns: Union[None, List[str]] = None):
        """

        :param cache_path: Store directory.
        :param columns: Which columns to read into the obs copy in the output, see pandas.read_parquet().
        """
        # Collect all data loaders from files in directory:
        adata_by_key = {}
        indices = {}
        for f in np.sort(os.listdir(cache_path)):
            adata = None
            trial_path = os.path.join(cache_path, f)
            if os.path.isdir(trial_path):
                # zarr-backed anndata are saved as directories with the elements of the array group as further sub
                # directories, e.g. a directory called "X", and a file ".zgroup" which identifies the zarr group.
                adata = read_dao(trial_path, use_dask=True, columns=columns, obs_separate=False)
            if adata is not None:
                adata_by_key[adata.uns["id"]] = adata
                indices[adata.uns["id"]] = np.arange(0, adata.n_obs)
        self._x_as_dask = True
        super(DistributedStoreDao, self).__init__(adata_by_key=adata_by_key, indices=indices, obs_by_key=None)

    @property
    def X(self) -> dask.array.Array:
        """
        One dask array of all cells.

        Requires feature dimension to be shared.
        """
        assert np.all([isinstance(self._adata_by_key[k].X, dask.array.Array) for k in self.indices.keys()])
        return dask.array.vstack([
            self._adata_by_key[k].X[v, :]
            for k, v in self.indices.items()
        ])

    @property
    def X_by_organism(self) -> Dict[str, dask.array.Array]:
        """
        One dask array of all cells per organism in store
        """
        assert np.all([isinstance(self._adata_by_key[k].X, dask.array.Array) for k in self.indices.keys()])
        return dict([
            (organism, dask.array.vstack([
                self._adata_by_key[k].X[v, :]
                for k, v in self.indices.items()
                if self.organisms_by_key[k] == organism
            ])) for organism in self.organisms
        ])

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

    @property
    def obs_by_organism(self) -> Dict[str, pd.DataFrame]:
        """
        Assemble .obs table of subset of selected data per organism.

        Resulting index is increasing integers starting with zero.

        :return: .obs data frame.
        """
        # TODO Using loc indexing here instead of iloc, this might be faster on larger tables?
        return dict([
            (organism, pd.concat([
                self.adata_by_key[k].obs.loc[self.adata_by_key[k].obs.index[v], :]
                for k, v in self.indices.items()
                if self.organisms_by_key[k] == organism
            ], axis=0, join="inner", ignore_index=True, copy=False))
            for organism in self.organisms
        ])

    def n_counts(self, idx: Union[Dict[str, Union[np.ndarray, list]], None] = None) -> Dict[str, np.ndarray]:
        """
        Compute sum over features for each observation in index.

        :param idx: Index vector over observations in object.
        :return: Array with sum per observations per organism: (number of observations in index,)
        """
        if idx is not None:
            return np.sum([np.asarray(x.sum(axis=1)).flatten() for x in self.X.values()])
        else:
            return dict([(x, np.asarray(x.sum(axis=1)).flatten()) for k, v in self.X_by_organism.items()])

    def generator(
            self,
            idx: Union[np.ndarray, None] = None,
            batch_size: int = 1,
            obs_keys: List[str] = [],
            return_dense: bool = True,
            randomized_batch_access: bool = False,
            random_access: bool = False,
    ) -> iter:
        """
        Yields an unbiased generator over observations in the contained data sets.

        :param idx: Global idx to query from store. These is an array with indicies corresponding to a contiuous index
            along all observations in self.adata_by_key, ordered along a hypothetical concatenation along the keys of
            self.adata_by_key.
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
            The generator returns a tuple of (.X, .obs) with types:

                - if store format is h5ad: (Union[scipy.sparse.csr_matrix, np.ndarray], pandas.DataFrame)
        """
        idx_dict, var_idx_dict = self._generator_helper(idx=idx)
        # Normalise cell indices such that each organism is indexed starting at zero:
        # This is required below because each organism is represented as its own dask array.
        idx_dict_organism = self.idx_by_organism
        for k, v in idx_dict.items():
            ref_idx = idx_dict_organism[k].tolist()
            idx_dict[k] = np.array([ref_idx.index(x) for x in v])
        if randomized_batch_access and random_access:
            raise ValueError("Do not use randomized_batch_access and random_access.")
        if batch_size > np.min([len(v) for v in idx_dict.values()]):
            batch_size_new = np.min([len(v) for v in idx_dict.values()])
            print(f"WARNING: reduing retieval batch size according to data availability in store "
                  f"from {batch_size} to {batch_size_new}")
            batch_size = batch_size_new
        x_dict = self.X_by_organism
        obs_dict = self.obs_by_organism

        def generator():
            # Can all data sets corresponding to one organism as a single array because they share the second dimension
            # and dask keeps expression data and obs out of memory.
            for organism in self.organisms:
                idx_o = idx_dict[organism]
                var_idx = var_idx_dict[organism]
                x = x_dict[organism][idx_o, :]
                obs = obs_dict[organism]
                obs = obs.loc[obs.index[idx_o], obs_keys]  # TODO better than iloc?
                # Redefine index so that .loc indexing can be used instead of .iloc indexing:
                obs.index = np.arange(0, obs.shape[0])
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
                epoch_indices = np.arange(0, n_obs)
                if random_access:
                    np.random.shuffle(epoch_indices)
                for i in batch_range:
                    s, e = batch_starts_ends[i]
                    # Feature indexing: Run in same operation as observation index so that feature chunking can be
                    # efficiently used if available. TODO does this make a difference in dask?
                    if random_access:
                        x_i = x[epoch_indices[s:e], :]
                        if var_idx is not None:
                            x_i = x_i[:, var_idx]
                    else:
                        # Use slicing because observations accessed in batch are ordered in data set:
                        # Note that epoch_indices[i] == i if not random_access.
                        x_i = x[s:e, :]
                        if var_idx is not None:
                            x_i = x_i[:, var_idx]
                    # Exploit fact that index of obs is just increasing list of integers, so we can use the .loc[]
                    # indexing instead of .iloc[]:
                    obs_i = obs.loc[epoch_indices[s:e], :]
                    yield x_i, obs_i

        return generator


def load_store(cache_path: Union[str, os.PathLike], store_format: str = "dao",
               columns: Union[None, List[str]] = None) -> Union[DistributedStoreH5ad, DistributedStoreDao]:
    """
    Instantiates a distributed store class.

    :param cache_path: Store directory.
    :param store_format: Format of store {"h5ad", "dao"}.

        - "h5ad": Returns instance of DistributedStoreH5ad.
        - "dao": Returns instance of DistributedStoreDoa (distributed access optimized).
    :param columns: Which columns to read into the obs copy in the output, see pandas.read_parquet().
        Only relevant if store_format is "dao".
    :return: Instances of a distributed store class.
    """
    if store_format == "h5ad":
        return DistributedStoreH5ad(cache_path=cache_path)
    elif store_format == "dao":
        return DistributedStoreDao(cache_path=cache_path, columns=columns)
    else:
        raise ValueError(f"Did not recognize store_format {store_format}.")
