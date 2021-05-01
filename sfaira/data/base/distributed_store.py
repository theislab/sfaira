import anndata
import numpy as np
import os
import pandas as pd
import pickle
import scipy.sparse
from typing import Dict, List, Union

from sfaira.consts import AdataIdsSfaira, OCS
from sfaira.data.base.dataset import is_child, UNS_STRING_META_IN_OBS
from sfaira.versions.genomes import GenomeContainer
from sfaira.versions.metadata import CelltypeUniverse


class DistributedStore:
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

    adatas: Dict[str, anndata.AnnData]
    indices: Dict[str, np.ndarray]

    def __init__(self, cache_path: Union[str, os.PathLike, None] = None):
        """
        This class is instantiated on a cache directory which contains pre-processed files in rapid access format.

        Supported and automatically identifed are the formats:

            - h5ad,
            - zarr

        :param cache_path: Directory in which pre-processed .h5ad files lie.
        :param genome_container: GenomeContainer with target features space defined.
        """
        # Collect all data loaders from files in directory:
        adatas = {}
        indices = {}
        for f in os.listdir(cache_path):
            if os.path.isfile(os.path.join(cache_path, f)):  # only files
                # Narrow down to supported file types:
                if f.split(".")[-1] == "h5ad":
                    adata = anndata.read_h5ad(
                        filename=os.path.join(cache_path, f),
                        backed="r",
                    )
                elif f.split(".")[-1] == "zarr":
                    # TODO this reads into memory! Might need to directly interface the zarr arrays to work with dask.
                    adata = anndata.read_zarr(os.path.join(cache_path, f))
                else:
                    adata = None
                if adata is not None:
                    adatas[adata.uns["id"]] = adata
                    indices[adata.uns["id"]] = np.arange(0, adata.n_obs)
        self.adatas = adatas
        self.indices = indices
        self.ontology_container = OCS
        self._genome_container = None
        self._adata_ids_sfaira = AdataIdsSfaira()
        self._celltype_universe = None

    @property
    def adata(self):
        return self.adatas[list(self.adatas.keys())[0]].concatenate(
            *[self.adatas[k] for k in list(self.adatas.keys())[1:]],
            batch_key="dataset_id",
            batch_categories=list(self.adatas.keys()),
        )

    @property
    def genome_container(self) -> Union[GenomeContainer, None]:
        return self._genome_container

    @genome_container.setter
    def genome_container(self, x: GenomeContainer):
        var_names = self.__validate_feature_space_homogeneity()
        # Validate genome container choice:
        # Make sure that all var names defined in genome container are also contained in loaded data sets.
        assert np.all([y in var_names for y in x.ensembl]), \
            "did not find variable names from genome container in store"
        self._genome_container = x

    def __validate_feature_space_homogeneity(self) -> List[str]:
        """
        Assert that the data sets which were kept have the same feature names.
        """
        var_names = self.adatas[list(self.adatas.keys())[0]].var_names.tolist()
        for k, v in self.adatas.items():
            assert len(var_names) == len(v.var_names), f"number of features in store differed in object {k}"
            assert np.all(var_names == v.var_names), f"var_names in store were not matched in object {k}"
        return var_names

    def generator(
            self,
            idx: Union[np.ndarray, None] = None,
            batch_size: int = 1,
            obs_keys: List[str] = [],
            return_dense: bool = True,
    ) -> iter:
        """
        Yields an unbiased generator over observations in the contained data sets.

        :param idx: Global idx to query from store. These is an array with indicies corresponding to a contiuous index
            along all observations in self.adatas, ordered along a hypothetical concatenation along the keys of
            self.adatas.
        :param batch_size: Number of observations in each batch (generator invocation). Increasing this may result in
            large speed-ups in query time but removes the ability of upstream generators to fully shuffle cells, as
            these batches are the smallest data unit that upstream generators can access.
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adatas.
        :param return_dense: Whether to return count data .X as dense batches. This allows more efficient feature
            indexing if the store is sparse (column indexing on csr matrices is slow).
        :return: Generator function which yields batch_size at every invocation.
            The generator returns a tuple of (.X, .obs) with types:

                - if store format is h5ad: (Union[scipy.sparse.csr_matrix, np.ndarray], pandas.DataFrame)
        """
        # Make sure that features are ordered in the same way in each object so that generator yields consistent cell
        # vectors.
        _ = self.__validate_feature_space_homogeneity()
        var_names_store = self.adatas[list(self.adatas.keys())[0]].var_names.tolist()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is not None:
            var_names_target = self.genome_container.ensembl
            var_idx = np.sort([var_names_store.index(x) for x in var_names_target])
            # Check if index vector is just full ordered list of indices, in this case, sub-setting is unnecessary.
            if len(var_idx) == len(var_names_store) and np.all(var_idx == np.arange(0, len(var_names_store))):
                var_idx = None
        else:
            var_idx = None

        def generator() -> tuple:
            global_index_set = self.indices_global
            for i, (k, v) in enumerate(self.adatas.items()):
                # Define batch partitions:
                # Get subset of target indices that fall into this data set.
                # Use indices relative to this data (via .index here).
                # continuous_slices is evaluated to establish whether slicing can be performed as the potentially
                # faster [start:end] or needs to tbe index wise [indices]
                if idx is not None:
                    idx_i = [global_index_set[k].tolist().index(x) for x in idx if x in global_index_set[k]]
                    idx_i = np.sort(idx_i)
                    continuous_slices = np.all(idx_i == np.arange(0, v.n_obs))
                else:
                    idx_i = np.arange(0, v.n_obs)
                    continuous_slices = True
                if len(idx_i) > 0:  # Skip data objects without matched cells.
                    n_obs = len(idx_i)
                    # Cells left over after batching to batch size, accounting for overhang:
                    remainder = n_obs % batch_size
                    batch_starts_ends = [
                        (int(x * batch_size), int(x * batch_size) + batch_size)
                        for x in np.arange(0, n_obs // batch_size + int(remainder > 0))
                    ]
                    # Iterate over batches:
                    for j, (s, e) in enumerate(batch_starts_ends):
                        if continuous_slices:
                            e = idx_i[e] if e < n_obs else n_obs
                            x = v.X[idx_i[s]:e, :]
                        else:
                            x = v.X[idx_i[s:e], :]
                        # Do dense conversion now so that col-wise indexing is not slow, often, dense conversion
                        # would be done later anyway.
                        if return_dense:
                            x = x.todense()
                        if var_idx is not None:
                            x = x[:, var_idx]
                        if continuous_slices:
                            e = idx_i[e] if e < n_obs else n_obs
                            obs = v.obs[obs_keys].iloc[idx_i[s]:e, :]
                        else:
                            obs = v.obs[obs_keys].iloc[idx_i[s:e], :]
                        assert isinstance(obs, pd.DataFrame), f"{type(obs)}"
                        # Yield current batch.
                        yield x, obs

        return generator

    @property
    def celltypes_universe(self) -> CelltypeUniverse:
        if self._celltype_universe is None:
            self._celltype_universe = CelltypeUniverse(
                cl=self.ontology_container.cellontology_class,
                uberon=self.ontology_container.organ,
                organism=None,  # TODO Does not load extensions!
            )
        return self._celltype_universe

    def get_subset_idx(self, attr_key, values: Union[str, List[str]]) -> dict:
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
        :param values: Classes to overlap to.
        :return dictionary of files and observation indices by file.
        """
        if not isinstance(values, list):
            values = [values]

        def get_idx(adata, k, dataset):
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
                    raise ValueError(f"did not find unique attribute {k} in data set {dataset}")
            values_found_unique = np.unique(values_found)
            try:
                ontology = getattr(self.ontology_container, k)
            except AttributeError:
                raise ValueError(f"{k} not a valid property of ontology_container object")
            # Test only unique elements found in ontology to save time.
            values_found_unique_matched = [
                x for x in values_found_unique if np.any([
                    is_child(query=x, ontology=ontology, ontology_parent=y)
                    for y in values
                ])
            ]
            # TODO keep this logging for now to catch undesired behaviour resulting from loaded edges in ontologies.
            print(f"matched cell-wise keys {str(values_found_unique_matched)} in data set {dataset}")
            idx = np.where([x in values_found_unique_matched for x in values_found])[0]
            return idx

        indices = {}
        for k, v in self.indices.items():
            idx_old = v.tolist()
            if k not in self.adatas.keys():
                raise ValueError(f"data set {k} queried by indices does not exist in store (.adatas)")
            idx_new = get_idx(adata=self.adatas[k], k=attr_key, dataset=k)
            # Keep intersection of old and new hits.
            idx_new = list(set(idx_old).intersection(set(idx_new)))
            if len(idx_new) > 0:
                indices[k] = np.asarray(idx_new, dtype="int32")
        return indices

    def subset(self, attr_key, values: Union[str, List[str]]):
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
        :param values: Classes to overlap to.
        """
        self.indices = self.get_subset_idx(attr_key=attr_key, values=values)

        for k in list(self.adatas.keys()):
            if k not in self.indices or self.indices[k].shape[0] == 0:  # No observations (cells) left.
                del self.adatas[k]

    def get_subset_idx_global(self, attr_key, values: Union[str, List[str]]) -> np.ndarray:
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
        idx_by_dataset = self.get_subset_idx(attr_key=attr_key, values=values)
        # Translate file-wise indices into global index list across all data sets.
        idx = []
        counter = 0
        for k, v in idx_by_dataset.items():
            idx.extend((v + counter).tolist())
            counter += self.adatas[k].n_obs
        return np.asarray(idx)

    @property
    def indices_global(self) -> dict:
        """
        Increasing indices across data sets which can be concatenated into a single index vector with unique entries
        for cells.
        """
        counter = 0
        indices = {}
        for k, v in self.indices.items():
            indices[k] = np.arange(counter, counter + len(v))
            counter += len(v)
        return indices

    def global_indices_to_dict(self, idx: np.ndarray) -> dict:
        assert np.max(idx) < self.n_obs
        indices = {}
        for k, v in self.indices_global.items():
            matched = [x in v for x in idx]
            if np.any(matched):
                indices[k] = np.array([v.tolist().index(x) for x in idx[matched]])[0]
        return indices

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
        var_names = self.__validate_feature_space_homogeneity()
        # Use feature space sub-selection based on assembly if provided, will use full feature space otherwise.
        if self.genome_container is None:
            return var_names
        else:
            return self.genome_container.ensembl

    @property
    def n_vars(self):
        var_names = self.__validate_feature_space_homogeneity()
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
    def obs(self) -> pd.DataFrame:
        """
        Assemble .obs table of subset of full data.

        :return: .obs data frame.
        """
        return pd.concat([self.adatas[k].obs.iloc[v, :] for k, v in self.indices.items()], axis=0)

    def n_counts(self, idx: Union[np.ndarray, None] = None) -> np.ndarray:
        """
        Compute sum over features for each observation in index.

        :param idx: Subset of index to run operation over. Can be either global index (list) or dictionary of data set
            wise indices.
        :return: Array with sum per observations: (number of observations in index,)
        """
        idx = self.indices if idx is None else idx
        if isinstance(idx, list) or isinstance(idx, np.ndarray):
            # Convert global index to dictionary of indices:
            idx = self.global_indices_to_dict(idx=idx)
        return np.concatenate([
            np.asarray(self.adatas[k].X[v, :].sum(axis=1)).flatten()
            for k, v in idx.items()
        ], axis=0)
