import anndata
import numpy as np
import os
import pandas as pd
import pickle
import scipy.sparse
from typing import Dict, List, Union

from sfaira.consts import AdataIdsSfaira, OCS
from sfaira.data.base.dataset import is_child
from sfaira.versions.metadata import CelltypeUniverse


class DistributedStore:
    """
    Data set group class tailored to data access requirements common in high-performance computing (HPC).

    This class does not inherit from DatasetGroup because it entirely relies on the cached objects.
    """

    indices: Dict[str, np.ndarray]

    def __init__(self, cache_path: Union[str, None] = None):
        """
        This class is instantiated on a cache directory which contains pre-processed files in rapid access format.

        Supported and automatically identifed are the formats:

            - h5ad,
            - zarr

        :param cache_path: Directory in which pre-processed .h5ad files lie.
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
                        backed=True,
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
        self._adata_ids_sfaira = AdataIdsSfaira()
        self._celltype_universe = None

    @property
    def adata(self):
        return list(self.adatas.values)[0].concatenate(
            list(self.adatas.values)[1:]
        )

    def generator(
            self,
            batch_size: int = 1,
            obs_keys: List[str] = [],
            continuous_batches: bool = True,
    ) -> iter:
        """
        Yields an unbiased generator over observations in the contained data sets.

        :param batch_size: Number of observations in each batch (generator invocation).
        :param obs_keys: .obs columns to return in the generator. These have to be a subset of the columns available
            in self.adatas.
        :param continuous_batches: Whether to build batches of batch_size across data set boundaries if end of one
            data set is reached.
        :return: Generator function which yields batch_size at every invocation.
            The generator returns a tuple of (.X, .obs) with types:

                - if store format is h5ad: (scipy.sparse.csr_matrix, pandas.DataFrame)
        """

        def generator() -> tuple:
            n_datasets = len(list(self.adatas.keys()))
            x_last = None
            obs_last = None
            for i, (k, v) in enumerate(self.adatas.items()):
                # Define batch partitions:
                if continuous_batches and x_last is not None:
                    # Prepend data set with residual data from last data set.
                    remainder_start = x_last.shape[0]
                    n_obs = v.n_obs + remainder_start
                else:
                    # Partition into equally sized batches up to last batch.
                    remainder_start = 0
                    n_obs = v.n_obs
                remainder = n_obs % batch_size
                batch_starts = [
                    np.min([0, int(x * batch_size - remainder_start)])
                    for x in np.arange(1, n_obs // batch_size + int(remainder > 0))
                ]
                n_batches = len(batch_starts)
                # Iterate over batches:
                for j, x in enumerate(batch_starts):
                    batch_end = int(x + batch_size)
                    x = v.X[x:batch_end, :]
                    obs = v.obs[obs_keys].iloc[x:batch_end, :]
                    assert isinstance(x, scipy.sparse.csr_matrix), f"{type(x)}"
                    assert isinstance(obs, pd.DataFrame), f"{type(obs)}"
                    if continuous_batches and remainder > 0 and i < (n_datasets - 1) and j == (n_batches - 1):
                        # Cache incomplete last batch to append to next first batch of next data set.
                        x_last = x
                        obs_last = obs
                    elif continuous_batches and x_last is not None:
                        # Append last incomplete batch current batch.
                        x = scipy.sparse.hstack(blocks=[x_last, x], format="csr")
                        obs = pd.concat(objs=[obs_last, obs], axis=0)
                        yield x, obs
                    else:
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

    def subset(self, attr_key, values):
        """
        Subset list of adata objects based on match to values in key property.

        Keys need to be available in adata.uns

        :param attr_key: Property to subset by.
        :param values: Classes to overlap to.
        :return:
        """
        if isinstance(values, np.ndarray):
            values = values.tolist()
        if isinstance(values, tuple):
            values = list(values)
        if not isinstance(values, list):
            values = [values]
        # Get ontology container to be able to do relational reasoning:
        ontology = getattr(self.ontology_container, attr_key)
        for k in list(self.adatas.keys()):
            if getattr(self._adata_ids_sfaira, attr_key) in self.adatas[k].uns.keys():
                values_found = self.adatas[k].uns[getattr(self._adata_ids_sfaira, attr_key)]
                if not isinstance(values_found, list):
                    values_found = [values_found]
                if not np.any([
                    np.any([
                        is_child(query=x, ontology=ontology, ontology_parent=y)
                        for y in values
                    ]) for x in values_found
                ]):
                    # Delete entries which a non-matching meta data value associated with this item.
                    del self.adatas[k]
            else:
                # Delete entries which did not have this key annotated.
                del self.adatas[k]

    def subset_cells_idx(self, attr_key, values: Union[str, List[str]]):
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

        def get_subset_idx(adata, k):
            values_found = adata.obs[getattr(self._adata_ids_sfaira, k)].values
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
            print(f"matched cell-wise keys {str(values_found_unique_matched)} in data set {self.id}")
            idx = np.where([x in values_found_unique_matched for x in values_found])[0]
            return idx

        indices = {}
        for k, v in self.adatas.items():
            idx_old = self.indices[k].tolist()
            idx_new = get_subset_idx(adata=v, k=attr_key)
            # Keep intersection of old and new hits.
            indices[k] = np.array(list(set(idx_old).intersection(set(idx_new))))
        return indices

    def subset_cells(self, attr_key, values: Union[str, List[str]]):
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
        self.indices = self.subset_cells_idx(attr_key=attr_key, values=values)

        for k, v in self.indices.items():
            if v.shape[0] == 0:  # No observations (cells) left.
                del self.adatas[k]

    def subset_cells_idx_global(self, attr_key, values: Union[str, List[str]]):
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
        idx_by_dataset = self.subset_cells_idx(attr_key=attr_key, values=values)
        # Translate file-wise indices into global index list across all data sets.
        idx = []
        counter = 0
        for k, v in self.adatas.items():
            idx_k = np.arange(counter, counter + v.n_obs)
            idx.extend(idx_k[idx_by_dataset[k]])
            counter += v.n_obs
        return idx

    def write_config(self, fn: Union[str, os.PathLike]):
        """
        Writes a config file that describes the current data sub-setting.

        This config file can be loaded later to recreate a sub-setting.
        This config file contains observation-wise subsetting information.

        :param fn: Output file without file type extension.
        """
        with open(fn + '.pickle', 'w') as f:
            pickle.dump(self.indices, f)

    def load_config(self, fn: Union[str, os.PathLike]):
        """
        Load a config file and recreates a data sub-setting.
        This config file contains observation-wise subsetting information.

        :param fn: Output file without file type extension.
        """
        with open(fn + '.pickle', 'rb') as f:
            self.indices = pickle.load(f)
        # Subset to described data sets:
        for x in self.indices.keys():
            if x not in self.adatas.keys():
                raise ValueError(f"did not find object with name {x} in currently loaded universe")
        # Only retain data sets with which are mentioned in config file.
        self.subset(attr_key="id", values=list(self.indices.keys()))

    @property
    def n_vars(self):
        # assumes that all adata
        return list(self.adatas.values())[0].n_vars

    @property
    def n_obs(self):
        return np.sum([len(v) for _, v in self.indices])
