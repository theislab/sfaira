import abc
import anndata
import numpy as np
import os
import pickle
from typing import Dict, List, Tuple, Union

from sfaira.consts import AdataIdsSfaira
from sfaira.data.store.single_store import DistributedStoreSingleFeatureSpace, \
    DistributedStoreDao, DistributedStoreH5ad
from sfaira.data.store.io_dao import read_dao
from sfaira.versions.genomes.genomes import GenomeContainer


class DistributedStoreMultipleFeatureSpaceBase(abc.ABC):

    """
    Umbrella class for a dictionary over multiple instances DistributedStoreSingleFeatureSpace.

    Allows for operations on data sets that are defined in different feature spaces.
    """

    _adata_ids_sfaira: AdataIdsSfaira
    _stores: Dict[str, DistributedStoreSingleFeatureSpace]

    def __init__(self, stores: Dict[str, DistributedStoreSingleFeatureSpace]):
        self._stores = stores

    @property
    def stores(self) -> Dict[str, DistributedStoreSingleFeatureSpace]:
        """
        Only expose stores that contain observations.
        """
        return dict([(k, v) for k, v in self._stores.items() if v.n_obs > 0])

    @stores.setter
    def stores(self, x: Dict[str, DistributedStoreSingleFeatureSpace]):
        raise NotImplementedError("cannot set this attribute, it s defined in constructor")

    @property
    def genome_containers(self) -> Dict[str, Union[GenomeContainer, None]]:
        return dict([(k, v.genome_container) for k, v in self._stores.items()])

    @genome_containers.setter
    def genome_containers(self, x: Union[GenomeContainer, Dict[str, GenomeContainer]]):
        if isinstance(x, GenomeContainer):
            # Transform into dictionary first.
            organisms = [k for k, v in self.stores.items()]
            if isinstance(organisms, list) and len(organisms) == 0:
                raise Warning("found empty organism lists in genome_container.setter")
            if len(organisms) > 1:
                raise ValueError(f"Gave a single GenomeContainer for a store instance that has mulitiple organism: "
                                 f"{organisms}, either further subset the store or give a dictionary of "
                                 f"GenomeContainers")
            else:
                x = {organisms[0]: x}
        for k, v in x.items():
            self.stores[k].genome_container = v

    @property
    def indices(self) -> Dict[str, np.ndarray]:
        """
        Dictionary of indices of selected observations contained in all stores.
        """
        return dict([(kk, vv) for k, v in self.stores.items() for kk, vv in v.indices.items()])

    @property
    def adata_by_key(self) -> Dict[str, anndata.AnnData]:
        """
        Dictionary of all anndata instances for each selected data set in store, sub-setted by selected cells, for each
        stores.
        """
        return dict([(kk, vv) for k, v in self.stores.items() for kk, vv in v.adata_by_key.items()])

    @property
    def data_by_key(self):
        """
        Data matrix for each selected data set in store, sub-setted by selected cells.
        """
        return dict([(kk, vv) for k, v in self.stores.items() for kk, vv in v.data_by_key.items()])

    @property
    def var_names(self) -> Dict[str, List[str]]:
        """
        Dictionary of variable names by store.
        """
        return dict([(k, v.var_names) for k, v in self.stores.items()])

    @property
    def n_vars(self) -> Dict[str, int]:
        """
        Dictionary of number of features by store.
        """
        return dict([(k, v.n_vars) for k, v in self.stores.items()])

    @property
    def n_obs(self) -> Dict[str, int]:
        """
        Dictionary of number of observations by store.
        """
        return dict([(k, v.n_obs) for k, v in self.stores.items()])

    @property
    def obs(self):
        """
        Dictionary of concatenated .obs tables by store, only including non-empty stores.
        """
        return dict([(k, v.obs) for k, v in self.stores.items()])

    @property
    def X(self):
        """
        Dictionary of concatenated data matrices by store, only including non-empty stores.
        """
        return dict([(k, v.X) for k, v in self.stores.items()])

    @property
    def shape(self) -> Dict[str, Tuple[int, int]]:
        """
        Dictionary of full selected data matrix shape by store.
        """
        return dict([(k, v.shape) for k, v in self.stores.items()])

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
        for k in self.stores.keys():
            self.stores[k].subset(attr_key=attr_key, values=values, excluded_values=excluded_values, verbose=0)
        if self.n_obs == 0 and verbose > 0:
            print("WARNING: multi store is now empty.")

    def write_config(self, fn: Union[str, os.PathLike]):
        """
        Writes a config file that describes the current data sub-setting.

        This config file can be loaded later to recreate a sub-setting.
        This config file contains observation-wise subsetting information.

        :param fn: Output file without file type extension.
        """
        indices = {}
        for v in self.stores.values():
            indices.update(v.indices)
        with open(fn + '.pickle', 'wb') as f:
            pickle.dump(indices, f)

    def load_config(self, fn: Union[str, os.PathLike]):
        """
        Load a config file and recreates a data sub-setting.
        This config file contains observation-wise subsetting information.

        :param fn: Output file without file type extension.
        """
        with open(fn, 'rb') as f:
            indices = pickle.load(f)
        # Distribute indices to corresponding stores by matched keys.
        keys_not_found = list(indices.keys())
        for k, v in self.stores.items():
            indices_k = {}
            for i, (kk, vv) in enumerate(indices.items()):
                if kk in v.adata_by_key.keys():
                    indices_k[kk] = vv
                    del keys_not_found[i]
            self.stores[k].indices = indices_k
        # Make sure all declared data were assigned to stores:
        if len(keys_not_found) > 0:
            raise ValueError(f"did not find object(s) with name(s) in store: {keys_not_found}")

    def generator(
            self,
            idx: Union[Dict[str, Union[np.ndarray, None]], None] = None,
            intercalated: bool = True,
            **kwargs
    ) -> Tuple[iter, int]:
        """
        Emission of batches from unbiased generators of all stores.

        See also DistributedStore*.generator().

        :param idx:
        :param intercalated: Whether to do sequential or intercalated emission.
        :param kwargs: See parameters of DistributedStore*.generator().
        """
        if idx is None:
            idx = dict([(k, None) for k in self.stores.keys()])
        for k in self.stores.keys():
            assert k in idx.keys(), (idx.keys(), self.stores.keys())
        generators = [
            v.generator(idx=idx[k], **kwargs)
            for k, v in self.stores.items()
        ]
        generator_fns = [x[0]() for x in generators]
        generator_len = [x[1] for x in generators]

        if intercalated:
            # Define relative drawing frequencies from iterators for intercalation.
            ratio = np.asarray(np.round(np.max(generator_len) / np.asarray(generator_len), 0), dtype="int64")

            def generator():
                # Document which generators are still yielding batches:
                yielding = np.ones((ratio.shape[0],)) == 1.
                while np.any(yielding):
                    # Loop over one iterator length adjusted cycle of emissions.
                    for i, (g, n) in enumerate(zip(generator_fns, ratio)):
                        for _ in range(n):
                            try:
                                x = next(g)
                                yield x
                            except StopIteration:
                                yielding[i] = False
        else:
            def generator():
                for g in generator_fns:
                    for x in g():
                        yield x

        return generator, int(np.sum(generator_len))


class DistributedStoresDao(DistributedStoreMultipleFeatureSpaceBase):

    _dataset_weights: Union[None, Dict[str, float]]

    def __init__(self, cache_path: Union[str, os.PathLike], columns: Union[None, List[str]] = None):
        """

        :param cache_path: Store directory.
        :param columns: Which columns to read into the obs copy in the output, see pandas.read_parquet().
        """
        # Collect all data loaders from files in directory:
        self._adata_ids_sfaira = AdataIdsSfaira()
        adata_by_key = {}
        x_by_key = {}
        indices = {}
        for f in np.sort(os.listdir(cache_path)):
            adata = None
            x = None
            trial_path = os.path.join(cache_path, f)
            if os.path.isdir(trial_path):
                # zarr-backed anndata are saved as directories with the elements of the array group as further sub
                # directories, e.g. a directory called "X", and a file ".zgroup" which identifies the zarr group.
                adata, x = read_dao(trial_path, use_dask=True, columns=columns, obs_separate=False, x_separate=True)
            if adata is not None:
                organism = adata.uns[self._adata_ids_sfaira.organism]
                if organism not in adata_by_key.keys():
                    adata_by_key[organism] = {}
                    x_by_key[organism] = {}
                    indices[organism] = {}
                adata_by_key[organism][adata.uns["id"]] = adata
                x_by_key[organism][adata.uns["id"]] = x
                indices[organism][adata.uns["id"]] = np.arange(0, adata.n_obs)
        self._x_by_key = x_by_key
        stores = dict([
            (k, DistributedStoreDao(adata_by_key=adata_by_key[k], x_by_key=x_by_key[k], indices=indices[k],
                                    obs_by_key=None))
            for k in adata_by_key.keys()
        ])
        super(DistributedStoresDao, self).__init__(stores=stores)


class DistributedStoresH5ad(DistributedStoreMultipleFeatureSpaceBase):

    def __init__(self, cache_path: Union[str, os.PathLike], in_memory: bool = False):
        # Collect all data loaders from files in directory:
        self._adata_ids_sfaira = AdataIdsSfaira()
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
                            backed="r" if in_memory else None,
                        )
                    except OSError as e:
                        adata = None
                        print(f"WARNING: for data set {f}: {e}")
            if adata is not None:
                organism = adata.uns[self._adata_ids_sfaira.organism]
                if organism not in adata_by_key.keys():
                    adata_by_key[organism] = {}
                    indices[organism] = {}
                adata_by_key[organism][adata.uns["id"]] = adata
                indices[organism][adata.uns["id"]] = np.arange(0, adata.n_obs)
        stores = dict([
            (k, DistributedStoreH5ad(adata_by_key=adata_by_key[k], indices=indices[k], in_memory=in_memory))
            for k in adata_by_key.keys()
        ])
        super(DistributedStoresH5ad, self).__init__(stores=stores)


def load_store(cache_path: Union[str, os.PathLike], store_format: str = "dao",
               columns: Union[None, List[str]] = None) -> Union[DistributedStoresH5ad, DistributedStoresDao]:
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
    if store_format == "anndata":
        return DistributedStoresH5ad(cache_path=cache_path, in_memory=True)
    elif store_format == "dao":
        return DistributedStoresDao(cache_path=cache_path, columns=columns)
    elif store_format == "h5ad":
        return DistributedStoresH5ad(cache_path=cache_path, in_memory=False)
    else:
        raise ValueError(f"Did not recognize store_format {store_format}.")
