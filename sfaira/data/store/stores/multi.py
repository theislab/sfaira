import anndata
import dask.dataframe
import numpy as np
import os
import pandas as pd
import pickle
from typing import Dict, List, Tuple, Union

from sfaira.consts import AdataIdsSfaira
from sfaira.data.store.stores.base import StoreBase
from sfaira.data.store.stores.single import StoreSingleFeatureSpace, \
    StoreDao, StoreAnndata
from sfaira.data.store.carts.multi import CartMulti
from sfaira.data.store.io.io_dao import read_dao
from sfaira.versions.genomes.genomes import GenomeContainer


class StoreMultipleFeatureSpaceBase(StoreBase):

    """
    Umbrella class for a dictionary over multiple instances DistributedStoreSingleFeatureSpace.

    Allows for operations on data sets that are defined in different feature spaces.
    """

    _adata_ids_sfaira: AdataIdsSfaira
    _stores: Dict[str, StoreSingleFeatureSpace]

    def __init__(self, stores: Dict[str, StoreSingleFeatureSpace]):
        self._stores = stores

    @property
    def stores(self) -> Dict[str, StoreSingleFeatureSpace]:
        """
        Only expose stores that contain observations.
        """
        return dict([(k, v) for k, v in self._stores.items() if v.n_obs > 0])

    @stores.setter
    def stores(self, x: Dict[str, StoreSingleFeatureSpace]):
        raise NotImplementedError("cannot set this attribute, it s defined in constructor")

    @property
    def genome_container(self) -> Dict[str, Union[GenomeContainer, None]]:
        return dict([(k, v.genome_container) for k, v in self._stores.items()])

    @genome_container.setter
    def genome_container(self, x: Union[GenomeContainer, Dict[str, GenomeContainer]]):
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
    def obs_by_key(self) -> Dict[str, Union[pd.DataFrame, dask.dataframe.DataFrame]]:
        """
        Dictionary of all anndata instances for each selected data set in store, sub-setted by selected cells, for each
        stores.
        """
        return dict([(k, v.obs) for k, v in self.adata_by_key.items()])

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
    def n_obs(self) -> int:
        """
        Dictionary of number of observations across stores.
        """
        return np.asarray(np.sum([v.n_obs for v in self.stores.values()]), dtype="int32")

    @property
    def n_obs_dict(self) -> Dict[str, int]:
        """
        Dictionary of number of observations by store.
        """
        return dict([(k, v.n_obs) for k, v in self.stores.items()])

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
            - "cell_type" points to self.cell_type_obs_key
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
        keys_found = []
        for k, v in self.stores.items():
            indices_k = {}
            for kk, vv in indices.items():
                if kk in v.adata_by_key.keys():
                    indices_k[kk] = vv
                    keys_found.append(kk)
            self.stores[k].indices = indices_k
        # Make sure all declared data were assigned to stores:
        keys_not_found = list(set(list(indices.keys())).difference(set(keys_found)))
        if len(keys_not_found) > 0:
            raise ValueError(f"did not find object(s) with name(s) in store: {keys_not_found}")

    def checkout(
            self,
            idx: Union[Dict[str, Union[np.ndarray, None]], None] = None,
            intercalated: bool = True,
            **kwargs
    ) -> CartMulti:
        """
        Carts per store.

        See also DistributedStore*.checkout().

        :param idx:
        :param intercalated: Whether to do sequential or intercalated emission.
        :param kwargs: See parameters of DistributedStore*.generator().
        :return: Generator function which yields batch_size at every invocation.
            The generator returns a tuple of (.X, .obs).
        """
        if idx is None:
            idx = dict([(k, None) for k in self.stores.keys()])
        for k in self.stores.keys():
            assert k in idx.keys(), (idx.keys(), self.stores.keys())
        carts = dict([(k, v.checkout(idx=idx[k], **kwargs)) for k, v in self.stores.items()])
        return CartMulti(carts=carts, intercalated=intercalated)


class StoresAnndata(StoreMultipleFeatureSpaceBase):

    def __init__(self, adatas: Union[anndata.AnnData, List[anndata.AnnData], Tuple[anndata.AnnData]]):
        # Collect all data loaders from files in directory:
        self._adata_ids_sfaira = AdataIdsSfaira()
        adata_by_key = {}
        indices = {}
        if isinstance(adatas, anndata.AnnData):
            adatas = [adatas]
        for i, adata in enumerate(adatas):
            # Check if adata has a unique ID, if not, add one:
            if self._adata_ids_sfaira.id not in adata.uns.keys():
                adata.uns[self._adata_ids_sfaira.id] = f"adata_{i}"
            if self._adata_ids_sfaira.organism in adata.uns.keys():
                organism = adata.uns[self._adata_ids_sfaira.organism]
            else:
                # Declare as unknown organism and genome and make a group of its own:
                organism = adata.uns[self._adata_ids_sfaira.id]
            if isinstance(organism, list):
                if len(organism) == 1:
                    organism = organism[0]
                    assert isinstance(organism, str), organism
                else:
                    raise ValueError(f"tried to register mixed organism data set ({organism})")
            adata_id = adata.uns[self._adata_ids_sfaira.id]
            # Make up a new merged ID for data set indexing if there is a list of IDs in .uns.
            if isinstance(adata_id, list):
                adata_id = "_".join(adata_id)
            if organism not in adata_by_key.keys():
                adata_by_key[organism] = {}
                indices[organism] = {}
            try:
                adata_by_key[organism][adata_id] = adata
                indices[organism][adata_id] = np.arange(0, adata.n_obs)
            except TypeError as e:
                raise TypeError(f"{e} for {organism} or {adata.uns[self._adata_ids_sfaira.id]}")
        stores = dict([
            (k, StoreAnndata(adata_by_key=adata_by_key[k], indices=indices[k], in_memory=True))
            for k in adata_by_key.keys()
        ])
        super(StoresAnndata, self).__init__(stores=stores)


class StoresDao(StoreMultipleFeatureSpaceBase):

    _dataset_weights: Union[None, Dict[str, float]]

    def __init__(self,
                 cache_path: Union[str, os.PathLike, List[str], List[os.PathLike]],
                 columns: Union[None, List[str]] = None):
        """

        :param cache_path: Store directory.
        :param columns: Which columns to read into the obs copy in the output, see pandas.read_parquet().
        """
        # Collect all data loaders from files in directory:
        self._adata_ids_sfaira = AdataIdsSfaira()
        adata_by_key = {}
        x_by_key = {}
        indices = {}
        if not isinstance(cache_path, list) or isinstance(cache_path, tuple) or isinstance(cache_path, np.ndarray):
            cache_path = [cache_path]
        for cache_path_i in cache_path:
            for f in np.sort(os.listdir(cache_path_i)):
                adata = None
                x = None
                trial_path = os.path.join(cache_path_i, f)
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
                    if adata.uns[self._adata_ids_sfaira.id] in adata_by_key[organism].keys():
                        print(f"WARNING: overwriting store entry in {adata.uns[self._adata_ids_sfaira.id]} in store "
                              f"{cache_path_i}.")
                    adata_by_key[organism][adata.uns[self._adata_ids_sfaira.id]] = adata
                    x_by_key[organism][adata.uns[self._adata_ids_sfaira.id]] = x
                    indices[organism][adata.uns[self._adata_ids_sfaira.id]] = np.arange(0, adata.n_obs)
        stores = dict([
            (k, StoreDao(adata_by_key=adata_by_key[k], x_by_key=x_by_key[k], indices=indices[k],
                         obs_by_key=None))
            for k in adata_by_key.keys()
        ])
        self._x_by_key = x_by_key
        super(StoresDao, self).__init__(stores=stores)


class StoresH5ad(StoreMultipleFeatureSpaceBase):

    def __init__(
            self,
            cache_path: Union[str, os.PathLike, List[str], List[os.PathLike]],
            in_memory: bool = False):
        # Collect all data loaders from files in directory:
        self._adata_ids_sfaira = AdataIdsSfaira()
        adata_by_key = {}
        indices = {}
        if not isinstance(cache_path, list) or isinstance(cache_path, tuple) or isinstance(cache_path, np.ndarray):
            cache_path = [cache_path]
        for cache_path_i in cache_path:
            for f in np.sort(os.listdir(cache_path_i)):
                adata = None
                trial_path = os.path.join(cache_path_i, f)
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
                    if adata.uns[self._adata_ids_sfaira.id] in adata_by_key[organism].keys():
                        print(f"WARNING: overwriting store entry in {adata.uns[self._adata_ids_sfaira.id]} in store "
                              f"{cache_path_i}.")
                    adata_by_key[organism][adata.uns[self._adata_ids_sfaira.id]] = adata
                    indices[organism][adata.uns[self._adata_ids_sfaira.id]] = np.arange(0, adata.n_obs)
        stores = dict([
            (k, StoreAnndata(adata_by_key=adata_by_key[k], indices=indices[k], in_memory=in_memory))
            for k in adata_by_key.keys()
        ])
        super(StoresH5ad, self).__init__(stores=stores)
