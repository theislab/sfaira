import abc
import anndata
import dask.dataframe
import numpy as np
import os
import pandas as pd
import sys
from typing import Dict, List, Union


class DistributedStoreBase(abc.ABC):
    """
    Base class for store API for attribute typing and shared methods.
    """

    @property
    @abc.abstractmethod
    def adata_by_key(self) -> Dict[str, anndata.AnnData]:
        pass

    @property
    @abc.abstractmethod
    def data_by_key(self):
        pass

    @property
    @abc.abstractmethod
    def indices(self) -> Dict[str, np.ndarray]:
        pass

    @property
    @abc.abstractmethod
    def genome_container(self):
        pass

    @property
    @abc.abstractmethod
    def n_obs(self) -> int:
        pass

    @property
    @abc.abstractmethod
    def n_vars(self):
        pass

    @property
    @abc.abstractmethod
    def obs(self):
        pass

    @property
    @abc.abstractmethod
    def obs_by_key(self) -> Dict[str, Union[pd.DataFrame, dask.dataframe.DataFrame]]:
        pass

    @property
    @abc.abstractmethod
    def var_names(self):
        pass

    @property
    @abc.abstractmethod
    def shape(self):
        pass

    @property
    @abc.abstractmethod
    def var(self):
        pass

    @property
    @abc.abstractmethod
    def X(self):
        pass

    @abc.abstractmethod
    def subset(self, attr_key, values: Union[str, List[str], None],
               excluded_values: Union[str, List[str], None], verbose: int):
        pass

    @abc.abstractmethod
    def write_config(self, fn: Union[str, os.PathLike]):
        pass

    @abc.abstractmethod
    def load_config(self, fn: Union[str, os.PathLike]):
        pass

    @abc.abstractmethod
    def generator(
            self,
            idx: Union[np.ndarray, None],
            batch_size: int,
            obs_keys: List[str],
            return_dense: bool,
            randomized_batch_access: bool,
            random_access: bool,
            **kwargs
    ) -> iter:
        pass

    @property
    def adata_memory_footprint(self) -> Dict[str, float]:
        """
        Memory foot-print of data set k in MB.
        """
        return dict([(k, sys.getsizeof(v) / np.power(1024, 2)) for k, v in self.adata_by_key.items()])
