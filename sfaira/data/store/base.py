import abc
import anndata
import dask.dataframe
import numpy as np
import os
import pandas as pd
from typing import Dict, List, Union


class DistributedStoreBase(abc.ABC):
    """
    Base class for store API for attribute typing.
    """

    @abc.abstractmethod
    def adata_by_key(self) -> Dict[str, anndata.AnnData]:
        pass

    @abc.abstractmethod
    def data_by_key(self):
        pass

    @abc.abstractmethod
    def indices(self) -> Dict[str, np.ndarray]:
        pass

    @abc.abstractmethod
    def obs_by_key(self) -> Dict[str, Union[pd.DataFrame, dask.dataframe.DataFrame]]:
        pass

    @abc.abstractmethod
    def genome_container(self):
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
    def var_names(self):
        pass

    @abc.abstractmethod
    def n_vars(self):
        pass

    @abc.abstractmethod
    def n_obs(self) -> int:
        pass

    @abc.abstractmethod
    def shape(self):
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

    @abc.abstractmethod
    def X(self):
        pass

    @abc.abstractmethod
    def obs(self):
        pass
