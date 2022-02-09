import os
from typing import List, Union

from sfaira.data.store.multi_store import DistributedStoresDao, DistributedStoresH5ad, \
    DistributedStoreMultipleFeatureSpaceBase


def load_store(cache_path: Union[str, os.PathLike], store_format: str = "dao",
               columns: Union[None, List[str]] = None) -> DistributedStoreMultipleFeatureSpaceBase:
    """
    Instantiates a distributed store class.

    Note that any store is instantiated as a DistributedStoreMultipleFeatureSpaceBase.
    This instances can be subsetted to the desired single feature space.

    :param cache_path: Store directory.
    :param store_format: Format of store {"h5ad", "dao"}.

        - "h5ad": Returns instance of DistributedStoreH5ad and keeps data in memory. See also "h5ad_backed".
        - "dao": Returns instance of DistributedStoreDoa (distributed access optimized).
        - "h5ad_backed": Returns instance of DistributedStoreH5ad and keeps data as backed (out of memory). See also
            "h5ad".
    :param columns: Which columns to read into the obs copy in the output, see pandas.read_parquet().
        Only relevant if store_format is "dao".
    :return: Instances of a distributed store class.
    """
    if store_format == "h5ad":
        return DistributedStoresH5ad(cache_path=cache_path, in_memory=True)
    elif store_format == "dao":
        return DistributedStoresDao(cache_path=cache_path, columns=columns)
    elif store_format == "h5ad_backed":
        return DistributedStoresH5ad(cache_path=cache_path, in_memory=False)
    else:
        raise ValueError(f"Did not recognize store_format {store_format}.")
