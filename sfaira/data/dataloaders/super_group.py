from typing import Union

from sfaira.data.dataloaders.loaders import DatasetSuperGroupLoaders
from sfaira.data.dataloaders.databases import DatasetSuperGroupDatabases
from sfaira.data import DatasetSuperGroup


class DatasetSuperGroupSfaira(DatasetSuperGroup):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
    ):
        """
        Nested super group of data loaders, unifying data set wise data loader SuperGroup and the database
        interface SuperGroup.

        :param path:
        :param meta_path:
        :param cache_path:
        """
        super().__init__(dataset_groups=[
            DatasetSuperGroupLoaders(
                path=path,
                meta_path=meta_path,
                cache_path=cache_path,
            ),
            DatasetSuperGroupDatabases(
                path=path,
                meta_path=meta_path,
                cache_path=cache_path,
            )
        ])
