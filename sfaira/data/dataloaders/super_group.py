from typing import Union

try:
    import sfaira_extension as sfairae
except ImportError:
    sfairae = None

from sfaira.data.dataloaders.loaders import DatasetSuperGroupLoaders
from sfaira.data.dataloaders.databases import DatasetSuperGroupDatabases
from sfaira.data.dataloaders.base.dataset_group import DatasetSuperGroup


class Universe(DatasetSuperGroup):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            exclude_databases: bool = True,
    ):
        """
        Nested super group of data loaders, unifying data set wise data loader SuperGroup and the database
        interface SuperGroup.

        :param data_path:
        :param meta_path:
        :param cache_path:
        """
        dsgs = [
            DatasetSuperGroupLoaders(
                data_path=data_path,
                meta_path=meta_path,
                cache_path=cache_path,
            ),
        ]
        if not exclude_databases:
            dsgs.append(
                DatasetSuperGroupDatabases(
                    data_path=data_path,
                    meta_path=meta_path,
                    cache_path=cache_path,
                )
            )
        if sfairae is not None:
            dsgs.append(sfairae.data.dataloaders.loaders.DatasetSuperGroupLoaders(
                data_path=data_path,
                meta_path=meta_path,
                cache_path=cache_path,
            ))
        super().__init__(dataset_groups=dsgs)
