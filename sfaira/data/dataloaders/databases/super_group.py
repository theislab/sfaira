from typing import Union

from sfaira.data.dataloaders.base.dataset_group import DatasetGroup, DatasetSuperGroup
from sfaira.data.dataloaders.databases.cellxgene import DatasetSuperGroupCellxgene


class DatasetSuperGroupDatabases(DatasetSuperGroup):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            cache_metadata: bool = False,
    ):
        dataset_super_groups = [
            DatasetSuperGroupCellxgene(
                data_path=data_path,
                meta_path=meta_path,
                cache_path=cache_path,
                cache_metadata=cache_metadata,
            ),
        ]
        super().__init__(dataset_groups=dataset_super_groups)
