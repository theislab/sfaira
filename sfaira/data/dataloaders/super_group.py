from typing import Union

from sfaira.data import DatasetSuperGroupSfairaBase
from sfaira.data.dataloaders.databases import DatasetSuperGroupDatabases


class DatasetSuperGroupSfaira(DatasetSuperGroupSfairaBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
    ):
        # Load all organisms:
        super().__init__(
            file_base=__file__,
            dir_prefix="",
            dir_exlcude=["databases"],
            path=path,
            meta_path=meta_path,
            cache_path=cache_path,
        )
        # Load structured data bases
        self.extend_dataset_groups([DatasetSuperGroupDatabases(
            path=path,
            meta_path=meta_path,
            cache_path=cache_path,
        )])
