from typing import Union

from sfaira.data import DatasetSuperGroupOrganismBase


class DatasetSuperGroupDirectoryOriented(DatasetSuperGroupOrganismBase):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
    ):
        super().__init__(
            file_base=__file__,
            path=path,
            meta_path=meta_path,
            cache_path=cache_path
        )
