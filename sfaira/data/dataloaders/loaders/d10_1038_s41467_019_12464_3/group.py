from typing import Union

from sfaira.data import DatasetGroupLoadingManyFiles

from .human_mixed_2019_10x_szabo_001 import Dataset


class Group(DatasetGroupLoadingManyFiles):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None
    ):
        sample_fns = [
            "GSM3589406_PP001swap.filtered.matrix.txt.gz",
            "GSM3589407_PP002swap.filtered.matrix.txt.gz",
            "GSM3589408_PP003swap.filtered.matrix.txt.gz",
            "GSM3589409_PP004swap.filtered.matrix.txt.gz",
            "GSM3589410_PP005swap.filtered.matrix.txt.gz",
            "GSM3589411_PP006swap.filtered.matrix.txt.gz",
            "GSM3589412_PP009swap.filtered.matrix.txt.gz",
            "GSM3589413_PP010swap.filtered.matrix.txt.gz",
            "GSM3589414_PP011swap.filtered.matrix.txt.gz",
            "GSM3589415_PP012swap.filtered.matrix.txt.gz",
            "GSM3589416_PP013swap.filtered.matrix.txt.gz",
            "GSM3589417_PP014swap.filtered.matrix.txt.gz",
            "GSM3589418_PP017swap.filtered.matrix.txt.gz",
            "GSM3589419_PP018swap.filtered.matrix.txt.gz",
            "GSM3589420_PP019swap.filtered.matrix.txt.gz",
            "GSM3589421_PP020swap.filtered.matrix.txt.gz",
        ]
        super().__init__(
            cls=Dataset,
            sample_fns=sample_fns,
            path=path, meta_path=meta_path, cache_path=cache_path
        )
