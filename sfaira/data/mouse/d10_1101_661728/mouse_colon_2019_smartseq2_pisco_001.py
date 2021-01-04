import anndata
import os
from typing import Union
from .external import DatasetTms


class Dataset(DatasetTms):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            source: str = "aws",
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, cache_path=cache_path, source=source, **kwargs)
        self.id = "mouse_colon_2019_smartseq2_pisco_001_10.1101/661728"
        self.organ = "colon"
        self.sub_tissue = "colon"

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.path, "mouse", "colon", "Large_Intestine_facs.h5ad")
            if self.source == "aws":
                fn = os.path.join(self.path, "mouse", "colon", "tabula-muris-senis-facs-processed-official-annotations-Large_Intestine.h5ad")
            elif self.source == "figshare":
                fn = os.path.join(self.path, "mouse", "colon", "Large_Intestine_facs.h5ad")
            else:
                raise ValueError("source %s not recognized" % self.source)
        self._load_tms(fn=fn)
