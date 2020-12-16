import anndata
import os
from typing import Union
from .external import DatasetTms


class Dataset(DatasetTms):

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            source: str = "aws",
            **kwargs
    ):
        super().__init__(path=path, meta_path=meta_path, source=source, **kwargs)
        self.id = "mouse_adipose_2019_smartseq2_pisco_001_10.1101/661728"
        self.organ = "adipose"
        self.sub_tissue = "adipose"        
        self.protocol = self._get_protocol_tms(self.id)        

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            if self.source == "aws":
                fn = os.path.join(self.path, "mouse", "adipose", "tabula-muris-senis-facs-processed-official-annotations-BAT.h5ad")
            elif self.source == "figshare":
                fn = os.path.join(self.path, "mouse", "adipose", "bat_facs.h5ad")
            else:
                raise ValueError("source %s not recognized" % self.source)
        self._load_tms(fn=fn)
