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
        self.id = "mouse_bone_2019_10x_pisco_001_10.1101/661728"
        self.organ = "bone"
        self.sub_tissue = "marrow"        
        self.protocol = self._get_protocol_tms(self.id)        

        self.class_maps = {
            "0": {},
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            if self.source == "aws":
                fn = os.path.join(self.path, "mouse", "bone", "tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad")
            elif self.source == "figshare":
                fn = os.path.join(self.path, "mouse", "bone", "Marrow_droplet.h5ad")
            else:
                raise ValueError("source %s not recognized" % self.source)
        self._load_tms(fn=fn)
        self.set_unkown_class_id(ids=["unknown cell-1", "unknown cell-2"])
