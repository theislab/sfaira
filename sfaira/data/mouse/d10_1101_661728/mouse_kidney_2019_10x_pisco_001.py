import anndata
import numpy as np
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
        DatasetTms.__init__(self=self, path=path, meta_path=meta_path, **kwargs)

        self.id = "mouse_kidney_2019_10x_pisco_001_10.1101/661728"
        self.source = source
        if self.source == "aws":
            self.download = "https://czb-tabula-muris-senis.s3-us-west-2.amazonaws.com/Data-objects/"
        elif self.source == "figshare":
            self.download = "https://ndownloader.figshare.com/articles/8273102/versions/2"
        else:
            raise ValueError("source %s not recognized" % self.source)
        self.organ = "kidney"
        self.sub_tissue = "kidney"

        self.class_maps = {
            "0": {
                "kidney capillary endothelial cell": "endothelial cell",
                "kidney mesangial cell": "mesangial cell",
                "kidney interstitial fibroblast": "interstitial fibroblast",
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            if self.source == "aws":
                fn = os.path.join(self.path, "mouse", "kidney", "tabula-muris-senis-droplet-processed-official-annotations-Kidney.h5ad")
            elif self.source == "figshare":
                fn = os.path.join(self.path, "mouse", "kidney", "Kidney_droplet.h5ad")
            else:
                raise ValueError("source %s not recognized" % self.source)
        self._load_tms(fn=fn)
        self.set_unkown_class_id(ids=["nan", "kidney cell"])
        