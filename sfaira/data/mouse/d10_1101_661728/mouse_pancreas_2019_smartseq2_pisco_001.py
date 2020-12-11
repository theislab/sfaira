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
        DatasetTms.__init__(self=self, path=path, meta_path=meta_path, **kwargs)

        self.id = "mouse_pancreas_2019_smartseq2_pisco_001_10.1101/661728"
        self.source = source
        if self.source == "aws":
            self.download = "https://czb-tabula-muris-senis.s3-us-west-2.amazonaws.com/Data-objects/"
        elif self.source == "figshare":
            self.download = "https://ndownloader.figshare.com/articles/8273102/versions/2"
        else:
            raise ValueError("source %s not recognized" % self.source)
        self.organ = "pancreas"
        self.sub_tissue = "pancreas"

        self.class_maps = {
            "0": {
                "pancreatic ductal cel": "pancreatic ductal cell"
            },
        }

    def _load(self, fn=None):
        if fn is None:
            if self.path is None:
                raise ValueError("provide either fn in load or path in constructor")
            if self.source == "aws":
                fn = os.path.join(self.path, "mouse", "pancreas", "tabula-muris-senis-facs-processed-official-annotations-Pancreas.h5ad")
            elif self.source == "figshare":
                fn = os.path.join(self.path, "mouse", "pancreas", "Pancreas_facs.h5ad")
            else:
                raise ValueError("source %s not recognized" % self.source)
        self._load_tms(fn=fn)
