import logging

import os
from rich import print
from typing import Dict

from sfaira.commands.utils import get_ds
from sfaira.consts.utils import clean_doi
from sfaira.data import DatasetBase

log = logging.getLogger(__name__)


class H5adExport:

    datasets: Dict[str, DatasetBase]
    doi: str
    doi_sfaira_repr: str
    path_cache: str
    path_data: str
    path_loader: str
    path_out: str
    schema: str

    def __init__(self, doi, path_cache, path_data, path_loader, path_out, schema):
        self.doi = doi
        self.doi_sfaira_repr = clean_doi(self.doi)
        self.path_cache = path_cache
        self.path_data = path_data
        self.path_loader = path_loader
        self.path_out = path_out
        if schema not in ["cellxgene"]:
            raise ValueError(f"Did not recognize schema {schema}")
        self.schema = schema

    def write(self):
        self._load_objects()
        self._write_h5ads()

    def _load_objects(self):
        dsg, _ = get_ds(doi_sfaira_repr=self.doi_sfaira_repr, path_cache=self.path_cache, path_data=self.path_data,
                        path_loader=self.path_loader)
        dsg.load(load_raw=False, allow_caching=True)
        if self.schema == "cellxgene":
            dsg.streamline_features(schema="cellxgene:" + "2.0.0")
            dsg.streamline_metadata(
                schema=self.schema.lower(),
                clean_obs=False,
                clean_var=False,
                clean_uns=True,
                clean_obs_names=False,
                keep_orginal_obs=False,
                keep_symbol_obs=True,
                keep_id_obs=True,
            )
            dsg.collapse_counts()
        self.datasets = dsg.datasets

    def _write_h5ads(self):
        counter = 0
        for k, v in self.datasets.items():
            fn = v.doi_cleaned_id + ".h5ad"
            dir_name = v.directory_formatted_doi
            if not os.path.exists(os.path.join(self.path_out, dir_name)):
                os.makedirs(os.path.join(self.path_out, dir_name))
            fn_out = os.path.join(self.path_out, dir_name, fn)
            print(f'[bold orange]Sfaira butler: "Preparing {fn_out} for you."')
            v.adata.write_h5ad(fn_out)
            counter += 1
        print(f'[bold orange]Sfaira butler: "I wrote a total of {counter} .h5ad files."')
