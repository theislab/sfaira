import anndata
from IPython.display import display_javascript, display_html
import json
import os
import pathlib
import pickle
import requests
from typing import List, Union
import uuid

from sfaira.data.dataloaders.base import DatasetBase
from sfaira.consts import AdataIdsCellxgene, AdataIdsCellxgeneHuman_v1_1_0, AdataIdsCellxgeneMouse_v1_1_0
from sfaira.consts.directories import CACHE_DIR_DATABASES_CELLXGENE
from sfaira.data.dataloaders.databases.cellxgene.rest_helpers import get_collection, get_data
from sfaira.data.dataloaders.databases.cellxgene.rest_helpers import CELLXGENE_PRODUCTION_ENDPOINT, DOWNLOAD_DATASET


def cellxgene_fn(dir, dataset_id):
    return os.path.join(dir, dataset_id + ".h5ad")


def clean_cellxgene_meta_obs(k, val, adata_ids) -> Union[str, List[str]]:
    """
    :param k: Found meta data name.
    :param val: Found meta data entry.
    :returns: Cleaned meta data entry.
    """
    if k == "disease":
        # TODO normal state label varies in disease annotation. This can be removed once streamlined.
        val = ["healthy" if (v.lower() == "normal" or v.lower() == "healthy") else v for v in val]
    elif k == "organ":
        # Organ labels contain labels on tissue type also, such as 'UBERON:0001911 (cell culture)'.
        val = [v.split(" ")[0] for v in val]
    elif k == "organism":
        organism_map = {
            "Homo sapiens": "human",
            "Mus musculus": "mouse",
        }
        val = [organism_map[v] if v in organism_map.keys() else v for v in val]
    return val


def clean_cellxgene_meta_uns(k, val, adata_ids) -> Union[str, List[str]]:
    """
    :param k: Found meta data name.
    :param val: Found meta data entry.
    :returns: Cleaned meta data entry.
    """
    x_clean = []
    for v in val:
        if k == "sex":
            v = v[0]
        else:
            # Decide if labels are read from name or ontology ID:
            if k == "disease" and (v["label"].lower() == "normal" or v["label"].lower() == "healthy"):
                # TODO normal state label varies in disease annotation. This can be removed once streamlined.
                v = "healthy"
            elif k in ["assay_sc", "disease", "organ"] and \
                    v["ontology_term_id"] != adata_ids.unknown_metadata_identifier:
                v = v["ontology_term_id"]
            else:
                v = v["label"]
            # Organ labels contain labels on tissue type also, such as 'UBERON:0001911 (cell culture)'.
            if k == "organ":
                v = v.split(" ")[0]
            if k == "organism":
                organism_map = {
                    "Homo sapiens": "human",
                    "Mus musculus": "mouse",
                }
                if v not in organism_map:
                    raise ValueError(f"value {v} not recognized")
                v = organism_map[v]
        if v != adata_ids.unknown_metadata_identifier and v != adata_ids.invalid_metadata_identifier:
            x_clean.append(v)
    return x_clean


class Dataset(DatasetBase):
    """
    This is a dataloader for downloaded h5ad from cellxgene.

    In contrast to the base class, each instance is coupled to a particular collection_id to allow query.
    In the base classes, collection_id are only defined on the group level.
    """

    collection_id: str

    def __init__(
            self,
            collection_id: str = "default",
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            load_func=None,
            dict_load_func_annotation=None,
            yaml_path: Union[str, None] = None,
            sample_fn: Union[str, None] = None,
            sample_fns: Union[List[str], None] = None,
            additional_annotation_key: Union[str, None] = None,
            verbose: int = 0,
            **kwargs
    ):
        super().__init__(
            data_path=data_path,
            meta_path=meta_path,
            cache_path=cache_path,
            load_func=load_func,
            sample_fn=sample_fn,
            sample_fns=sample_fns,
        )
        # General keys are defined in the shared IDs object. Further down, the species specific one is loaded to
        # disambiguate species-dependent differences.
        self._adata_ids_cellxgene = AdataIdsCellxgene()
        self._collection = None
        self.collection_id = collection_id
        self.supplier = "cellxgene"
        doi = [x['link_url'] for x in self.collection["links"] if x['link_type'] == 'DOI']
        self.doi_journal = collection_id if len(doi) == 0 else doi[0]  # TODO access journal DOI explicitly.
        self.id = sample_fn
        # Set h5ad download URLs:
        download_url_data = []
        for asset in self._collection_dataset['dataset_assets']:
            if asset['filetype'].lower() == "h5ad":
                download_url_data.append(CELLXGENE_PRODUCTION_ENDPOINT + DOWNLOAD_DATASET + asset['dataset_id'] +
                                         "/asset/" + asset['id'])
        self.download_url_data = download_url_data

        # Set dataset-wise attributes based on object preview from REST API (without h5ad download):
        # Set organism first so that other terms can access this attribute (e.g. developmental_stage ontology).
        reordered_keys = ["organism"] + [x for x in self._adata_ids_cellxgene.dataset_keys if x != "organism"]
        for k in reordered_keys:
            val = self._collection_dataset[getattr(self._adata_ids_cellxgene, k)]
            # Unique label if list is length 1:
            # Otherwise do not set property and resort to cell-wise labels.
            if isinstance(val, dict) or k == "sex":
                val = [val]
            v_clean = clean_cellxgene_meta_uns(k=k, val=val, adata_ids=self._adata_ids_cellxgene)
            try:
                # Set as single element or list if multiple entries are given.
                if len(v_clean) == 1:
                    setattr(self, k, v_clean[0])
                else:
                    setattr(self, k, v_clean)
            except ValueError as e:
                if verbose > 0:
                    print(f"WARNING: {e} in {self.collection_id} and data set {self.id}")

        if self.organism == "human":
            self._adata_ids_cellxgene = AdataIdsCellxgeneHuman_v1_1_0()
        elif self.organism == "mouse":
            self._adata_ids_cellxgene = AdataIdsCellxgeneMouse_v1_1_0()
        else:
            assert False, self.organism
        # Add author information.  # TODO need to change this to contributor?
        setattr(self, "author", "cellxgene")
        # The h5ad objects from cellxgene follow a particular structure and the following attributes are guaranteed to
        # be in place. Note that these point at the anndata instance and will only be available for evaluation after
        # download. See below for attributes that are lazily available
        self.cell_type_obs_key = self._adata_ids_cellxgene.cell_type
        self.development_stage_obs_key = self._adata_ids_cellxgene.development_stage
        self.disease_obs_key = self._adata_ids_cellxgene.disease
        self.ethnicity_obs_key = self._adata_ids_cellxgene.ethnicity
        self.sex_obs_key = self._adata_ids_cellxgene.sex
        self.organ_obs_key = self._adata_ids_cellxgene.organism
        self.state_exact_obs_key = self._adata_ids_cellxgene.state_exact

        self.gene_id_symbols_var_key = self._adata_ids_cellxgene.feature_symbol

        self._unknown_celltype_identifiers = self._adata_ids_cellxgene.unknown_metadata_identifier

    @property
    def _collection_cache_dir(self):
        """
        The cache dir is in a cache directory in the sfaira installation that is excempt from git versioning.
        """
        cache_dir_path = pathlib.Path(CACHE_DIR_DATABASES_CELLXGENE)
        cache_dir_path.mkdir(parents=True, exist_ok=True)
        return CACHE_DIR_DATABASES_CELLXGENE

    @property
    def _collection_cache_fn(self):
        return os.path.join(self._collection_cache_dir, self.collection_id + ".pickle")

    @property
    def collection(self):
        if self._collection is None:
            # Check if cached:
            if os.path.exists(self._collection_cache_fn):
                with open(self._collection_cache_fn, "rb") as f:
                    self._collection = pickle.load(f)
            else:
                # Download and cache:
                self._collection = get_collection(collection_id=self.collection_id)
                with open(self._collection_cache_fn, "wb") as f:
                    pickle.dump(obj=self._collection, file=f)
        return self._collection

    @property
    def _collection_dataset(self):
        return self.collection['datasets'][self._sample_fns.index(self.sample_fn)]

    @property
    def directory_formatted_doi(self) -> str:
        return self.collection_id

    @property
    def doi_cleaned_id(self):
        return self.id

    def load(
            self,
            remove_gene_version: bool = True,
            match_to_reference: Union[str, bool, None] = None,
            load_raw: bool = False,
            allow_caching: bool = True,
            set_metadata: bool = True,
            **kwargs
    ):
        # Invoke load with cellxgene adapted parameters:
        #  - Never cache as the cellxgene objects already fast to read.
        super().load(
            remove_gene_version=False,
            match_to_reference=match_to_reference,
            load_raw=True,
            allow_caching=False,
            set_metadata=set_metadata,
            adata_ids=self._adata_ids_cellxgene,
            **kwargs
        )

    def download(self, filetype: str = "h5ad", verbose: int = 0):
        """

        Only download if file does not already exist.

        :param filetype: File type to download.

            - "h5ad"
            - "rds"
            - "loom"
        """
        counter = 0
        if not os.path.exists(os.path.join(self.data_dir_base, self.directory_formatted_doi)):
            os.makedirs(os.path.join(self.data_dir_base, self.directory_formatted_doi))
        for asset in self._collection_dataset['dataset_assets']:
            if asset['filetype'].lower() == filetype:
                # Only download if file does not already exist:
                fn = cellxgene_fn(dir=self.data_dir, dataset_id=self.sample_fn)
                if not os.path.isfile(fn):
                    counter += 1
                    assert counter < 2, f"found more than one {filetype} for data set {self.sample_fn}"
                    url = CELLXGENE_PRODUCTION_ENDPOINT + DOWNLOAD_DATASET + asset['dataset_id'] + "/asset/" + \
                        asset['id']
                    r = requests.post(url)
                    r.raise_for_status()
                    presigned_url = r.json()['presigned_url']
                    # Report:
                    headers = {'range': 'bytes=0-0'}
                    r1 = requests.get(presigned_url, headers=headers)
                    if r1.status_code == requests.codes.partial:
                        if verbose > 0:
                            print(f"Downloading {r1.headers['Content-Range']} from {r1.headers['Server']}")
                    # Download:
                    open(fn, 'wb').write(get_data(presigned_url=presigned_url))

    def show_summary(self):
        uuid_session = str(uuid.uuid4())
        display_html('<div id="{}" style="height: 600px; width:100%;"></div>'.format(uuid_session), raw=True)
        display_javascript("""
        require(["https://rawgit.com/caldwell/renderjson/master/renderjson.js"], function() {
          document.getElementById('%s').appendChild(renderjson(%s))
        });
        """ % (uuid_session, json.dumps(self._collection_dataset)), raw=True)


def load(data_dir, sample_fn, adata_ids: AdataIdsCellxgene, **kwargs):
    """
    Generalised load function for cellxgene-provided data sets.

    This function corresponds to the dataset-wise load() functions defined in standard sfaira data loaders.
    """
    fn = cellxgene_fn(dir=data_dir, dataset_id=sample_fn)
    adata = anndata.read_h5ad(fn)
    if adata.raw is not None:  # TODO still need this?
        adata.X = adata.raw.X
        del adata.raw
    for k in adata_ids.ontology_constrained:
        col_name = getattr(adata_ids, k)
        if col_name in adata.obs.columns:
            adata.obs[col_name] = clean_cellxgene_meta_obs(k=k, val=adata.obs[col_name].values, adata_ids=adata_ids)
    return adata
