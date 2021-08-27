from __future__ import annotations

import abc
import anndata
from anndata.utils import make_index_unique
import h5py
import numpy as np
import pandas as pd
import os
from os import PathLike
import pandas
import scipy.sparse
from typing import Dict, List, Tuple, Union
import warnings
import urllib.request
import urllib.parse
import urllib.error
import cgi
import ssl

from sfaira.versions.genomes import GenomeContainer
from sfaira.versions.metadata import Ontology, OntologyHierarchical, CelltypeUniverse
from sfaira.consts import AdataIds, AdataIdsCellxgeneGeneral, AdataIdsCellxgeneHuman_v1_1_0, AdataIdsCellxgeneMouse_v1_1_0, \
    AdataIdsSfaira, META_DATA_FIELDS, OCS
from sfaira.data.dataloaders.export_adaptors import cellxgene_export_adaptor
from sfaira.data.store.io_dao import write_dao
from sfaira.data.dataloaders.base.utils import is_child, get_directory_formatted_doi
from sfaira.data.utils import collapse_matrix, read_yaml
from sfaira.consts.utils import clean_id_str


load_doc = \
    """
    :param remove_gene_version: Remove gene version string from ENSEMBL ID so that different versions in different data sets are superimposed.
    :param match_to_reference: Reference genomes name or False to keep original feature space.
    :param load_raw: Loads unprocessed version of data if available in data loader.
    :param allow_caching: Whether to allow method to cache adata object for faster re-loading.
    """


class DatasetBase(abc.ABC):
    adata: Union[None, anndata.AnnData]
    class_maps: dict
    _meta: Union[None, pandas.DataFrame]
    data_dir_base: Union[None, str]
    meta_path: Union[None, str]
    cache_path: Union[None, str]
    id: Union[None, str]
    genome: Union[None, str]
    supplier: str

    _assay_sc: Union[None, str]
    _assay_differentiation: Union[None, str]
    _assay_type_differentiation: Union[None, str]
    _author: Union[None, str]
    _bio_sample: Union[None, str]
    _cell_line: Union[None, str]
    _cell_type: Union[None, str]
    _default_embedding: Union[None, str]
    _development_stage: Union[None, str]
    _disease: Union[None, str]
    _doi_journal: Union[None, str]
    _doi_preprint: Union[None, str]
    _download_url_data: Union[Tuple[List[None]], Tuple[List[str]], None]
    _download_url_meta: Union[Tuple[List[None]], Tuple[List[str]], None]
    _ethnicity: Union[None, str]
    _id: Union[None, str]
    _individual: Union[None, str]
    _ncells: Union[None, int]
    _normalization: Union[None, str]
    _organ: Union[None, str]
    _organism: Union[None, str]
    _primary_data: Union[None, bool]
    _sex: Union[None, str]
    _source: Union[None, str]
    _sample_source: Union[None, str]
    _state_exact: Union[None, str]
    _title: Union[None, str]
    _bio_sample: Union[None, str]
    _year: Union[None, int]

    assay_sc_obs_key: Union[None, str]
    assay_differentiation_obs_key: Union[None, str]
    assay_type_differentiation_obs_key: Union[None, str]
    assay_cell_line_obs_key: Union[None, str]
    bio_sample_obs_key: Union[None, str]
    cell_type_obs_key: Union[None, str]
    development_stage_obs_key: Union[None, str]
    disease_obs_key: Union[None, str]
    ethnicity_obs_key: Union[None, str]
    individual_obs_key: Union[None, str]
    organ_obs_key: Union[None, str]
    organism_obs_key: Union[None, str]
    sample_source_obs_key: Union[None, str]
    sex_obs_key: Union[None, str]
    state_exact_obs_key: Union[None, str]
    tech_sample_obs_key: Union[None, str]

    gene_id_symbols_var_key: Union[None, str]
    gene_id_ensembl_var_key: Union[None, str]

    _celltype_universe: Union[None, CelltypeUniverse]
    _ontology_class_map: Union[None, dict]

    load_raw: Union[None, bool]
    mapped_features: Union[None, str, bool]
    remove_gene_version: Union[None, bool]
    subset_gene_type: Union[None, str]
    streamlined_meta: bool

    sample_fn: Union[None, str]
    _sample_fns: Union[None, List[str]]

    _additional_annotation_key: Union[None, str]

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            load_func=None,
            dict_load_func_annotation=None,
            yaml_path: Union[str, None] = None,
            sample_fn: Union[str, None] = None,
            sample_fns: Union[List[str], None] = None,
            additional_annotation_key: Union[str, None] = None,
            **kwargs
    ):
        """

        :param data_path:
        :param meta_path:
        :param cache_path:
        :param load_func: Function to load data from disk into memory.

            Signature: load(data_dir, sample_fn, **kwargs)
        :param dict_load_func_annotation: Dictionary of functions to load additional observatino-wise annotation. The
            functions in the values of the dictionary can be selected via  self.additional_annotation_key which needs
            to correspond to a key of the dictionary.

            Signature: Dict[str, load_annotation(data_dir, sample_fn, additional_annotation_key, **kwargs)]
        :param yaml_path:
        :param sample_fn:
        :param sample_fns:
        :param additional_annotation_key: Key used by dict_load_func_annotation to identify which additional annotation
            is to be loaded.
        :param kwargs:
        """
        self._adata_ids = AdataIdsSfaira()
        self.ontology_container_sfaira = OCS  # Using a pre-instantiated version of this yields drastic speed-ups.

        self.adata = None
        self.meta = None
        self.genome = None
        self.data_dir_base = data_path
        self.meta_path = meta_path
        self.cache_path = cache_path

        self._author = None
        self._assay_sc = None
        self._assay_differentiation = None
        self._assay_type_differentiation = None
        self._bio_sample = None
        self._cell_line = None
        self._cell_type = None
        self._default_embedding = None
        self._development_stage = None
        self._disease = None
        self._doi_journal = None
        self._doi_preprint = None
        self._download_url_data = None
        self._download_url_meta = None
        self._ethnicity = None
        self._id = None
        self._individual = None
        self._ncells = None
        self._normalization = None
        self._organ = None
        self._organism = None
        self._primary_data = None
        self._sample_source = None
        self._sex = None
        self._source = None
        self._state_exact = None
        self._tech_sample = None
        self._title = None
        self._year = None

        self.assay_sc_obs_key = None
        self.assay_differentiation_obs_key = None
        self.assay_type_differentiation_obs_key = None
        self.bio_sample_obs_key = None
        self.cell_line_obs_key = None
        self.cell_type_obs_key = None
        self.development_stage_obs_key = None
        self.disease_obs_key = None
        self.ethnicity_obs_key = None
        self.individual_obs_key = None
        self.organ_obs_key = None
        self.organism_obs_key = None
        self.sample_source_obs_key = None
        self.sex_obs_key = None
        self.state_exact_obs_key = None
        self.tech_sample_obs_key = None

        self.gene_id_symbols_var_key = None
        self.gene_id_ensembl_var_key = None

        self.class_maps = {"0": {}}

        self._celltype_universe = None
        self._ontology_class_map = None

        self.load_raw = None
        self.mapped_features = None
        self.remove_gene_version = None
        self.subset_gene_type = None
        self.streamlined_meta = False

        self.sample_fn = sample_fn
        self._sample_fns = sample_fns

        # Check if YAML files exists, read meta data from there if available:
        if yaml_path is not None:
            assert os.path.exists(yaml_path), f"did not find yaml {yaml_path}"
            yaml_vals = read_yaml(fn=yaml_path)
            # Set organism first as this is required to disambiguate valid entries for other meta data.
            k = "organism"
            v = yaml_vals["attr"]["organism"]
            setattr(self, k, v)
            for k, v in yaml_vals["attr"].items():
                if v is not None and k not in ["organism", "sample_fns", "dataset_index"]:
                    if isinstance(v, dict):  # v is a dictionary over file-wise meta-data items
                        assert self.sample_fn in v.keys(), f"did not find key {self.sample_fn} in yamls keys for {k}"
                        v = v[self.sample_fn]
                    # Catches spelling errors in meta data definition (yaml keys).
                    if not hasattr(self, k) and not hasattr(self, "_" + k):
                        raise ValueError(f"Tried setting unavailable property {k}.")
                    try:
                        setattr(self, k, v)
                    except AttributeError as e:
                        raise ValueError(f"An error occured when setting {k} as {v}: {e}")
            # ID can be set now already because YAML was used as input instead of child class constructor.
            self.set_dataset_id(idx=yaml_vals["meta"]["dataset_index"])

        self.load_func = load_func
        self.dict_load_func_annotation = dict_load_func_annotation
        self._additional_annotation_key = additional_annotation_key

        self.supplier = "sfaira"

    @property
    def _directory_formatted_id(self) -> str:
        return "_".join("_".join(self.id.split("/")).split("."))

    def clear(self):
        """
        Remove loaded .adata to reduce memory footprint.

        :return:
        """
        import gc
        self.adata = None
        gc.collect()

    def download(self, **kwargs):
        assert self.download_url_data is not None, f"The `download_url_data` attribute of dataset {self.id} " \
                                                   f"is not set, cannot download dataset."
        assert self.data_dir_base is not None, "No path was provided when instantiating the dataset container, " \
                                               "cannot download datasets."

        if not os.path.exists(os.path.join(self.data_dir_base, self.directory_formatted_doi)):
            os.makedirs(os.path.join(self.data_dir_base, self.directory_formatted_doi))

        urls = self.download_url_data[0] + self.download_url_meta[0]

        for url in urls:
            if url is None:
                continue
            # Special case for data that is not publically available
            if url.split(",")[0] == 'private':
                if "," in url:
                    fn = ','.join(url.split(',')[1:])
                    if os.path.isfile(os.path.join(self.data_dir, fn)):
                        print(f"File {fn} already found on disk, skipping download.")
                    else:
                        warnings.warn(f"Dataset {self.id} is not available for automatic download, please manually "
                                      f"copy the file {fn} to the following location: "
                                      f"{self.data_dir}")
                else:
                    warnings.warn(f"A file for dataset {self.id} is not available for automatic download, please"
                                  f"manually copy the associated file to the following location: {self.data_dir}")
            # Special case for data from the synapse portal
            elif url.split(",")[0].startswith('syn'):
                fn = ",".join(url.split(",")[1:])
                if os.path.isfile(os.path.join(self.data_dir, fn)):
                    print(f"File {fn} already found on disk, skipping download.")
                else:
                    self._download_synapse(url.split(",")[0], fn, **kwargs)
            # Special case for public data that is labelled as not automatically downloadable
            elif url.split(",")[0] == 'manual':
                u = ",".join(url.split(",")[2:])
                fn = url.split(",")[1]
                if os.path.isfile(os.path.join(self.data_dir, fn)):
                    print(f"File {fn} already found on disk, skipping download.")
                else:
                    print(f"Data file {fn} for dataset {self.id} cannot be retrieved automatically. "
                          f"Please download it from {u} and copy to {os.path.join(self.data_dir, fn)}")
            # All other cases
            else:
                url = urllib.parse.unquote(url)
                try:
                    urllib.request.urlopen(url)
                except urllib.error.HTTPError as err:
                    # modify headers if urllib useragent is blocked (eg.10x datasets)
                    if err.code == 403:
                        opener = urllib.request.build_opener()
                        opener.addheaders = [('User-Agent', 'Mozilla/5.0 (Windows NT 6.1; WOW64)')]
                        urllib.request.install_opener(opener)
                except urllib.error.URLError:
                    # Catch SSLCertVerificationError: [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable
                    # to get local issuer certificate (_ssl.c:1124)
                    ssl._create_default_https_context = ssl._create_unverified_context

                if 'Content-Disposition' in urllib.request.urlopen(url).info().keys():
                    fn = cgi.parse_header(urllib.request.urlopen(url).info()['Content-Disposition'])[1]["filename"]
                else:
                    fn = url.split("/")[-1]
                # Only download if file not already downloaded:
                if not os.path.isfile(os.path.join(self.data_dir, fn)):
                    print(f"Downloading: {fn}")
                    urllib.request.urlretrieve(url, os.path.join(self.data_dir, fn))

    def _download_synapse(self, synapse_entity, fn, **kwargs):
        try:
            import synapseclient
        except ImportError:
            warnings.warn("synapseclient python package not found. This package is required to download some of the "
                          "selected datasets. Run `pip install synapseclient` to install it. Skipping download of the "
                          f"following dataset: {self.id}")
            return
        import shutil
        import logging
        logging.captureWarnings(False)  # required to properly display warning messages below with sypaseclient loaded

        if "synapse_user" not in kwargs.keys():
            warnings.warn(f"No synapse username provided, skipping download of synapse dataset {fn}."
                          f"Provide your synapse username as the `synapse_user` argument to the download method.")
            return
        if "synapse_pw" not in kwargs.keys():
            warnings.warn(f"No synapse password provided, skipping download of synapse dataset {fn}."
                          f"Provide your synapse password as the `synapse_pw` argument to the download method.")
            return

        print(f"Downloading from synapse: {fn}")
        syn = synapseclient.Synapse()
        syn.login(kwargs['synapse_user'], kwargs['synapse_pw'])
        dataset = syn.get(entity=synapse_entity)
        shutil.move(dataset.data_dir_base, os.path.join(self.data_dir, fn))

    @property
    def cache_fn(self):
        if self.directory_formatted_doi is None or self._directory_formatted_id is None:
            # TODO is this case necessary?
            warnings.warn("Caching enabled, but Dataset.id or Dataset.doi not set. Disabling caching for now.")
            return None
        else:
            if self.cache_path is None:
                cache = self.data_dir
            else:
                cache = os.path.join(self.cache_path, self.directory_formatted_doi)
            return os.path.join(cache, "cache", self._directory_formatted_id + ".h5ad")

    def load(
            self,
            load_raw: bool = False,
            allow_caching: bool = True,
            **kwargs
    ):
        """
        Load the selected datasets into memory.
        Cache is written into director named after doi and h5ad named after data set id.
        Cache is not over-written.

        :param load_raw: Force reading the raw object even when a cached one is present.
        :param allow_caching: Write the object to cache after loading if no cache exists yet.
        """
        # Sanity checks
        if self.adata is not None:
            raise ValueError(f"adata of {self.id} already loaded.")
        if self.data_dir is None:
            raise ValueError("No sfaira data repo path provided in constructor.")

        def _error_buffered_reading(**load_kwargs):
            self.adata = self.load_func(data_dir=self.data_dir, sample_fn=self.sample_fn, **load_kwargs)

        # Run data set-specific loading script:
        def _assembly_wrapper():
            if self.load_func is None:
                raise ValueError(f"Tried to access load_func for {self.id} but did not find any.")
            _error_buffered_reading(**kwargs)
            # Enable loading of additional annotation, e.g. secondary cell type annotation
            # The additional annotation `obs2 needs to be on a subset of the original annotation `self.adata.obs`.
            if self.dict_load_func_annotation is not None:
                obs2 = self.dict_load_func_annotation[self.additional_annotation_key](
                    data_dir=self.data_dir, sample_fn=self.sample_fn)
                assert np.all([x in self.adata.obs.index for x in obs2.index]), \
                    "index mismatch between additional annotation and original"
                self.adata = self.adata[obs2.index, :]
                # Overwrite annotation
                for k, v in obs2.items():
                    self.adata.obs[k] = v

        def _cached_reading(filename):
            if filename is not None:
                if os.path.exists(filename):
                    self.adata = anndata.read_h5ad(filename)
                else:
                    _error_buffered_reading()
            else:
                _error_buffered_reading()

        def _cached_writing(filename):
            if filename is not None:
                dir_cache = os.path.dirname(filename)
                if not os.path.exists(dir_cache):
                    os.makedirs(dir_cache)
                if not os.path.exists(filename) and self.adata is not None:
                    self.adata.write_h5ad(filename)

        if load_raw and allow_caching:
            _assembly_wrapper()
            _cached_writing(self.cache_fn)
        elif load_raw and not allow_caching:
            _assembly_wrapper()
        elif not load_raw and allow_caching:
            _cached_reading(self.cache_fn)
            _cached_writing(self.cache_fn)
        else:  # not load_raw and not allow_caching
            _cached_reading(self.cache_fn)

        # Set loading-specific metadata:
        self.load_raw = load_raw

    load.__doc__ = load_doc

    def _add_missing_featurenames(
            self,
            match_to_reference: Union[str, bool, None],
    ):
        if self.gene_id_symbols_var_key is None and self.gene_id_ensembl_var_key is None:
            raise ValueError("Either gene_id_symbols_var_key or gene_id_ensembl_var_key needs to be provided in the"
                             " dataloader")
        elif self.gene_id_symbols_var_key is None and self.gene_id_ensembl_var_key:
            # Convert ensembl ids to gene symbols
            id_dict = self.genome_container.id_to_symbols_dict
            ensids = self.adata.var.index if self.gene_id_ensembl_var_key == "index" else self.adata.var[self.gene_id_ensembl_var_key]
            self.adata.var[self._adata_ids.feature_symbol] = [
                id_dict[n.split(".")[0]] if n.split(".")[0] in id_dict.keys() else 'n/a'
                for n in ensids
            ]
            self.gene_id_symbols_var_key = self._adata_ids.feature_symbol
        elif self.gene_id_symbols_var_key and self.gene_id_ensembl_var_key is None:
            # Convert gene symbols to ensembl ids
            id_dict = self.genome_container.symbol_to_id_dict
            id_strip_dict = self.genome_container.strippednames_to_id_dict
            # Matching gene names to ensembl ids in the following way: if the gene is present in the ensembl dictionary,
            # match it straight away, if it is not in there we try to match everything in front of the first period in
            # the gene name with a dictionary that was modified in the same way, if there is still no match we append na
            ensids = []
            symbs = self.adata.var.index if self.gene_id_symbols_var_key == "index" else \
                self.adata.var[self.gene_id_symbols_var_key]
            for n in symbs:
                if n in id_dict.keys():
                    ensids.append(id_dict[n])
                elif n.split(".")[0] in id_strip_dict.keys():
                    ensids.append(id_strip_dict[n.split(".")[0]])
                else:
                    ensids.append('n/a')
            self.adata.var[self._adata_ids.feature_id] = ensids
            self.gene_id_ensembl_var_key = self._adata_ids.feature_id

    def _collapse_ensembl_gene_id_versions(self):
        """
        Remove version tag on ensembl gene ID so that different versions are superimposed downstream.

        :return:
        """
        if not self.gene_id_ensembl_var_key:
            raise ValueError(
                "Cannot remove gene version when gene_id_ensembl_var_key is not set in dataloader and "
                "match_to_reference is False"
            )
        elif self.gene_id_ensembl_var_key == "index":
            self.adata.index = [
                x.split(".")[0] for x in self.adata.var.index
            ]
        else:
            self.adata.var[self.gene_id_ensembl_var_key] = [
                x.split(".")[0] for x in self.adata.var[self.gene_id_ensembl_var_key].values
            ]
        # Collapse if necessary:
        self.adata = collapse_matrix(adata=self.adata, var_column=self.gene_id_ensembl_var_key)

    def collapse_counts(self):
        """
        Collapse count matrix along duplicated index.
        """
        if len(np.unique(self.adata.var.index)) < self.adata.var.shape[0]:
            self.adata = collapse_matrix(adata=self.adata, var_column="index")

    def streamline_features(
            self,
            match_to_reference: Union[str, Dict[str, str], None],
            remove_gene_version: bool = True,
            subset_genes_to_type: Union[None, str, List[str]] = None,
    ):
        """
        Subset and sort genes to genes defined in an assembly or genes of a particular type, such as protein coding.
        This also adds missing ensid or gene symbol columns if match_to_reference is not set to False and removes all
        adata.var columns that are not defined as gene_id_ensembl_var_key or gene_id_symbol_var_key in the dataloader.

        :param match_to_reference: Which annotation to map the feature space to. Can be:
                - str: Provide the name of the annotation in the format Organism.Assembly.Release
                - dict: Mapping of organism to name of the annotation (see str format). Chooses annotation for each
                    data set based on organism annotation.
        :param remove_gene_version: Whether to remove the version number after the colon sometimes found in ensembl
            gene ids.
        :param subset_genes_to_type: Type(s) to subset to. Can be a single type or a list of types or None.
            Types can be:
                - None: All genes in assembly.
                - "protein_coding": All protein coding genes in assembly.
        """
        self.__assert_loaded()

        # Set genome container if mapping of gene labels is requested
        if isinstance(match_to_reference, dict):
            match_to_reference = match_to_reference[self.organism]
        self._set_genome(assembly=match_to_reference)
        self.mapped_features = self.genome_container.assembly

        self.remove_gene_version = remove_gene_version
        self.subset_gene_type = subset_genes_to_type
        # Streamline feature space:
        self._add_missing_featurenames(match_to_reference=match_to_reference)
        for key in [self.gene_id_ensembl_var_key, self.gene_id_symbols_var_key]:
            # Make features unique (to avoid na-matches in converted columns to be collapsed by
            # _collapse_ensembl_gene_id_versions() below.
            if not key:
                pass
            elif key == "index":
                self.adata.var.index = make_index_unique(self.adata.var.index).tolist()
            else:
                self.adata.var[key] = make_index_unique(pd.Index(self.adata.var[key].values.tolist())).tolist()
        if remove_gene_version:
            self._collapse_ensembl_gene_id_versions()

        # Convert data matrix to csc matrix
        if isinstance(self.adata.X, np.ndarray):
            # Change NaN to zero. This occurs for example in concatenation of anndata instances.
            if np.any(np.isnan(self.adata.X)):
                self.adata.X[np.isnan(self.adata.X)] = 0
            x = scipy.sparse.csc_matrix(self.adata.X)
        elif isinstance(self.adata.X, scipy.sparse.spmatrix):
            x = self.adata.X.tocsc()
        else:
            raise ValueError(f"Data type {type(self.adata.X)} not recognized.")

        # Compute indices of genes to keep
        data_ids_ensg = self.adata.var.index.values if self.gene_id_ensembl_var_key == "index" \
            else self.adata.var[self.gene_id_ensembl_var_key].values
        if subset_genes_to_type is None:
            subset_ids_ensg = self.genome_container.ensembl
            subset_ids_symbol = self.genome_container.symbols
        else:
            if isinstance(subset_genes_to_type, str):
                subset_genes_to_type = [subset_genes_to_type]
            keys = np.unique(self.genome_container.biotype)
            if subset_genes_to_type not in keys:
                raise ValueError(f"subset type {subset_genes_to_type} not available in list {keys}")
            subset_ids_ensg = [
                x.upper() for x, y in zip(self.genome_container.ensembl, self.genome_container.biotype)
                if y in subset_genes_to_type
            ]
            subset_ids_symbol = [
                x.upper() for x, y in zip(self.genome_container.symbols, self.genome_container.biotype)
                if y in subset_genes_to_type
            ]

        # Remove unmapped genes
        idx_feature_kept = np.where([x.upper() in subset_ids_ensg for x in data_ids_ensg])[0]
        data_ids_kept = data_ids_ensg[idx_feature_kept]
        x = x[:, idx_feature_kept]
        # Build map of subset_ids to features in x:
        idx_feature_map = np.array([subset_ids_ensg.index(x) for x in data_ids_kept])
        # Create reordered feature matrix based on reference and convert to csr
        x_new = scipy.sparse.csc_matrix((x.shape[0], len(subset_ids_ensg)), dtype=x.dtype)
        # copying this over to the new matrix in chunks of size `steps` prevents a strange scipy error:
        # ... scipy/sparse/compressed.py", line 922, in _zero_many i, j, offsets)
        # ValueError: could not convert integer scalar
        step = 500
        if step < len(idx_feature_map):
            i = 0
            for i in range(0, len(idx_feature_map), step):
                x_new[:, idx_feature_map[i:i + step]] = x[:, i:i + step]
            x_new[:, idx_feature_map[i + step:]] = x[:, i + step:]
        else:
            x_new[:, idx_feature_map] = x
        x_new = x_new.tocsr()

        # Create new var dataframe
        if self.gene_id_symbols_var_key == "index":
            var_index = subset_ids_symbol
            var_data = {self.gene_id_ensembl_var_key: subset_ids_ensg}
        elif self.gene_id_ensembl_var_key == "index":
            var_index = subset_ids_ensg
            var_data = {self.gene_id_symbols_var_key: subset_ids_symbol}
        else:
            var_index = None
            var_data = {self.gene_id_symbols_var_key: subset_ids_symbol,
                        self.gene_id_ensembl_var_key: subset_ids_ensg}
        var_new = pd.DataFrame(data=var_data, index=var_index)

        self.adata = anndata.AnnData(
            X=x_new,
            obs=self.adata.obs,
            obsm=self.adata.obsm,
            var=var_new,
            uns=self.adata.uns
        )
        self.adata.uns[self._adata_ids.mapped_features] = match_to_reference

    def streamline_metadata(
            self,
            schema: str = "sfaira",
            clean_obs: bool = True,
            clean_var: bool = True,
            clean_uns: bool = True,
            clean_obs_names: bool = True,
            keep_orginal_obs: bool = False,
            keep_symbol_obs: bool = True,
            keep_id_obs: bool = True,
    ):
        """
        Streamline the adata instance to a defined output schema.

        Output format are saved in ADATA_FIELDS* classes.

        Note on ontology-controlled meta data:
        These are defined for a given format in `ADATA_FIELDS*.ontology_constrained`.
        They may appear in three different formats:
            - original (free text) annotation
            - ontology symbol
            - ontology ID
        During streamlining, these ontology-controlled meta data are projected to all of these three different formats.
        The initially annotated column may be any of these and is defined as "{attr}_obs_col".
        The resulting three column per meta data item are named:
            - ontology symbol: "{ADATA_FIELDS*.attr}"
            - ontology ID: {ADATA_FIELDS*.attr}_{ADATA_FIELDS*.onto_id_suffix}"
            - original (free text) annotation: "{ADATA_FIELDS*.attr}_{ADATA_FIELDS*.onto_original_suffix}"

        :param schema: Export format.
            - "sfaira"
            - "cellxgene"
        :param clean_obs: Whether to delete non-streamlined fields in .obs, .obsm and .obsp.
        :param clean_var: Whether to delete non-streamlined fields in .var, .varm and .varp.
        :param clean_uns: Whether to delete non-streamlined fields in .uns.
        :param clean_obs_names: Whether to replace obs_names with a string comprised of dataset id and an increasing
            integer.
        :param keep_orginal_obs: For ontology-constrained .obs columns, whether to keep a column with original
            annotation.
        :param keep_symbol_obs: For ontology-constrained .obs columns, whether to keep a column with ontology symbol
            annotation.
        :param keep_id_obs: For ontology-constrained .obs columns, whether to keep a column with ontology ID annotation.
        :return:
        """
        schema_version = schema.split(":")[-1] if ":" in schema else None
        self.__assert_loaded()

        # Set schema as provided by the user
        if schema.startswith("sfaira"):
            adata_target_ids = AdataIdsSfaira()
        elif schema.startswith("cellxgene"):
            if self.organism == "human":
                adata_target_ids = AdataIdsCellxgeneHuman_v1_1_0()
            elif self.organism == "human":
                adata_target_ids = AdataIdsCellxgeneHuman_v1_1_0()
            else:
                adata_target_ids = AdataIdsCellxgeneGeneral()
        else:
            raise ValueError(f"did not recognize schema {schema}")

        if hasattr(adata_target_ids, "gene_id_ensembl") and not hasattr(self._adata_ids, "gene_id_ensembl"):
            raise ValueError(f"Cannot convert this object to schema {schema}, as the currently applied schema does not "
                             f"have an ensembl gene ID annotation. Please run .streamline_features() first.")
        experiment_batch_labels = [getattr(self._adata_ids, x) for x in self._adata_ids.batch_keys]

        # Creating new var annotation
        var_new = pd.DataFrame()
        for k in adata_target_ids.var_keys:
            if k == "feature_id":
                if not self.gene_id_ensembl_var_key:
                    raise ValueError("feature_id not set in dataloader despite being required by the "
                                     "selected meta data schema. please run streamline_features() first to create the "
                                     "missing annotation")
                elif self.gene_id_ensembl_var_key == "index":
                    var_new[getattr(adata_target_ids, k)] = self.adata.var.index.tolist()
                else:
                    var_new[getattr(adata_target_ids, k)] = self.adata.var[self.gene_id_ensembl_var_key].tolist()
                    del self.adata.var[self.gene_id_ensembl_var_key]
                self.gene_id_ensembl_var_key = getattr(adata_target_ids, k)
            elif k == "feature_symbol":
                if not self.gene_id_symbols_var_key:
                    raise ValueError("gene_id_symbols_var_key not set in dataloader despite being required by the "
                                     "selected meta data schema. please run streamline_features() first to create the "
                                     "missing annotation")
                elif self.gene_id_symbols_var_key == "index":
                    var_new[getattr(adata_target_ids, k)] = self.adata.var.index.tolist()
                else:
                    var_new[getattr(adata_target_ids, k)] = self.adata.var[self.gene_id_symbols_var_key].tolist()
                    del self.adata.var[self.gene_id_symbols_var_key]
                self.gene_id_symbols_var_key = getattr(adata_target_ids, k)
            else:
                val = getattr(self, k)
                while hasattr(val, '__len__') and not isinstance(val, str) and len(val) == 1:  # unpack nested lists/tuples
                    val = val[0]
                var_new[getattr(adata_target_ids, k)] = val
        # set var index
        var_new.index = var_new[adata_target_ids.feature_index].tolist()
        if clean_var:
            if self.adata.varm is not None:
                del self.adata.varm
            if self.adata.varp is not None:
                del self.adata.varp
            self.adata.var = var_new
            if "feature_id" not in adata_target_ids.var_keys:
                self.gene_id_ensembl_var_key = None
            if "feature_symbol" not in adata_target_ids.var_keys:
                self.gene_id_symbols_var_key = None
        else:
            index_old = self.adata.var.index.copy()
            # Add old columns in if they are not duplicated:
            self.adata.var = pd.concat([
                var_new,
                pd.DataFrame(dict([(k, v) for k, v in self.adata.var.items() if k not in var_new.columns]))
            ], axis=1)
            self.adata.var.index = index_old

        # Prepare new .uns dict:
        uns_new = {}
        for k in adata_target_ids.uns_keys:
            if hasattr(self, k) and getattr(self, k) is not None:
                val = getattr(self, k)
            elif hasattr(self, f"{k}_obs_key") and getattr(self, f"{k}_obs_key") is not None:
                val = np.sort(np.unique(self.adata.obs[getattr(self, f"{k}_obs_key")].values)).tolist()
            elif getattr(self._adata_ids, k) in self.adata.obs.columns:
                val = np.sort(np.unique(self.adata.obs[getattr(self._adata_ids, k)].values)).tolist()
            else:
                val = None
            while hasattr(val, '__len__') and not isinstance(val, str) and len(val) == 1:  # Unpack nested lists/tuples.
                val = val[0]
            uns_new[getattr(adata_target_ids, k)] = val
        if clean_uns:
            self.adata.uns = uns_new
        else:
            self.adata.uns.update(uns_new)

        # Prepare new .obs dataframe
        # Queried meta data may be:
        # 1) in .obs
        #   a) for an ontology-constrained meta data item
        #       I) as free annotation with a term map to an ontology
        #       II) as column with ontology symbols
        #       III) as column with ontology IDs
        #   b) for a non-ontology-constrained meta data item:
        #       I) as free annotation
        # 2) in .uns
        #   b) as elements that are ontology symbols
        #   c) as elements that are ontology IDs
        # .obs annotation takes priority over .uns annotation if both are present.
        # The output columns are:
        # - for an ontology-constrained meta data item "attr":
        #   *  symbols:         "attr"
        #   *  IDs:             "attr" + self._adata_ids.onto_id_suffix
        #   *  original labels: "attr" + self._adata_ids.onto_original_suffix
        # - for a non-ontology-constrained meta data item "attr":
        #   *  original labels: "attr" + self._adata_ids.onto_original_suffix
        obs_new = pd.DataFrame(index=self.adata.obs.index)
        for k in [x for x in adata_target_ids.obs_keys]:
            if k in experiment_batch_labels and getattr(self, f"{k}_obs_key") is not None and \
                    "*" in getattr(self, f"{k}_obs_key"):
                # Handle batch-annotation columns which can be provided as a combination of columns separated by an
                # asterisk.
                # The queried meta data are always:
                # 1b-I) a combination of existing columns in .obs
                old_cols = getattr(self, f"{k}_obs_key")
                batch_cols = []
                for batch_col in old_cols.split("*"):
                    if batch_col in self.adata.obs_keys():
                        batch_cols.append(batch_col)
                    else:
                        # This should not occur in single data set loaders (see warning below) but can occur in
                        # streamlined data loaders if not all instances of the streamlined data sets have all columns
                        # in .obs set.
                        print(f"WARNING: attribute {batch_col} of data set {self.id} was not found in columns.")
                # Build a combination label out of all columns used to describe this group.
                val = [
                    "_".join([str(xxx) for xxx in xx])
                    for xx in zip(*[self.adata.obs[batch_col].values.tolist() for batch_col in batch_cols])
                ]
            else:
                # Locate annotation.
                if hasattr(self, f"{k}_obs_key") and getattr(self, f"{k}_obs_key") is not None and \
                        getattr(self, f"{k}_obs_key") in self.adata.obs.columns:
                    # Last and-clause to check if this column is included in data sets. This may be violated if data
                    # is obtained from a database which is not fully streamlined.
                    # Look for 1a-* and 1b-I
                    val = self.adata.obs[getattr(self, f"{k}_obs_key")].values.tolist()
                else:
                    # Look for 2a, 2b
                    val = getattr(self, k)
                    if val is None:
                        val = self._adata_ids.unknown_metadata_identifier
                    # Unpack nested lists/tuples:
                    while hasattr(val, '__len__') and not isinstance(val, str) and len(val) == 1:
                        val = val[0]
                    val = [val] * self.adata.n_obs
            # Identify annotation: disambiguate 1a-I, 1a-II, 1a-III, 1b-I.
            if k in self._adata_ids.ontology_constrained:
                # 1a-*.
                if isinstance(self.get_ontology(k=k), OntologyHierarchical) and np.all([
                    self.get_ontology(k=k).is_a_node_name(x) or x == self._adata_ids.unknown_metadata_identifier
                    for x in np.unique(val)
                ]):  # 1a-II)
                    new_col = getattr(adata_target_ids, k)
                    validation_ontology = self.get_ontology(k=k)
                elif isinstance(self.get_ontology(k=k), OntologyHierarchical) and np.all([
                    self.get_ontology(k=k).is_a_node_id(x) or x == self._adata_ids.unknown_metadata_identifier
                    for x in np.unique(val)
                ]):  # 1a-III)
                    new_col = getattr(adata_target_ids, k) + self._adata_ids.onto_id_suffix
                    validation_ontology = None
                else:  # 1a-I)
                    new_col = getattr(adata_target_ids, k) + self._adata_ids.onto_original_suffix
                    validation_ontology = None
            else:
                # 1b-I.
                new_col = getattr(adata_target_ids, k)
                validation_ontology = self.get_ontology(k=k)
            # Check values for validity:
            self._value_protection(attr=new_col, allowed=validation_ontology, attempted=[
                x for x in np.unique(val)
                if x not in [
                    self._adata_ids.unknown_metadata_identifier,
                ]
            ])
            obs_new[new_col] = val
            # For ontology-constrained meta data, the remaining columns are added after .obs cleaning below.
        if clean_obs:
            if self.adata.obsm is not None:
                del self.adata.obsm
            if self.adata.obsp is not None:
                del self.adata.obsp
            self.adata.obs = obs_new
        else:
            index_old = self.adata.obs.index.copy()
            # Add old columns in if they are not duplicated in target obs column space, even if this column is not
            # defined. This would result in the instance accessing this column assuming it was streamlined.
            self.adata.obs = pd.concat([
                obs_new,
                pd.DataFrame(dict([(k, v) for k, v in self.adata.obs.items()
                                   if k not in adata_target_ids.controlled_meta_keys]))
            ], axis=1)
            self.adata.obs.index = index_old
        for k in [x for x in adata_target_ids.obs_keys if x in adata_target_ids.ontology_constrained]:
            # Add remaining output columns for ontology-constrained meta data.
            self.__impute_ontology_cols_obs(attr=k, adata_ids=adata_target_ids)
            # Delete attribute-specific columns that are not desired.
            col_name = getattr(self._adata_ids, k) + self._adata_ids.onto_id_suffix
            if not keep_id_obs and col_name in self.adata.obs.columns:
                del self.adata.obs[col_name]
            col_name = getattr(self._adata_ids, k) + self._adata_ids.onto_original_suffix
            if not keep_orginal_obs and col_name in self.adata.obs.columns:
                del self.adata.obs[col_name]
            col_name = getattr(self._adata_ids, k)
            if not keep_symbol_obs and col_name in self.adata.obs.columns:
                del self.adata.obs[col_name]
        if clean_obs_names:
            self.adata.obs.index = [f"{self.id}_{i}" for i in range(1, self.adata.n_obs + 1)]

        # Make sure that correct unknown_metadata_identifier is used in .uns, .obs and .var metadata
        unknown_old = self._adata_ids.unknown_metadata_identifier
        unknown_new = adata_target_ids.unknown_metadata_identifier
        self.adata.obs = self.adata.obs.replace({None: unknown_new})
        self.adata.obs = self.adata.obs.replace({unknown_old: unknown_new})
        self.adata.var = self.adata.var.replace({None: unknown_new})
        self.adata.var = self.adata.var.replace({unknown_old: unknown_new})
        for k in self.adata.uns_keys():
            if self.adata.uns[k] is None or self.adata.uns[k] == unknown_old:
                self.adata.uns[k] = unknown_new

        self._adata_ids = adata_target_ids  # set new adata fields to class after conversion
        self.streamlined_meta = True
        # Add additional hard-coded description changes for cellxgene schema:
        if schema.startswith("cellxgene"):
            self.adata = cellxgene_export_adaptor(adata=self.adata, adata_ids=self._adata_ids, version=schema_version)

    def write_distributed_store(
            self,
            dir_cache: Union[str, os.PathLike],
            store_format: str = "dao",
            dense: bool = False,
            compression_kwargs: dict = {},
            chunks: Union[int, None] = None,
    ):
        """
        Write data set into a format that allows distributed access to data set on disk.

        Stores are useful for distributed access to data sets, in many settings this requires some streamlining of the
        data sets that are accessed. Use .streamline_* before calling this method to streamline the data sets.

        :param dir_cache: Directory to write cache in.
        :param store_format: Disk format for objects in cache. Recommended is "dao".

            - "h5ad": Allows access via backed .h5ad.
                Note on compression: .h5ad supports sparse data with is a good compression that gives fast row-wise
                    access if the files are csr, so further compression potentially not necessary.
            - "dao": Distributed access optimised format, recommended for batched access in optimisation, for example.
        :param dense: Whether to write sparse or dense store, this will be homogenously enforced.
        :param compression_kwargs: Compression key word arguments to give to h5py or zarr
            For store_format=="h5ad", see also anndata.AnnData.write_h5ad:
                - compression,
                - compression_opts.
            For store_format=="dao", see also sfaira.data.write_dao which relays kwargs to
            zarr.hierarchy.create_dataset:
                - compressor
                - overwrite
                - order
                and others.
        :param chunks: Observation axes of chunk size of zarr array, see anndata.AnnData.write_zarr documentation.
            Only relevant for store=="dao". The feature dimension of the chunks is always is the full feature space.
            Uses zarr default chunking across both axes if None.
        """
        self.__assert_loaded()
        if store_format == "h5ad":
            if not isinstance(self.adata.X, scipy.sparse.csr_matrix):
                print(f"WARNING: high-perfomances caches based on .h5ad work better with .csr formatted expression "
                      f"data, found {type(self.adata.X)}")
            fn = os.path.join(dir_cache, self.doi_cleaned_id + ".h5ad")
            as_dense = ("X",) if dense else ()
            print(f"writing {self.adata.shape} into {fn}")
            self.adata.write_h5ad(filename=fn, as_dense=as_dense, **compression_kwargs)
        elif store_format == "dao":
            # Convert data object to sparse / dense as required:
            if not dense:
                raise ValueError("WARNING: sparse zarr array performance is not be optimal and not supported yet, "
                                 "consider writing as dense and consider that zarr arrays are compressed on disk!")
            fn = os.path.join(dir_cache, self.doi_cleaned_id)
            chunks = (chunks, self.adata.X.shape[1]) if chunks is not None else True
            write_dao(store=fn, adata=self.adata, chunks=chunks, compression_kwargs=compression_kwargs)
        else:
            raise ValueError()

    def write_backed(
            self,
            adata_backed: anndata.AnnData,
            genome: str,
            idx: np.ndarray,
            load_raw: bool = False,
            allow_caching: bool = True
    ):
        """
        Loads data set into slice of backed anndata object.

        Note: scatter updates to backed sparse arrays are not yet supported by anndata. Accordingly, we need to work
        around below using .append() of the backed matrix.

        :param adata_backed:
        :param genome: Genome name to use as refernce.
        :param idx: Indices in adata_backed to write observations to. This can be used to immediately create a
            shuffled object.
        :param load_raw: See .load().
        :param allow_caching: See .load().
        :return: New row index for next element to be written into backed anndata.
        """
        self.load(load_raw=load_raw, allow_caching=allow_caching)
        # Check if writing to sparse or dense matrix:
        if isinstance(adata_backed.X, np.ndarray) or \
                isinstance(adata_backed.X, h5py._hl.dataset.Dataset):  # backed dense
            if isinstance(self.adata.X, scipy.sparse.csr_matrix) or \
                    isinstance(self.adata.X, scipy.sparse.csc_matrix) or \
                    isinstance(self.adata.X, scipy.sparse.lil_matrix):
                # map to dense array
                x_new = self.adata.X.toarray()
            else:
                x_new = self.adata.X

            adata_backed.X[np.sort(idx), :] = x_new[np.argsort(idx), :]
            for k in adata_backed.obs.columns:
                if k == self._adata_ids.dataset:
                    adata_backed.obs.loc[np.sort(idx), self._adata_ids.dataset] = [
                        self.id for _ in range(len(idx))]
                elif k in self.adata.obs.columns:
                    adata_backed.obs.loc[np.sort(idx), k] = self.adata.obs[k].values[np.argsort(idx)]
                elif k in list(self.adata.uns.keys()):
                    adata_backed.obs.loc[np.sort(idx), k] = [self.adata.uns[k] for i in range(len(idx))]
                else:
                    # Need to fill this instead of throwing an exception as this condition can trigger for one element
                    # within a loop over multiple data sets (ie in data set human).
                    adata_backed.obs.loc[idx, k] = ["key_not_found" for i in range(len(idx))]
        elif isinstance(adata_backed.X, anndata._core.sparse_dataset.SparseDataset):  # backed sparse
            # cannot scatter update on backed sparse yet! assert that updated block is meant to be appended:
            assert np.all(idx == np.arange(adata_backed.shape[0], adata_backed.shape[0] + len(idx)))
            if not isinstance(self.adata.X, scipy.sparse.csr_matrix):
                x_new = self.adata.X.tocsr()
            else:
                x_new = self.adata.X
            adata_backed.X.append(x_new[np.argsort(idx)])
            adata_backed._n_obs = adata_backed.X.shape[0]  # not automatically updated after append
            adata_backed.obs = adata_backed.obs.append(  # .obs was not broadcasted to the right shape!
                pandas.DataFrame(dict([
                    (k, [self.id for i in range(len(idx))]) if k == self._adata_ids.dataset
                    else (k, self.adata.obs[k].values[np.argsort(idx)]) if k in self.adata.obs.columns
                    else (k, [self.adata.uns[k] for _ in range(len(idx))]) if k in list(self.adata.uns.keys())
                    else (k, ["key_not_found" for _ in range(len(idx))])
                    for k in adata_backed.obs.columns
                ]))
            )
            self.clear()
        else:
            raise ValueError(f"Did not recognize backed AnnData.X format {type(adata_backed.X)}")

    def _set_genome(self, assembly: Union[str, None]):
        self.genome_container = GenomeContainer(
            assembly=assembly,
        )

    @property
    def doi_cleaned_id(self):
        return "_".join(self.id.split("_")[:-1])

    def get_ontology(self, k) -> OntologyHierarchical:
        x = getattr(self.ontology_container_sfaira, k) if hasattr(self.ontology_container_sfaira, k) else None
        if isinstance(x, dict):
            assert isinstance(self.organism, str)
            x = x[self.organism]
        return x

    @property
    def fn_ontology_class_map_tsv(self):
        """Standardised file name under which cell type conversion tables are saved."""
        return self.doi_cleaned_id + ".tsv"

    def _write_ontology_class_map(self, fn, tab: pd.DataFrame):
        """
        Write class map to file.

        Helper to allow direct interaction with written table instead of using table from instance.

        :param fn: File name of csv to write class maps to.
        :param tab: Class map table.
        """
        tab.to_csv(fn, index=False, sep="\t")

    def write_ontology_class_map(
            self,
            fn,
            protected_writing: bool = True,
            **kwargs
    ):
        """
        Load class maps of free text cell types to ontology classes.

        :param fn: File name of tsv to write class maps to.
        :param protected_writing: Only write if file was not already found.
        :return:
        """
        if not self.annotated:
            warnings.warn(f"attempted to write ontology class maps for data set {self.id} without annotation")
        else:
            labels_original = np.sort(np.unique(self.adata.obs[self.cell_type_obs_key].values))
            tab = self.celltypes_universe.prepare_celltype_map_tab(
                source=labels_original,
                match_only=False,
                anatomical_constraint=self.organ,
                include_synonyms=True,
                omit_list=[self._adata_ids.unknown_metadata_identifier],
                **kwargs
            )
            if not os.path.exists(fn) or not protected_writing:
                self._write_ontology_class_map(fn=fn, tab=tab)

    def _read_ontology_class_map(self, fn) -> pd.DataFrame:
        """
        Read class map.

        Helper to allow direct interaction with resulting table instead of loading into instance.

        :param fn: File name of csv to load class maps from.
        :return:
        """
        try:
            # Need dtype="str" to force numeric cell type identifiers, e.g. cluster numbers to be in string format.
            tab = pd.read_csv(fn, header=0, index_col=None, sep="\t", dtype="str")
        except pandas.errors.ParserError as e:
            print(f"{self.id}")
            raise pandas.errors.ParserError(e)
        return tab

    def read_ontology_class_map(self, fn):
        """
        Load class maps of free text cell types to ontology classes.

        :param fn: File name of csv to load class maps from.
        :return:
        """
        if os.path.exists(fn):
            self.cell_type_map = self._read_ontology_class_map(fn=fn)
        else:
            if self.cell_type_obs_key is not None:
                warnings.warn(f"file {fn} does not exist but cell_type_obs_key {self.cell_type_obs_key} is given")

    def project_free_to_ontology(self, attr: str, copy: bool = False):
        """
        Project free text cell type names to ontology based on mapping table.

        ToDo: add ontology ID setting here.
        ToDo: only for cell type right now, extend to other meta data in the future.

        :param copy: If True, a dataframe with the celltype annotation is returned, otherwise self.adata.obs is updated
            inplace.

        :return:
        """
        ontology_map = attr + "_map"
        if hasattr(self, ontology_map):
            ontology_map = getattr(self, ontology_map)
        else:
            ontology_map = None
            print(f"WARNING: did not find ontology map for {attr} which was only defined by free annotation")
        adata_fields = self._adata_ids
        results = {}
        col_original = attr + adata_fields.onto_original_suffix
        labels_original = self.adata.obs[col_original].values
        if ontology_map is not None:  # only if this was defined
            labels_mapped = [
                ontology_map[x] if x in ontology_map.keys()
                else x for x in labels_original
            ]
            # Convert unknown celltype placeholders (needs to be hardcoded here as placeholders are also hardcoded in
            # conversion tsv files
            placeholder_conversion = {
                "UNKNOWN": adata_fields.unknown_metadata_identifier,
                "NOT_A_CELL": adata_fields.not_a_cell_celltype_identifier,
            }
            labels_mapped = [
                placeholder_conversion[x] if x in placeholder_conversion.keys()
                else x for x in labels_mapped
            ]
            map_exceptions = [adata_fields.unknown_metadata_identifier]
            if attr == "cell_type":
                map_exceptions.append(adata_fields.not_a_cell_celltype_identifier)
            # Validate mapped IDs based on ontology:
            # This aborts with a readable error if there was a target in the mapping file that doesnt match the ontology
            # This protection blocks progression in the unit test if not deactivated.
            self._value_protection(
                attr=attr,
                allowed=getattr(self.ontology_container_sfaira, attr),
                attempted=[x for x in list(set(labels_mapped)) if x not in map_exceptions],
            )
            # Add cell type IDs into object:
            # The IDs are not read from a source file but inferred based on the class name.
            # TODO this could be changed in the future, this allows this function to be used both on cell type name
            #  mapping files with and without the ID in the third column.
            # This mapping blocks progression in the unit test if not deactivated.
            results[getattr(adata_fields, attr)] = labels_mapped
            self.__project_ontology_ids_obs(attr=attr, map_exceptions=map_exceptions, from_id=False,
                                            adata_ids=adata_fields)
        else:
            results[getattr(adata_fields, attr)] = labels_original
            results[getattr(adata_fields, attr) + adata_fields.onto_id_suffix] = \
                [adata_fields.unknown_metadata_identifier] * self.adata.n_obs
        results[getattr(adata_fields, attr) + adata_fields.onto_original_suffix] = labels_original
        if copy:
            return pd.DataFrame(results, index=self.adata.obs.index)
        else:
            for k, v in results.items():
                self.adata.obs[k] = v

    def __impute_ontology_cols_obs(
            self,
            attr: str,
            adata_ids: AdataIds,
    ):
        """
        Add missing ontology defined columns (symbol, ID, original) for a given ontology.

        1) If original column is non-empty and symbol and ID are empty:
            orginal column is projected to ontology and both symbol and ID are inferred.
            Note that in this case, a label map is required.
        2) If ID column is non-empty or symbol is non-empty, an error is thrown.
            a) If ID column is non-empty and symbol is empty, symbol is inferred.
            b) If ID column is empty and symbol is non-empty, ID is inferred.
            c) If ID column is non-empty and non-symbol is empty, symbol is inferred and over-written.
                Note that this setting allows usage of data sets which were streamlined with a different ontology
                version.
            In all cases original is kept if it is set and is set to symbol otherwise.
        3) If original, ID and symbol columns are empty, no action is taken (meta data item was not set).
        """
        ontology = self.get_ontology(k=attr)
        col_symbol = getattr(adata_ids, attr)
        col_id = getattr(adata_ids, attr) + self._adata_ids.onto_id_suffix
        col_original = getattr(adata_ids, attr) + self._adata_ids.onto_original_suffix
        if ontology is None:
            # Fill with invalid ontology identifiers if no ontology was found.
            self.adata.obs[col_id] = \
                [self._adata_ids.invalid_metadata_identifier for _ in range(self.adata.n_obs)]
            self.adata.obs[col_original] = \
                [self._adata_ids.invalid_metadata_identifier for _ in range(self.adata.n_obs)]
            self.adata.obs[col_symbol] = \
                [self._adata_ids.invalid_metadata_identifier for _ in range(self.adata.n_obs)]
        else:
            # Note that for symbol and ID, the columns may be filled but not streamlined according to the ontology,
            # in that case the corresponding meta data is defined as absent.
            # Check which level of meta data annotation is present.
            # Symbols:
            symbol_col_present = col_symbol in self.adata.obs.columns
            symbol_col_streamlined = np.all([
                ontology.is_a_node_name(x) or x == self._adata_ids.unknown_metadata_identifier
                for x in np.unique(self.adata.obs[col_symbol].values)]) if symbol_col_present else False
            symbol_present = symbol_col_present and symbol_col_streamlined
            # IDs:
            id_col_present = col_id in self.adata.obs.columns
            id_col_streamlined = np.all([
                ontology.is_a_node_id(x) or x == self._adata_ids.unknown_metadata_identifier
                for x in np.unique(self.adata.obs[col_id].values)]) if id_col_present else False
            id_present = id_col_present and id_col_streamlined
            # Original annotation (free text):
            original_present = col_original in self.adata.obs.columns
            if original_present and not symbol_present and not id_present:  # 1)
                self.project_free_to_ontology(attr=attr, copy=False)
            if symbol_present or id_present:  # 2)
                if symbol_present and not id_present:  # 2a)
                    self.__project_ontology_ids_obs(attr=attr, from_id=False, adata_ids=adata_ids)
                if not symbol_present and id_present:  # 2b)
                    self.__project_ontology_ids_obs(attr=attr, from_id=True, adata_ids=adata_ids)
                if symbol_present and id_present:  # 2c)
                    self.__project_ontology_ids_obs(attr=attr, from_id=True, adata_ids=adata_ids)
                if not original_present:
                    val = self.adata.obs[col_symbol]
                    self.adata.obs[col_original] = val

    def __project_ontology_ids_obs(
            self,
            attr: str,
            adata_ids: AdataIds,
            map_exceptions: Union[None, List[str]] = None,
            map_exceptions_value=None,
            from_id: bool = False,
    ):
        """
        Project ontology names to IDs for a given ontology in .obs entries.

        :param ontology: ontology to use when converting to IDs
        :param attr: name of obs_column containing names to convert or python list containing these values
        :param map_exceptions: list of values that should not be mapped.
            Defaults to unknown meta data identifier defined in ID object if None.
        :param map_exceptions_value: placeholder target value for values excluded from mapping.
            Defaults to unknown meta data identifier defined in ID object if None.
        :param from_id: Whether to output ontology symbol or ID.
        :return:
        """
        ontology = self.get_ontology(k=attr)
        assert ontology is not None, f"cannot project value for {attr} because ontology is None"
        assert isinstance(attr, (str, list)), f"argument key_in needs to be of type str or list. Supplied" \
                                              f"type: {type(attr)}"
        map_exceptions = map_exceptions if map_exceptions is not None else [adata_ids.unknown_metadata_identifier]
        map_exceptions = [x.lower() for x in map_exceptions]
        if map_exceptions_value is None:
            # TODO this may be simplified in the future once all unknown meta data labels are the same.
            if attr == "cell_type":
                map_exceptions_value = adata_ids.unknown_metadata_identifier
            else:
                map_exceptions_value = adata_ids.unknown_metadata_identifier
        col_name = getattr(adata_ids, attr)
        if from_id:
            col_name += adata_ids.onto_id_suffix
        input_values = self.adata.obs[col_name].values
        map_vals = dict([
            (x, ontology.convert_to_name(x)) if from_id else
            (x, ontology.convert_to_id(x))
            for x in np.unique([
                xx for xx in input_values
                if (xx.lower() not in map_exceptions and xx is not None)
            ])
        ])
        output_values = [
            map_vals[x] if x in map_vals.keys() else map_exceptions_value
            for x in input_values
        ]
        key_out = getattr(adata_ids, attr) if from_id else getattr(adata_ids, attr) + adata_ids.onto_id_suffix
        self.adata.obs[key_out] = output_values

    @property
    def citation(self):
        """
        Return all information necessary to cite data set.

        :return:
        """
        return [self.author, self.year, self.doi_journal]

    # Meta data handling code: Reading, writing and selected properties. Properties are either set in constructor
    # (and saved in self._somename) or accessed in self.meta.

    @property
    def meta_fn(self):
        if self.meta_path is None:
            meta = self.data_dir
        else:
            meta = os.path.join(self.meta_path, self.directory_formatted_doi)
        if meta is None:
            return None
        else:
            return os.path.join(meta, "meta", self.doi_cleaned_id + "_meta.csv")

    def load_meta(self, fn: Union[PathLike, str, None]):
        if fn is None:
            if self.meta_fn is not None:
                fn = self.meta_fn
        else:
            if isinstance(fn, str):
                fn = os.path.normpath(fn)
        # Only load meta data if file exists:
        if fn is not None and os.path.isfile(fn):
            meta = pandas.read_csv(
                fn,
                usecols=list(META_DATA_FIELDS.keys()),
            )
            # using dtype in read_csv through errors some times.
            for k, v in META_DATA_FIELDS.items():
                if k in meta.columns:
                    if meta[k].values[0] is not None:
                        meta[k] = np.asarray(meta[k].values, dtype=v)
            self.meta = meta.fillna("None").replace({"None": None})

    def write_meta(
            self,
            fn_meta: Union[None, str] = None,
            dir_out: Union[None, str] = None,
    ):
        """
        Write meta data object for data set.

        Does not cache data and attempts to load raw data.

        :param fn_meta: File to write to, selects automatically based on self.meta_path and self.id otherwise.
        :param dir_out: Path to write to, file name is selected automatically based on self.id.
        :return:
        """
        if fn_meta is not None and dir_out is not None:
            raise ValueError("supply either fn_meta or dir_out but not both")
        elif fn_meta is None and dir_out is None:
            if self.meta_fn is None:
                raise ValueError("provide either fn in load or via constructor (meta_path)")
            fn_meta = self.meta_fn
        elif fn_meta is None and dir_out is not None:
            fn_meta = os.path.join(dir_out, self.doi_cleaned_id + "_meta.csv")
        elif fn_meta is not None and dir_out is None:
            pass  # fn_meta is used
        else:
            assert False, "bug in switch"

        if self.adata is None:
            self.load(load_raw=True, allow_caching=False)
        # Add data-set wise meta data into table:
        meta = pandas.DataFrame(index=range(1))
        # Expand table by variably cell-wise or data set-wise meta data:
        for x in self._adata_ids.controlled_meta_fields:
            if hasattr(self, f"{x}_obs_key") and getattr(self, f"{x}_obs_key") is not None:
                col = getattr(self._adata_ids, x)
                meta[col] = (self.adata.obs[col].unique(), )
                if x in self._adata_ids.ontology_constrained:
                    col = getattr(self._adata_ids, x + self._adata_ids.onto_id_suffix)
                    meta[col] = (self.adata.obs[col].unique(), )
            else:
                meta[getattr(self._adata_ids, x)] = getattr(self, x)
        meta.to_csv(fn_meta)

    def set_dataset_id(
            self,
            idx: int = 1
    ):
        if self.sample_fn is not None:
            idx += self._sample_fns.index(self.sample_fn)
        idx = str(idx).zfill(3)

        if isinstance(self.author, List):
            author = self.author[0]
        else:
            author = self.author

        # Note: access private attributes here, e.g. _organism, to avoid loading of content via meta data, which would
        # invoke call to self.id before it is set.
        self.id = f"{clean_id_str(self._organism)}_" \
                  f"{clean_id_str(self._organ)}_" \
                  f"{self._year}_" \
                  f"{clean_id_str(self._assay_sc)}_" \
                  f"{clean_id_str(author)}_" \
                  f"{idx}_" \
                  f"{self.doi_main}"

    # Properties:

    @property
    def additional_annotation_key(self) -> Union[None, str]:
        return self._additional_annotation_key

    @additional_annotation_key.setter
    def additional_annotation_key(self, x: str):
        self._additional_annotation_key = x

    @property
    def annotated(self) -> Union[bool, None]:
        if self.cell_type_obs_key is not None:
            return True
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.annotated in self.meta.columns:
                return self.meta[self._adata_ids.annotated].values[0]
            elif self.loaded:
                # If data set was loaded and there is still no annotation indicated, it is declared unannotated.
                return False
            else:
                # If data set was not yet loaded, it is unclear if annotation would be loaded in ._load(),
                # if also no meta data is available, we do not know the status of the data set.
                return None

    @property
    def assay_sc(self) -> Union[None, str]:
        if self._assay_sc is not None:
            return self._assay_sc
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.assay_sc in self.meta.columns:
                return self.meta[self._adata_ids.assay_sc]
            else:
                return None

    @assay_sc.setter
    def assay_sc(self, x: str):
        x = self._value_protection(attr="assay_sc", allowed=self.ontology_container_sfaira.assay_sc, attempted=x)
        self._assay_sc = x

    @property
    def assay_differentiation(self) -> Union[None, str]:
        if self._assay_differentiation is not None:
            return self._assay_differentiation
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.assay_differentiation in self.meta.columns:
                return self.meta[self._adata_ids.assay_differentiation]
            else:
                return None

    @assay_differentiation.setter
    def assay_differentiation(self, x: str):
        x = self._value_protection(attr="assay_differentiation",
                                   allowed=self.ontology_container_sfaira.assay_differentiation, attempted=x)
        self._assay_differentiation = x

    @property
    def assay_type_differentiation(self) -> Union[None, str]:
        if self._assay_type_differentiation is not None:
            return self._assay_type_differentiation
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.assay_type_differentiation in self.meta.columns:
                return self.meta[self._adata_ids.assay_type_differentiation]
            else:
                return None

    @assay_type_differentiation.setter
    def assay_type_differentiation(self, x: str):
        x = self._value_protection(attr="assay_type_differentiation",
                                   allowed=self.ontology_container_sfaira.assay_type_differentiation, attempted=x)
        self._assay_type_differentiation = x

    @property
    def bio_sample(self) -> Union[None, str]:
        if self._bio_sample is not None:
            return self._bio_sample
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.bio_sample in self.meta.columns:
                return self.meta[self._adata_ids.bio_sample]
            else:
                return None

    @bio_sample.setter
    def bio_sample(self, x: str):
        self._bio_sample = x

    @property
    def cell_line(self) -> Union[None, str]:
        if self._cell_line is not None:
            return self._cell_line
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.cell_line in self.meta.columns:
                return self.meta[self._adata_ids.cell_line]
            else:
                return None

    @cell_line.setter
    def cell_line(self, x: str):
        self._cell_line = x

    @property
    def cell_type(self) -> Union[None, str]:
        if self._cell_type is not None:
            return self._cell_type
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.cell_type in self.meta.columns:
                return self.meta[self._adata_ids.cell_type]
            else:
                return None

    @cell_type.setter
    def cell_type(self, x: str):
        self._cell_type = x

    @property
    def data_dir(self):
        # Data is either directly in user supplied directory or in a sub directory if the overall directory is managed
        # by sfaira: In this case, the sub directory is named after the doi of the data set.
        if self.data_dir_base is None:
            return None
        else:
            sfaira_path = os.path.join(self.data_dir_base, self.directory_formatted_doi)
            # Allow checking in secondary path, named after second DOI associated with study.
            # This allows association of raw data already downloaded even after DOI is updated.
            if self.doi_preprint is not None:
                sfaira_path_secondary = os.path.join(self.data_dir_base,
                                                     get_directory_formatted_doi(x=self.doi_preprint))
            else:
                sfaira_path_secondary = None
            if os.path.exists(sfaira_path):
                return sfaira_path
            elif self.doi_preprint is not None and os.path.exists(sfaira_path_secondary):
                return sfaira_path_secondary
            else:
                return self.data_dir_base

    @property
    def default_embedding(self) -> Union[None, str]:
        if self._default_embedding is not None:
            return self._default_embedding
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.default_embedding in self.meta.columns:
                return self.meta[self._adata_ids.default_embedding]
            else:
                return None

    @default_embedding.setter
    def default_embedding(self, x: str):
        x = self._value_protection(attr="default_embedding", allowed=self.ontology_container_sfaira.default_embedding,
                                   attempted=x)
        self._default_embedding = x

    @property
    def development_stage(self) -> Union[None, str]:
        if self._development_stage is not None:
            return self._development_stage
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.development_stage in self.meta.columns:
                return self.meta[self._adata_ids.development_stage]
            else:
                return None

    @development_stage.setter
    def development_stage(self, x: str):
        x = self._value_protection(attr="development_stage",
                                   allowed=self.ontology_container_sfaira.development_stage[self.organism],
                                   attempted=x)
        self._development_stage = x

    @property
    def disease(self) -> Union[None, str]:
        if self._disease is not None:
            return self._disease
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.disease in self.meta.columns:
                return self.meta[self._adata_ids.disease]
            else:
                return None

    @disease.setter
    def disease(self, x: str):
        x = self._value_protection(attr="disease", allowed=self.ontology_container_sfaira.disease,
                                   attempted=x)
        self._disease = x

    @property
    def doi_journal(self) -> str:
        """
        The prepring publication (secondary) DOI associated with the study.
        See also `.doi_journal`.
        """
        return self._doi_journal

    @doi_journal.setter
    def doi_journal(self, x: str):
        self._doi_journal = x

    @property
    def doi_preprint(self) -> str:
        """
        The journal publication (main) DOI associated with the study.
        See also `.doi_preprint`.
        """
        return self._doi_preprint

    @doi_preprint.setter
    def doi_preprint(self, x: str):
        self._doi_preprint = x

    @property
    def doi(self) -> List[str]:
        """
        All publication DOI associated with the study which are the journal publication and the preprint.
        See also `.doi_preprint`, `.doi_journal`.
        """
        dois = []
        if self.doi_journal is not None:
            dois.append(self.doi_journal)
        if self.doi_preprint is not None:
            dois.append(self.doi_preprint)
        return dois

    @property
    def doi_main(self) -> str:
        """
        The main DOI associated with the study which is the journal publication if available, otherwise the preprint.
        See also `.doi_preprint`, `.doi_journal`.
        """
        return self.doi_preprint if self.doi_journal is None else self.doi_journal

    @property
    def directory_formatted_doi(self) -> str:
        return get_directory_formatted_doi(x=self.doi_main)

    @property
    def download_url_data(self) -> Union[Tuple[List[str]], Tuple[List[None]]]:
        """
        Data download website(s).

        Save as tuple with single element, which is a list of all download websites relevant to dataset.
        :return:
        """
        x = self._download_url_data
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        return x

    @download_url_data.setter
    def download_url_data(self, x: Union[str, None, List[str], Tuple[str], List[None], Tuple[None]]):
        # Formats to tuple with single element, which is a list of all download websites relevant to dataset,
        # which can be used as a single element column in a pandas data frame.
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        self._download_url_data = x

    @property
    def download_url_meta(self) -> Union[Tuple[List[str]], Tuple[List[None]]]:
        """
        Meta data download website(s).

        Save as tuple with single element, which is a list of all download websites relevant to dataset.
        :return:
        """
        x = self._download_url_meta
        # if self._download_url_meta is not None:  # TODO add this back in once download_meta is set in all datasets
        #    x = self._download_url_meta
        # else:
        #    if self.meta is None:
        #        self.load_meta(fn=None)
        #    x = self.meta[self._adata_ids.download_url_meta]
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        return x

    @download_url_meta.setter
    def download_url_meta(self, x: Union[str, None, List[str], Tuple[str], List[None], Tuple[None]]):
        # Formats to tuple with single element, which is a list of all download websites relevant to dataset,
        # which can be used as a single element column in a pandas data frame.
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        self._download_url_meta = x

    @property
    def ethnicity(self) -> Union[None, str]:
        if self._ethnicity is not None:
            return self._ethnicity
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.ethnicity in self.meta.columns:
                return self.meta[self._adata_ids.ethnicity]
            else:
                return None

    @ethnicity.setter
    def ethnicity(self, x: str):
        x = self._value_protection(attr="ethnicity", allowed=self.ontology_container_sfaira.ethnicity[self.organism],
                                   attempted=x)
        self._ethnicity = x

    @property
    def id(self) -> str:
        if self._id is not None:
            return self._id
        else:
            raise AttributeError(f"Dataset ID was not set in dataloader in {self.doi_main}, please ensure the "
                                 f"dataloader constructor of this dataset contains a call to self.set_dataset_id()")

    @id.setter
    def id(self, x: str):
        self._id = x

    @property
    def individual(self) -> Union[None, str]:
        if self._individual is not None:
            return self._individual
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.individual in self.meta.columns:
                return self.meta[self._adata_ids.individual]
            else:
                return None

    @individual.setter
    def individual(self, x: str):
        self._individual = x

    @property
    def loaded(self) -> bool:
        """
        :return: Whether DataSet was loaded into memory.
        """
        return self.adata is not None

    @property
    def meta(self) -> Union[None, pd.DataFrame]:
        return self._meta

    @meta.setter
    def meta(self, x: Union[None, pd.DataFrame]):
        # Make sure formatting is correct:
        if x is not None:
            for k, v in x.items():
                v = v.tolist()  # avoid numpy data types
                if k not in META_DATA_FIELDS.keys():
                    raise ValueError(f"did not find {k} in format look up table")
                else:
                    if x[k] is not None:  # None is always allowed.
                        if not isinstance(v[0], META_DATA_FIELDS[k]):
                            raise ValueError(f"key '{k}' of value `{v[0]}` and signature `{type(v[0])}` "
                                             f"in meta data table did not match signature "
                                             f"{str(META_DATA_FIELDS[k])}")
        self._meta = x

    @property
    def ncells(self) -> int:
        # ToDo cache this if it was loaded from meta?
        if self.adata is not None:
            x = self.adata.n_obs
        elif self._ncells is not None:
            x = self._ncells
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            x = self.meta[self._adata_ids.ncells]
        return int(x)

    @property
    def normalization(self) -> Union[None, str]:
        if self._normalization is not None:
            return self._normalization
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.normalization in self.meta.columns:
                return self.meta[self._adata_ids.normalization]
            else:
                return None

    @normalization.setter
    def normalization(self, x: str):
        x = self._value_protection(attr="normalization", allowed=self.ontology_container_sfaira.normalization,
                                   attempted=x)
        self._normalization = x

    @property
    def primary_data(self) -> Union[None, bool]:
        if self._primary_data is not None:
            return self._primary_data
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.primary_data in self.meta.columns:
                return self.meta[self._adata_ids.primary_data]
            else:
                return None

    @primary_data.setter
    def primary_data(self, x: bool):
        x = self._value_protection(attr="primary_data", allowed=self.ontology_container_sfaira.primary_data,
                                   attempted=x)
        self._primary_data = x

    @property
    def organ(self) -> Union[None, str]:
        if self._organ is not None:
            return self._organ
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.organ in self.meta.columns:
                return self.meta[self._adata_ids.organ]
            else:
                return None

    @organ.setter
    def organ(self, x: str):
        x = self._value_protection(attr="organ", allowed=self.ontology_container_sfaira.organ, attempted=x)
        self._organ = x

    @property
    def organism(self) -> Union[None, str]:
        if self._organism is not None:
            return self._organism
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.organism in self.meta.columns:
                return self.meta[self._adata_ids.organism]
            else:
                return None

    @organism.setter
    def organism(self, x: str):
        x = self._value_protection(attr="organism", allowed=self.ontology_container_sfaira.organism, attempted=x)
        # Update ontology container so that correct ontologies are queried:
        self.ontology_container_sfaira.organism_cache = x
        self._organism = x

    @property
    def sample_source(self) -> Union[None, str]:
        if self._sample_source is not None:
            return self._sample_source
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.sample_source in self.meta.columns:
                return self.meta[self._adata_ids.sample_source]
            else:
                return None

    @sample_source.setter
    def sample_source(self, x: str):
        x = self._value_protection(attr="sample_source", allowed=self.ontology_container_sfaira.sample_source,
                                   attempted=x)
        self._sample_source = x

    @property
    def sex(self) -> Union[None, str]:
        if self._sex is not None:
            return self._sex
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.sex in self.meta.columns:
                return self.meta[self._adata_ids.sex]
            else:
                return None

    @sex.setter
    def sex(self, x: str):
        x = self._value_protection(attr="sex", allowed=self.ontology_container_sfaira.sex, attempted=x)
        self._sex = x

    @property
    def source(self) -> str:
        return self._source

    @source.setter
    def source(self, x: Union[str, None]):
        self._source = x

    @property
    def state_exact(self) -> Union[None, str]:
        if self._state_exact is not None:
            return self._state_exact
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.state_exact in self.meta.columns:
                return self.meta[self._adata_ids.state_exact]
            else:
                return None

    @state_exact.setter
    def state_exact(self, x: str):
        self._state_exact = x

    @property
    def tech_sample(self) -> Union[None, str]:
        if self._tech_sample is not None:
            return self._tech_sample
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.tech_sample in self.meta.columns:
                return self.meta[self._adata_ids.tech_sample]
            else:
                return None

    @tech_sample.setter
    def tech_sample(self, x: str):
        self._tech_sample = x

    @property
    def year(self) -> Union[None, int]:
        if self._year is not None:
            return self._year
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.year in self.meta.columns:
                return self.meta[self._adata_ids.year]
            else:
                return None

    @year.setter
    def year(self, x: int):
        x = self._value_protection(attr="year", allowed=self.ontology_container_sfaira.year, attempted=x)
        self._year = x

    @property
    def celltypes_universe(self):
        if self._celltype_universe is None:
            self._celltype_universe = CelltypeUniverse(
                cl=getattr(self.ontology_container_sfaira, "cell_type"),
                uberon=getattr(self.ontology_container_sfaira, "organ"),
            )
        return self._celltype_universe

    @property
    def cell_type_map(self) -> dict:
        return self._ontology_class_map

    @cell_type_map.setter
    def cell_type_map(self, x: pd.DataFrame):
        assert x.shape[1] in [2, 3], f"{x.shape} in {self.id}"
        assert x.columns[0] == self._adata_ids.classmap_source_key
        assert x.columns[1] == self._adata_ids.classmap_target_key
        # Check for weird entries:
        # nan arises if columns was empty in that row.
        nan_vals = np.where([
            False if isinstance(x, str) else (np.isnan(x) or x is None)
            for x in x[self._adata_ids.classmap_target_key].values.tolist()
        ])[0]
        assert len(nan_vals) == 0, \
            f"Found nan target values in {self.id} for {x[self._adata_ids.classmap_target_key].values[nan_vals]}"
        # Transform data frame into a mapping dictionary:
        self._ontology_class_map = dict(list(zip(
            x[self._adata_ids.classmap_source_key].values.tolist(),
            x[self._adata_ids.classmap_target_key].values.tolist()
        )))

    def __crossref_query(self, k):
        """
        Queries cross REST API via package crossref_commons.

        :param k: Key to extract from crossref query container.
        :return:
        """
        from crossref_commons.retrieval import get_entity
        from crossref_commons.types import EntityType, OutputType
        try:
            attempt_counter = 0
            while True:
                x = None
                try:
                    attempt_counter += 1
                    x = get_entity(self.doi_main, EntityType.PUBLICATION, OutputType.JSON)[k]
                except ConnectionError as e:
                    # Terminate trial after 5 attempts with ConnectionError:
                    if attempt_counter > 5:
                        raise ConnectionError(e)
                finally:
                    if k == "author":
                        pass
                    return x
        except ValueError as e:
            print(f"ValueError: {e}")
            return None
        except ConnectionError as e:
            print(f"ConnectionError: {e}")
            return None

    @property
    def author(self) -> str:
        if self._author is not None:
            return self._author
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is None or self._adata_ids.author not in self.meta.columns:
                raise ValueError("author must be set but was neither set in constructor nor in meta data")
            return self.meta[self._adata_ids.author]

    @author.setter
    def author(self, x: str):
        self._author = x

    @property
    def title(self):
        if self._title is not None:
            return self._title
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids.title in self.meta.columns:
                return self.meta[self._adata_ids.title]
            else:
                return self.__crossref_query(k="title")

    def _value_protection(
            self,
            attr: str,
            allowed: Union[Ontology, None],
            attempted
    ):
        """
        Check whether value is from set of allowed values.

        Does not check if allowed is None.
        Cleans entry to term name if ontology ID is provided.

        :param attr: Attribute to set.
        :param allowed: Constraint for values of `attr`.
            Either ontology instance used to constrain entries, or list of allowed values.
        :param attempted: Value(s) to attempt to set in `attr`.
        :return:
        """
        if not isinstance(attempted, list):
            if isinstance(attempted, np.ndarray):
                attempted_ls = attempted.tolist()
            elif isinstance(attempted, tuple):
                attempted_ls = list(attempted)
            else:
                attempted_ls = [attempted]
        else:
            attempted_ls = attempted
        attempted_clean = []
        for x in attempted_ls:
            if allowed is None:
                attempted_clean.append(x)
            elif isinstance(allowed, Ontology):
                if attr == "disease" and (x.lower() == "normal" or x.lower() == "healthy"):
                    # TODO required because of missing streamlining between sfaira and 10x, remove in future.
                    attempted_clean.append("healthy")
                elif x in allowed.node_names:
                    attempted_clean.append(x)
                else:
                    if isinstance(allowed, OntologyHierarchical) and x in allowed.node_ids:
                        attempted_clean.append(allowed.convert_to_name(x))
                    else:
                        raise ValueError(f"'{x}' is not a valid entry for {attr} in data set {self.doi}.")
            else:
                raise ValueError(f"argument allowed of type {type(allowed)} is not a valid entry for {attr}.")
        # Flatten attempts if only one was made:
        if len(attempted_clean) == 1:
            attempted_clean = attempted_clean[0]
        return attempted_clean

    def subset_cells(self, key, values):
        """
        Subset list of adata objects based on cell-wise properties.

        These keys are properties that are not available in lazy model and require loading first because the
        subsetting works on the cell-level: .adata are maintained but reduced to matches.

        :param key: Property to subset by. Options:

            - "assay_sc" points to self.assay_sc_obs_key
            - "assay_differentiation" points to self.assay_differentiation_obs_key
            - "assay_type_differentiation" points to self.assay_type_differentiation_obs_key
            - "cell_line" points to self.cell_line
            - "cell_type" points to self.cell_type_obs_key
            - "developmental_stage" points to self.developmental_stage_obs_key
            - "ethnicity" points to self.ethnicity_obs_key
            - "organ" points to self.organ_obs_key
            - "organism" points to self.organism_obs_key
            - "sample_source" points to self.sample_source_obs_key
            - "sex" points to self.sex_obs_key
            - "state_exact" points to self.state_exact_obs_key
        :param values: Classes to overlap to.
        :return:
        """
        if not isinstance(values, list):
            values = [values]

        def get_subset_idx(samplewise_key, cellwise_key):
            try:
                sample_attr = getattr(self, samplewise_key)
                if not isinstance(sample_attr, list):
                    sample_attr = [sample_attr]
            except AttributeError:
                sample_attr = None
            obs_key = getattr(self, cellwise_key)
            if sample_attr is not None and len(sample_attr) == 1:
                # Only use sample-wise subsetting if the sample-wise attribute is unique (not mixed).
                if np.any([x in values for x in sample_attr]):
                    idx = np.arange(1, self.ncells)
                else:
                    idx = np.array([])
            elif obs_key is not None:
                assert self.adata is not None, "adata was not yet loaded"
                values_found = self.adata.obs[obs_key].values
                values_found_unique = np.unique(values_found)
                try:
                    ontology = getattr(self.ontology_container_sfaira, samplewise_key)
                except AttributeError:
                    raise ValueError(f"{key} not a valid property of ontology_container object")
                # Test only unique elements  found in ontology to save time.
                values_found_unique_matched = [
                    x for x in values_found_unique if np.any([
                        is_child(query=x, ontology=ontology, ontology_parent=y)
                        for y in values
                    ])
                ]
                # TODO keep this logging for now to catch undesired behaviour resulting from loaded edges in ontologies.
                print(f"matched cell-wise keys {str(values_found_unique_matched)} in data set {self.id}")
                idx = np.where([x in values_found_unique_matched for x in values_found])[0]
            else:
                assert False, "no subset chosen"
            return idx

        idx_keep = get_subset_idx(samplewise_key=key, cellwise_key=key + "_obs_key")
        self.adata = self.adata[idx_keep, :].copy()  # if len(idx_keep) > 0 else None

    def show_summary(self):
        print(f"{(self.supplier, self.organism, self.organ, self.assay_sc, self.disease)}")

    def __assert_loaded(self):
        if self.adata is None:
            raise ValueError("adata was not loaded, this is necessary for this operation")
