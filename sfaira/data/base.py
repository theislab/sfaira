from __future__ import annotations

import abc
import anndata
import h5py
import multiprocessing
import numpy as np
import pandas as pd
import os
from os import PathLike
import pandas
import pydoc
import scipy.sparse
from typing import Dict, List, Tuple, Union
import warnings
import urllib.request
import urllib.parse
import urllib.error
import cgi
import ssl

from sfaira.versions.genome_versions import SuperGenomeContainer
from sfaira.versions.metadata import Ontology, CelltypeUniverse
from sfaira.consts import AdataIdsExtended, AdataIdsSfaira, META_DATA_FIELDS, OCS
from sfaira.data.utils import read_yaml

UNS_STRING_META_IN_OBS = "__obs__"


def map_fn(inputs):
    """
    Functional to load data set with predefined additional actions.

    :param inputs:
    :return: None if function ran, error report otherwise
    """
    ds, remove_gene_version, match_to_reference, load_raw, allow_caching, func, kwargs_func = inputs
    try:
        ds.load(
            remove_gene_version=remove_gene_version,
            match_to_reference=match_to_reference,
            load_raw=load_raw,
            allow_caching=allow_caching,
        )
        if func is not None:
            x = func(ds, **kwargs_func)
            ds.clear()
            return x
        else:
            return None
    except FileNotFoundError as e:
        return ds.id, e,


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

    _age: Union[None, str]
    _assay: Union[None, str]
    _author: Union[None, str]
    _bio_sample: Union[None, str]
    _development_stage: Union[None, str]
    _doi: Union[None, str]
    _download_url_data: Union[Tuple[List[None]], Tuple[List[str]], None]
    _download_url_meta: Union[Tuple[List[None]], Tuple[List[str]], None]
    _ethnicity: Union[None, str]
    _healthy: Union[None, bool]
    _id: Union[None, str]
    _individual: Union[None, str]
    _ncells: Union[None, int]
    _normalization: Union[None, str]
    _organ: Union[None, str]
    _organism: Union[None, str]
    _sex: Union[None, str]
    _source: Union[None, str]
    _state_exact: Union[None, str]
    _bio_sample: Union[None, str]
    _year: Union[None, int]

    _age_obs_key: Union[None, str]
    _assay_obs_key: Union[None, str]
    _cellontology_id_obs_key: Union[None, str]
    _cellontology_original_obs_key: Union[None, str]
    _development_stage_obs_key: Union[None, str]
    _ethnicity_obs_key: Union[None, str]
    _healthy_obs_key: Union[None, str]
    _healthy_obs_key: Union[None, str]
    _individual: Union[None, str]
    _organ_obs_key: Union[None, str]
    _organism_obs_key: Union[None, str]
    _bio_sample_obs_key: Union[None, str]
    _sex_obs_key: Union[None, str]
    _state_exact_obs_key: Union[None, str]
    _tech_sample_obs_key: Union[None, str]

    _healthy_state_healthy: Union[None, str]

    _var_symbol_col: Union[None, str]
    _var_ensembl_col: Union[None, str]

    _celltype_universe: Union[None, CelltypeUniverse]
    _ontology_class_map: Union[None, dict]

    sample_fn: Union[None, str]
    _sample_fns: Union[None, List[str]]

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            load_func=None,
            yaml_path: Union[str, None] = None,
            sample_fn: Union[str, None] = None,
            sample_fns: Union[List[str], None] = None,
            **kwargs
    ):
        self._adata_ids_sfaira = AdataIdsSfaira()
        self._ontology_container_sfaira = OCS  # Using a pre-instantiated version of this yields drastic speed-ups.

        self.adata = None
        self.meta = None
        self.genome = None
        self.data_dir_base = data_path
        self.meta_path = meta_path
        self.cache_path = cache_path

        self._age = None
        self._author = None
        self._bio_sample = None
        self._development_stage = None
        self._doi = None
        self._download_url_data = None
        self._download_url_meta = None
        self._ethnicity = None
        self._healthy = None
        self._id = None
        self._individual = None
        self._ncells = None
        self._normalization = None
        self._organ = None
        self._organism = None
        self._assay = None
        self._sex = None
        self._source = None
        self._state_exact = None
        self._tech_sample = None
        self._year = None

        self._age_obs_key = None
        self._cellontology_id_obs_key = None
        self._cellontology_original_obs_key = None
        self._development_stage_obs_key = None
        self._ethnicity_obs_key = None
        self._healthy_obs_key = None
        self._individual_obs_key = None
        self._organ_obs_key = None
        self._organism_obs_key = None
        self._assay_obs_key = None
        self._bio_sample_obs_key = None
        self._sex_obs_key = None
        self._state_exact_obs_key = None
        self._tech_sample_obs_key = None

        self._healthy_state_healthy = None

        self._var_symbol_col = None
        self._var_ensembl_col = None

        self.class_maps = {"0": {}}
        self._unknown_celltype_identifiers = self._adata_ids_sfaira.unknown_celltype_identifier

        self._celltype_universe = None
        self._ontology_class_map = None

        self.sample_fn = sample_fn
        self._sample_fns = sample_fns

        # Check if YAML files exists, read meta data from there if available:
        if yaml_path is not None:
            assert os.path.exists(yaml_path), f"did not find yaml {yaml_path}"
            yaml_vals = read_yaml(fn=yaml_path)
            for k, v in yaml_vals["attr"].items():
                if v is not None and k not in ["sample_fns", "sample_ids", "dataset_index"]:
                    if isinstance(v, dict):  # v is a dictionary over file-wise meta-data items
                        assert self.sample_fn in v.keys(), f"did not find key {self.sample_fn} in yamls keys for {k}"
                        setattr(self, k, v[self.sample_fn])
                    else:  # v is a meta-data item
                        setattr(self, k, v)
            # ID can be set now already because YAML was used as input instead of child class constructor.
            self.set_dataset_id(idx=yaml_vals["meta"]["dataset_index"])

        self.load_func = load_func

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

        urls = self.download_url_data[0][0] + self.download_url_meta[0][0]

        for url in urls:
            if url is None:
                continue
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

            elif url.split(",")[0].startswith('syn'):
                fn = ",".join(url.split(",")[1:])
                if os.path.isfile(os.path.join(self.data_dir, fn)):
                    print(f"File {fn} already found on disk, skipping download.")
                else:
                    self._download_synapse(url.split(",")[0], fn, **kwargs)

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
                if os.path.isfile(os.path.join(self.data_dir, fn)):
                    print(f"File {fn} already found on disk, skipping download.")
                else:
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
            warnings.warn("Caching enabled, but Dataset.id or Dataset.doi not set. Disabling caching for now.")
            return None
        else:
            if self.cache_path is None:
                cache = self.data_dir
            else:
                cache = os.path.join(self.cache_path, self.directory_formatted_doi)
            return os.path.join(cache, "cache", self._directory_formatted_id + ".h5ad")

    def _load_cached(
            self,
            load_raw: bool,
            allow_caching: bool,
    ):
        """
        Wraps data set specific load and allows for caching.

        Cache is written into director named after doi and h5ad named after data set id.

        :param load_raw: Loads unprocessed version of data if available in data loader.
        :param allow_caching: Whether to allow method to cache adata object for faster re-loading.
        :return:
        """

        def _cached_reading(filename):
            if filename is not None:
                if os.path.exists(filename):
                    self.adata = anndata.read_h5ad(filename)
                else:
                    warnings.warn(f"Cached loading enabled, but cache file {filename} not found. "
                                  f"Loading from raw files.")
                    self.adata = self.load_func(self.data_dir, self.sample_fn)
            else:
                self.adata = self.load_func(self.data_dir, self.sample_fn)

        def _cached_writing(filename):
            if filename is not None:
                dir_cache = os.path.dirname(filename)
                if not os.path.exists(dir_cache):
                    os.makedirs(dir_cache)
                self.adata.write_h5ad(filename)

        if load_raw and allow_caching:
            self.adata = self.load_func(self.data_dir, self.sample_fn)
            _cached_writing(self.cache_fn)
        elif load_raw and not allow_caching:
            self.adata = self.load_func(self.data_dir, self.sample_fn)
        elif not load_raw and allow_caching:
            _cached_reading(self.cache_fn)
            _cached_writing(self.cache_fn)
        else:  # not load_raw and not allow_caching
            _cached_reading(self.cache_fn)

    def load(
            self,
            remove_gene_version: bool = True,
            match_to_reference: Union[str, bool, None] = None,
            load_raw: bool = False,
            allow_caching: bool = True,
    ):
        if match_to_reference and not remove_gene_version:
            warnings.warn("it is not recommended to enable matching the feature space to a genomes reference"
                          "while not removing gene versions. this can lead to very poor matching results")

        # Set default genomes per organism if none provided:
        if isinstance(match_to_reference, str):
            genome = match_to_reference
        elif match_to_reference is None or (isinstance(match_to_reference, bool) and match_to_reference):
            if self.organism == "human":
                genome = "Homo_sapiens_GRCh38_97"
                warnings.warn(f"using default genome {genome}")
            elif self.organism == "mouse":
                genome = "Mus_musculus_GRCm38_97"
                warnings.warn(f"using default genome {genome}")
            else:
                raise ValueError(f"genome was not supplied and no default genome found for organism {self.organism}")
        elif not match_to_reference:
            genome = None
        else:
            raise ValueError(f"invalid choice for match_to_reference={match_to_reference}")
        self._set_genome(genome=genome)

        # Set path to dataset directory
        if self.data_dir is None:
            raise ValueError("No sfaira data repo path provided in constructor.")

        # Run data set-specific loading script:
        self._load_cached(load_raw=load_raw, allow_caching=allow_caching)
        # Set data-specific meta data in .adata:
        self._set_metadata_in_adata(adata_ids=self._adata_ids_sfaira)
        # Set loading hyper-parameter-specific meta data:
        self.adata.uns[self._adata_ids_sfaira.load_raw] = load_raw
        self.adata.uns[self._adata_ids_sfaira.mapped_features] = match_to_reference
        self.adata.uns[self._adata_ids_sfaira.remove_gene_version] = remove_gene_version
        # Streamline feature space:
        self._convert_and_set_var_names(match_to_reference=match_to_reference)
        self._collapse_gene_versions(remove_gene_version=remove_gene_version)
        if match_to_reference:
            self._match_features_to_reference()

    load.__doc__ = load_doc

    def _convert_and_set_var_names(
            self,
            match_to_reference: Union[str, bool, None],
            symbol_col: str = None,
            ensembl_col: str = None,
    ):
        # Use defaults defined in data loader if none given to this function.
        if symbol_col is None:
            symbol_col = self.var_symbol_col
        if ensembl_col is None:
            ensembl_col = self.var_ensembl_col
        if not ensembl_col and not symbol_col:
            raise ValueError('Please provide the name of at least the name of the var column containing ensembl ids or'
                             'the name of the var column containing gene symbols')
        # Process given gene names: Full gene names ("symbol") or ENSEMBL IDs ("ensembl").
        # Below the .var column that contain the target IDs are renamed to follow streamlined naming.
        # If the IDs were contained in the index, a new column is added to .var.
        if symbol_col:
            if symbol_col == 'index':
                self.adata.var[self._adata_ids_sfaira.gene_id_names] = self.adata.var.index.values.tolist()
            else:
                assert symbol_col in self.adata.var.columns, f"symbol_col {symbol_col} not found in .var"
                self.adata.var = self.adata.var.rename(
                    {symbol_col: self._adata_ids_sfaira.gene_id_names},
                    axis='columns'
                )
        if ensembl_col:
            if ensembl_col == 'index':
                self.adata.var[self._adata_ids_sfaira.gene_id_ensembl] = self.adata.var.index.values.tolist()
            else:
                assert ensembl_col in self.adata.var.columns, f"ensembl_col {ensembl_col} not found in .var"
                self.adata.var = self.adata.var.rename(
                    {ensembl_col: self._adata_ids_sfaira.gene_id_ensembl},
                    axis='columns'
                )
        # If only symbol or ensembl was supplied, the other one is inferred from a genome mapping dictionary.
        if not ensembl_col and match_to_reference:
            id_dict = self.genome_container.names_to_id_dict
            id_strip_dict = self.genome_container.strippednames_to_id_dict
            # Matching gene names to ensembl ids in the following way: if the gene is present in the ensembl dictionary,
            # match it straight away, if it is not in there we try to match everything in front of the first period in
            # the gene name with a dictionary that was modified in the same way, if there is still no match we append na
            ensids = []
            for n in self.adata.var[self._adata_ids_sfaira.gene_id_names]:
                if n in id_dict.keys():
                    ensids.append(id_dict[n])
                elif n.split(".")[0] in id_strip_dict.keys():
                    ensids.append(id_strip_dict[n.split(".")[0]])
                else:
                    ensids.append('n/a')
            self.adata.var[self._adata_ids_sfaira.gene_id_ensembl] = ensids

        if not symbol_col and match_to_reference:
            id_dict = self.genome_container.id_to_names_dict
            self.adata.var[self._adata_ids_sfaira.gene_id_names] = [
                id_dict[n.split(".")[0]] if n.split(".")[0] in id_dict.keys() else 'n/a'
                for n in self.adata.var[self._adata_ids_sfaira.gene_id_ensembl]
            ]

        if match_to_reference:
            # Lastly, the index of .var is set to ensembl IDs.
            try:  # debugging
                self.adata.var.index = self.adata.var[self._adata_ids_sfaira.gene_id_index].values.tolist()
            except KeyError as e:
                raise KeyError(e)
            self.adata.var_names_make_unique()

    def _collapse_gene_versions(self, remove_gene_version):
        """
        Remove version tag on ensembl gene ID so that different versions are superimposed downstream.

        :param remove_gene_version:
        :return:
        """
        if remove_gene_version:
            new_index = [x.split(".")[0] for x in self.adata.var_names.tolist()]
            # Collapse if necessary:
            new_index_collapsed = list(np.unique(new_index))
            if len(new_index_collapsed) < self.adata.n_vars:
                print("WARNING: duplicate features detected after removing gene versions. "
                      "the code to collapse these features is implemented but not tested.")
                idx_map = np.array([new_index_collapsed.index(x) for x in new_index])
                # Need reverse sorting to find index of last element in sorted list to split array using list index().
                idx_map_sorted_fwd = np.argsort(idx_map)
                idx_map_sorted_rev = idx_map_sorted_fwd[::-1].tolist()
                n_genes = len(idx_map_sorted_rev)
                # 1. Sort array in non-reversed order: idx_map_sorted_rev[::-1]
                # 2. Split into chunks based on blocks of identical entries in idx_map, using the occurrence of the
                # last element of each block as block boundaries:
                # n_genes - 1 - idx_map_sorted_rev.index(x)
                # Note that the blocks are named as positive integers starting at 1, without gaps.
                counts = np.concatenate([np.sum(x, axis=1, keepdims=True)
                                         for x in np.split(self.adata[:, idx_map_sorted_fwd].X,  # forward ordered data
                                                           indices_or_sections=[
                                                               n_genes - 1 - idx_map_sorted_rev.index(x)  # last occurrence of element in forward order
                                                               for x in np.arange(0, len(new_index_collapsed) - 1)],  # -1: do not need end of last partition
                                                           axis=1
                                                           )
                                         ][::-1], axis=1)
                # Remove varm and populate var with first occurrence only:
                obs_names = self.adata.obs_names
                self.adata = anndata.AnnData(
                    X=counts,
                    obs=self.adata.obs,
                    obsm=self.adata.obsm,
                    var=self.adata.var.iloc[[new_index.index(x) for x in new_index_collapsed]],
                    uns=self.adata.uns
                )
                self.adata.obs_names = obs_names
                self.adata.var_names = new_index_collapsed
                new_index = new_index_collapsed
            self.adata.var[self._adata_ids_sfaira.gene_id_ensembl] = new_index
            self.adata.var.index = self.adata.var[self._adata_ids_sfaira.gene_id_ensembl].values

    def _match_features_to_reference(self):
        """
        Match feature space to a genomes provided with sfaira
        """
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
        data_ids = self.adata.var[self._adata_ids_sfaira.gene_id_ensembl].values
        idx_feature_kept = np.where([x in self.genome_container.ensembl for x in data_ids])[0]
        idx_feature_map = np.array([self.genome_container.ensembl.index(x)
                                    for x in data_ids[idx_feature_kept]])
        # Remove unmapped genes
        x = x[:, idx_feature_kept]

        # Create reordered feature matrix based on reference and convert to csr
        x_new = scipy.sparse.csc_matrix((x.shape[0], self.genome_container.ngenes), dtype=x.dtype)
        # copying this over to the new matrix in chunks of size `steps` prevents a strange scipy error:
        # ... scipy/sparse/compressed.py", line 922, in _zero_many i, j, offsets)
        # ValueError: could not convert integer scalar
        step = 2000
        if step < len(idx_feature_map):
            i = 0
            for i in range(0, len(idx_feature_map), step):
                x_new[:, idx_feature_map[i:i + step]] = x[:, i:i + step]
            x_new[:, idx_feature_map[i + step:]] = x[:, i + step:]
        else:
            x_new[:, idx_feature_map] = x

        x_new = x_new.tocsr()

        self.adata = anndata.AnnData(
            X=x_new,
            obs=self.adata.obs,
            obsm=self.adata.obsm,
            var=pd.DataFrame(data={'names': self.genome_container.names,
                                   self._adata_ids_sfaira.gene_id_ensembl: self.genome_container.ensembl},
                             index=self.genome_container.ensembl),
            uns=self.adata.uns
        )

    def _set_metadata_in_adata(self, adata_ids: AdataIdsExtended):
        """
        Copy meta data from dataset class in .anndata.

        :return:
        """
        # Set data set-wide attributes (.uns):
        self.adata.uns[adata_ids.annotated] = self.annotated
        self.adata.uns[adata_ids.author] = self.author
        self.adata.uns[adata_ids.doi] = self.doi
        self.adata.uns[adata_ids.download_url_data] = self.download_url_data
        self.adata.uns[adata_ids.download_url_meta] = self.download_url_meta
        self.adata.uns[adata_ids.id] = self.id
        self.adata.uns[adata_ids.normalization] = self.normalization
        self.adata.uns[adata_ids.year] = self.year

        # Set cell-wise or data set-wide attributes (.uns / .obs):
        # These are saved in .uns if they are data set wide to save memory.
        for x, y, z, v in (
            [self.age, adata_ids.age, self.age_obs_key,
             self._ontology_container_sfaira.ontology_age],
            [self.assay, adata_ids.assay, self.assay_obs_key,
             self._ontology_container_sfaira.ontology_protocol],
            [self.bio_sample, adata_ids.bio_sample, self.bio_sample_obs_key, None],
            [self.development_stage, adata_ids.development_stage, self.development_stage_obs_key,
             self._ontology_container_sfaira.ontology_dev_stage],
            [self.ethnicity, adata_ids.ethnicity, self.ethnicity_obs_key,
             self._ontology_container_sfaira.ontology_ethnicity],
            [self.healthy, adata_ids.healthy, self.healthy_obs_key,
             self._ontology_container_sfaira.ontology_healthy],
            [self.individual, adata_ids.individual, self.individual_obs_key, None],
            [self.organ, adata_ids.organ, self.organ_obs_key,
             self._ontology_container_sfaira.ontology_organism],
            [self.organism, adata_ids.organism, self.organism_obs_key,
             self._ontology_container_sfaira.ontology_organism],
            [self.sex, adata_ids.sex, self.sex_obs_key, self._ontology_container_sfaira.ontology_sex],
            [self.state_exact, adata_ids.state_exact, self.state_exact_obs_key, None],
            [self.tech_sample, adata_ids.tech_sample, self.tech_sample_obs_key, None],
        ):
            if x is None and z is None:
                self.adata.uns[y] = None
            elif x is not None and z is None:
                # Attribute supplied per data set: Write into .uns.
                self.adata.uns[y] = x
            elif z is not None:
                # Attribute supplied per cell: Write into .obs.
                # Search for direct match of the sought-after column name or for attribute specific obs key.
                if z not in self.adata.obs.keys():
                    # This should not occur in single data set loaders (see warning below) but can occur in
                    # streamlined data loaders if not all instances of the streamlined data sets have all columns
                    # in .obs set.
                    self.adata.uns[y] = None
                    print(f"WARNING: attribute {y} of data set {self.id} was not found in column {z}")  # debugging
                else:
                    # Include flag in .uns that this attribute is in .obs:
                    self.adata.uns[y] = UNS_STRING_META_IN_OBS
                    # Remove potential pd.Categorical formatting:
                    self._value_protection(
                        attr="obs", allowed=v, attempted=np.unique(self.adata.obs[z].values).tolist())
                    self.adata.obs[y] = self.adata.obs[z].values.tolist()
            else:
                assert False, "switch option should not occur"
        # Set cell-wise attributes (.obs):
        # None so far other than celltypes.
        # Set cell types:
        # Map cell type names from raw IDs to ontology maintained ones::
        if self.cellontology_original_obs_key is not None:
            self.project_celltypes_to_ontology()

    def streamline(self, format: str = "sfaira", clean: bool = False):
        """
        Streamline the adata instance to output format.

        Output format are saved in ADATA_FIELDS* classes.

        :param format: Export format.

            - "sfaira"
            - "cellxgene"
        :param clean: Whether to delete non-streamlined fields.
        :return:
        """
        if format == "sfaira":
            adata_fields = self._adata_ids_sfaira
        elif format == "sfaira":
            from sfaira.consts import AdataIdsCellxgene
            adata_fields = AdataIdsCellxgene()
        self._set_metadata_in_adata(adata_ids=adata_fields)
        if clean:
            if self.adata.varm is not None:
                del self.adata.varm
            if self.adata.obsm is not None:
                del self.adata.obsm
            if self.adata.varm is not None:
                del self.adata.varp
            if self.adata.obsp is not None:
                del self.adata.obsp
            # Only retain target elements in adata.uns:
            self.adata.obs = self.adata.uns[[
                adata_fields.annotated,
                adata_fields.author,
                adata_fields.doi,
                adata_fields.download_url_data,
                adata_fields.download_url_meta,
                adata_fields.id,
                adata_fields.normalization,
                adata_fields.year,
            ]]
            # Only retain target elements in adata.var:
            self.adata.obs = self.adata.var[[
                adata_fields.gene_id_names,
                adata_fields.gene_id_ensembl,
            ]]
            # Only retain target columns in adata.obs:
            self.adata.obs = self.adata.obs[[
                adata_fields.age,
                adata_fields.bio_sample,
                adata_fields.development_stage,
                adata_fields.ethnicity,
                adata_fields.healthy,
                adata_fields.individual,
                adata_fields.organ,
                adata_fields.organism,
                adata_fields.sex,
                adata_fields.state_exact,
                adata_fields.tech_sample,
            ]]

    def load_tobacked(
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
        self.load(
            remove_gene_version=True,
            match_to_reference=genome,
            load_raw=load_raw,
            allow_caching=allow_caching
        )
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
                if k == self._adata_ids_sfaira.dataset:
                    adata_backed.obs.loc[np.sort(idx), self._adata_ids_sfaira.dataset] = [
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
                    (k, [self.id for i in range(len(idx))]) if k == self._adata_ids_sfaira.dataset
                    else (k, self.adata.obs[k].values[np.argsort(idx)]) if k in self.adata.obs.columns
                    else (k, [self.adata.uns[k] for _ in range(len(idx))]) if k in list(self.adata.uns.keys())
                    else (k, ["key_not_found" for _ in range(len(idx))])
                    for k in adata_backed.obs.columns
                ]))
            )
            self.clear()
        else:
            raise ValueError(f"Did not recognize backed AnnData.X format {type(adata_backed.X)}")

    def _set_genome(self, genome: Union[str, None]):
        if genome is not None:
            if genome.lower().startswith("homo_sapiens"):
                g = SuperGenomeContainer(
                    organism="human",
                    genome=genome
                )
            elif genome.lower().startswith("mus_musculus"):
                g = SuperGenomeContainer(
                    organism="mouse",
                    genome=genome
                )
            else:
                raise ValueError(f"Genome {genome} not recognised. Needs to start with 'Mus_Musculus' or "
                                 f"'Homo_Sapiens'.")
        else:
            g = None
        self.genome_container = g

    @property
    def doi_cleaned_id(self):
        return "_".join(self.id.split("_")[:-1])

    @property
    def fn_ontology_class_map_tsv(self):
        """Standardised file name under which cell type conversion tables are saved."""
        return self.doi_cleaned_id + ".tsv"

    def write_ontology_class_map(
            self,
            fn,
            protected_writing: bool = True,
            **kwargs
    ):
        """
        Load class maps of free text cell types to ontology classes.

        :param fn: File name of csv to load class maps from.
        :param protected_writing: Only write if file was not already found.
        :return:
        """
        if not self.annotated:
            warnings.warn(f"attempted to write ontology classmaps for data set {self.id} without annotation")
        else:
            labels_original = np.sort(np.unique(self.adata.obs[self._adata_ids_sfaira.cell_types_original].values))
            tab = self.celltypes_universe.prepare_celltype_map_tab(
                source=labels_original,
                match_only=False,
                anatomical_constraint=self.organ,
                include_synonyms=True,
                omit_list=self._unknown_celltype_identifiers,
                **kwargs
            )
            if not os.path.exists(fn) or not protected_writing:
                self._write_class_map(fn=fn, tab=tab)

    def _write_class_map(self, fn, tab):
        """
        Write class map.

        :param fn: File name of csv to write class maps to.
        :param tab: Table to write
        :return:
        """
        tab.to_csv(fn, index=False, sep="\t")

    def _read_class_map(self, fn) -> pd.DataFrame:
        """
        Read class map.

        :param fn: File name of csv to load class maps from.
        :return:
        """
        try:
            tab = pd.read_csv(fn, header=0, index_col=None, sep="\t")
        except pandas.errors.ParserError as e:
            print(f"{self.id}")
            raise pandas.errors.ParserError(e)
        return tab

    def load_ontology_class_map(self, fn):
        """
        Load class maps of free text cell types to ontology classes.

        :param fn: File name of csv to load class maps from.
        :return:
        """
        if os.path.exists(fn):
            self.cell_ontology_map = self._read_class_map(fn=fn)
        else:
            warnings.warn(f"file {fn} does not exist")

    def project_celltypes_to_ontology(self):
        """
        Project free text cell type names to ontology based on mapping table.

        ToDo: add ontology ID setting here.

        :return:
        """
        labels_original = self.adata.obs[self.obs_key_cellontology_original].values
        if self.cell_ontology_map is not None:  # only if this was defined
            labels_mapped = [
                self.cell_ontology_map[x] if x in self.cell_ontology_map.keys()
                else x for x in labels_original
            ]
        else:
            labels_mapped = labels_original
        # Validate mapped IDs based on ontology:
        # This aborts with a readable error if there was a target in the mapping file that does not match the
        # ontology.
        self._value_protection(
            attr="celltypes",
            allowed=self.ontology_celltypes,
            attempted=np.unique(labels_mapped).tolist()
        )
        self.adata.obs[self._adata_ids_sfaira.cell_ontology_class] = labels_mapped
        self.adata.obs[self._adata_ids_sfaira.cell_types_original] = labels_original
        # Add cell type IDs into object:
        # The IDs are not read from a source file but inferred based on the class name.
        # TODO this could be changed in the future, this allows this function to be used both on cell type name mapping
        #  files with and without the ID in the third column.
        ids_mapped = [
            self._ontology_container_sfaira.ontology_cell_types.id_from_name(x)
            if x not in [
                self._adata_ids_sfaira.unknown_celltype_identifier,
                self._adata_ids_sfaira.not_a_cell_celltype_identifier
            ] else x
            for x in labels_mapped
        ]
        self.adata.obs[self._adata_ids_sfaira.cell_ontology_id] = ids_mapped

    @property
    def citation(self):
        """
        Return all information necessary to cite data set.

        :return:
        """
        return [self.author, self.year, self.doi]

    # Meta data handling code: Reading, writing and selected properties. Properties are either set in constructor
    # (and saved in self._somename) or accessed in self.meta.

    @property
    def meta_fn(self):
        if self.meta_path is None:
            meta = self.data_dir
        else:
            meta = os.path.join(self.meta_path, self.directory_formatted_doi)

        return os.path.join(meta, "meta", self.doi_cleaned_id + "_meta.csv")

    def load_meta(self, fn: Union[PathLike, str, None]):
        if fn is None:
            if self.meta_fn is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = self.meta_fn
        else:
            if isinstance(fn, str):
                fn = os.path.normpath(fn)
        # Only load meta data if file exists:
        if os.path.isfile(fn):
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
            self.load(
                remove_gene_version=False,
                match_to_reference=None,
                load_raw=True,
                allow_caching=False,
            )
        # Add data-set wise meta data into table:
        meta = pandas.DataFrame({
            self._adata_ids_sfaira.annotated: self.adata.uns[self._adata_ids_sfaira.annotated],
            self._adata_ids_sfaira.author: self.adata.uns[self._adata_ids_sfaira.author],
            self._adata_ids_sfaira.doi: self.adata.uns[self._adata_ids_sfaira.doi],
            self._adata_ids_sfaira.download_url_data: self.adata.uns[self._adata_ids_sfaira.download_url_data],
            self._adata_ids_sfaira.download_url_meta: self.adata.uns[self._adata_ids_sfaira.download_url_meta],
            self._adata_ids_sfaira.id: self.adata.uns[self._adata_ids_sfaira.id],
            self._adata_ids_sfaira.ncells: self.adata.n_obs,
            self._adata_ids_sfaira.normalization: self.adata.uns[self._adata_ids_sfaira.normalization],
            self._adata_ids_sfaira.year: self.adata.uns[self._adata_ids_sfaira.year],
        }, index=range(1))
        # Expand table by variably cell-wise or data set-wise meta data:
        for x in [
            self._adata_ids_sfaira.age,
            self._adata_ids_sfaira.assay,
            self._adata_ids_sfaira.bio_sample,
            self._adata_ids_sfaira.development_stage,
            self._adata_ids_sfaira.ethnicity,
            self._adata_ids_sfaira.healthy,
            self._adata_ids_sfaira.individual,
            self._adata_ids_sfaira.organ,
            self._adata_ids_sfaira.organism,
            self._adata_ids_sfaira.sex,
            self._adata_ids_sfaira.state_exact,
            self._adata_ids_sfaira.tech_sample,
        ]:
            if self.adata.uns[x] == UNS_STRING_META_IN_OBS:
                meta[x] = (np.sort(np.unique(self.adata.obs[x].values)),)
            else:
                meta[x] = self.adata.uns[x]
        # Add cell types into table if available:
        if self._adata_ids_sfaira.cell_ontology_class in self.adata.obs.keys():
            meta[self._adata_ids_sfaira.cell_ontology_class] = str((
                np.sort(np.unique(self.adata.obs[self._adata_ids_sfaira.cell_ontology_class].values)),
            ))
        else:
            meta[self._adata_ids_sfaira.cell_ontology_class] = " "
        meta.to_csv(fn_meta)

    def set_dataset_id(
            self,
            idx: int = 1
    ):
        def clean(s):
            if s is not None:
                s = s.replace(' ', '').replace('-', '').replace('_', '').lower()
            return s

        if self.sample_fn is not None:
            idx += self._sample_fns.index(self.sample_fn)
        idx = str(idx).zfill(3)

        if isinstance(self.author, List):
            author = self.author[0]
        else:
            author = self.author

        # Note: access private attributes here, e.g. _organism, to avoid loading of content via meta data, which would
        # invoke call to self.id before it is set.
        self.id = f"{clean(self._organism)}_" \
                  f"{clean(self._organ)}_" \
                  f"{self._year}_" \
                  f"{clean(self._assay)}_" \
                  f"{clean(author)}_" \
                  f"{idx}_" \
                  f"{self.doi}"

    # Properties:

    @property
    def age(self) -> Union[None, str]:
        if self._age is not None:
            return self._age
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.age in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.age]
            else:
                return None

    @age.setter
    def age(self, x: str):
        self.__erasing_protection(attr="age", val_old=self._age, val_new=x)
        self._value_protection(attr="age", allowed=self._ontology_container_sfaira.ontology_age, attempted=x)
        self._age = x

    @property
    def annotated(self) -> Union[bool, None]:
        if self.cellontology_id_obs_key is not None or self.cellontology_original_obs_key is not None:
            return True
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.annotated in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.annotated].values[0]
            elif self.loaded:
                # If data set was loaded and there is still no annotation indicated, it is declared unannotated.
                return False
            else:
                # If data set was not yet loaded, it is unclear if annotation would be loaded in ._load(),
                # if also no meta data is available, we do not know the status of the data set.
                return None

    @property
    def assay(self) -> Union[None, str]:
        if self._assay is not None:
            return self._assay
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.assay in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.assay]
            else:
                return None

    @assay.setter
    def assay(self, x: str):
        self.__erasing_protection(attr="protocol", val_old=self._assay, val_new=x)
        self._value_protection(attr="protocol", allowed=self._ontology_container_sfaira.ontology_protocol,
                               attempted=x)
        self._assay = x

    @property
    def author(self) -> str:
        if self._author is not None:
            return self._author
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is None or self._adata_ids_sfaira.author not in self.meta.columns:
                raise ValueError("author must be set but was neither set in constructor nor in meta data")
            return self.meta[self._adata_ids_sfaira.author]

    @author.setter
    def author(self, x: str):
        self.__erasing_protection(attr="author", val_old=self._author, val_new=x)
        self._author = x

    @property
    def bio_sample(self) -> Union[None, str]:
        if self._bio_sample is not None:
            return self._bio_sample
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.bio_sample in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.bio_sample]
            else:
                return None

    @bio_sample.setter
    def bio_sample(self, x: str):
        self.__erasing_protection(attr="bio_sample", val_old=self._bio_sample, val_new=x)
        self._bio_sample = x

    @property
    def data_dir(self):
        # Data is either directly in user supplied directory or in a sub directory if the overall directory is managed
        # by sfaira: In this case, the sub directory is named after the doi of the data set.
        sfaira_path = os.path.join(self.data_dir_base, self.directory_formatted_doi)
        if os.path.exists(sfaira_path):
            return sfaira_path
        else:
            return self.data_dir_base

    @property
    def development_stage(self) -> Union[None, str]:
        if self._development_stage is not None:
            return self._development_stage
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.development_stage in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.development_stage]
            else:
                return None

    @development_stage.setter
    def development_stage(self, x: str):
        self.__erasing_protection(attr="dev_stage", val_old=self._development_stage, val_new=x)
        self._value_protection(attr="dev_stage", allowed=self._ontology_container_sfaira.ontology_dev_stage,
                               attempted=x)
        self._development_stage = x

    @property
    def doi(self) -> str:
        if self._doi is not None:
            return self._doi
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is None or self._adata_ids_sfaira.doi not in self.meta.columns:
                raise ValueError("doi must be set but was neither set in constructor nor in meta data")
            return self.meta[self._adata_ids_sfaira.doi]

    @doi.setter
    def doi(self, x: str):
        self.__erasing_protection(attr="doi", val_old=self._doi, val_new=x)
        self._doi = x

    @property
    def directory_formatted_doi(self) -> str:
        return "d" + "_".join("_".join("_".join(self.doi.split("/")).split(".")).split("-"))

    @property
    def download_url_data(self) -> Union[Tuple[List[str]], Tuple[List[None]]]:
        """
        Data download website(s).

        Save as tuple with single element, which is a list of all download websites relevant to dataset.
        :return:
        """
        if self._download_url_data is not None:
            x = self._download_url_data
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            x = self.meta[self._adata_ids_sfaira.download_url_data]
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        return x

    @download_url_data.setter
    def download_url_data(self, x: Union[str, None, List[str], Tuple[str], List[None], Tuple[None]]):
        self.__erasing_protection(attr="download_url_data", val_old=self._download_url_data, val_new=x)
        # Formats to tuple with single element, which is a list of all download websites relevant to dataset,
        # which can be used as a single element column in a pandas data frame.
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        self._download_url_data = (x,)

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
        #    x = self.meta[self._ADATA_IDS_SFAIRA.download_url_meta]
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        return x

    @download_url_meta.setter
    def download_url_meta(self, x: Union[str, None, List[str], Tuple[str], List[None], Tuple[None]]):
        self.__erasing_protection(attr="download_url_meta", val_old=self._download_url_meta, val_new=x)
        # Formats to tuple with single element, which is a list of all download websites relevant to dataset,
        # which can be used as a single element column in a pandas data frame.
        if isinstance(x, str) or x is None:
            x = [x]
        if isinstance(x, list):
            x = (x,)
        self._download_url_meta = (x,)

    @property
    def ethnicity(self) -> Union[None, str]:
        if self._ethnicity is not None:
            return self._ethnicity
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.ethnicity in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.ethnicity]
            else:
                return None

    @ethnicity.setter
    def ethnicity(self, x: str):
        self.__erasing_protection(attr="ethnicity", val_old=self._ethnicity, val_new=x)
        self._value_protection(attr="ethnicity", allowed=self._adata_ids_sfaira.ontology_ethnicity, attempted=x)
        self._ethnicity = x

    @property
    def healthy(self) -> Union[None, bool]:
        if self._healthy is not None:
            return self._healthy
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.healthy in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.healthy]
            else:
                return None

    @healthy.setter
    def healthy(self, x: bool):
        self.__erasing_protection(attr="healthy", val_old=self._healthy, val_new=x)
        self._healthy = x

    @property
    def healthy_state_healthy(self) -> str:
        return self._healthy_state_healthy

    @healthy_state_healthy.setter
    def healthy_state_healthy(self, x: str):
        self.__erasing_protection(attr="healthy_state_healthy", val_old=self._healthy_state_healthy, val_new=x)
        self._healthy_state_healthy = x

    @property
    def id(self) -> str:
        if self._id is not None:
            return self._id
        else:
            raise AttributeError(f"Dataset ID was not set in dataloader in {self.doi}, please ensure the dataloader "
                                 f"constructor of this dataset contains a call to self.set_dataset_id()")

    @id.setter
    def id(self, x: str):
        self.__erasing_protection(attr="id", val_old=self._id, val_new=x)
        self._id = x

    @property
    def individual(self) -> Union[None, str]:
        if self._individual is not None:
            return self._individual
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.individual in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.individual]
            else:
                return None

    @individual.setter
    def individual(self, x: str):
        self.__erasing_protection(attr="bio_sample", val_old=self._individual, val_new=x)
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
            x = self.meta[self._adata_ids_sfaira.ncells]
        return int(x)

    @property
    def normalization(self) -> Union[None, str]:
        if self._normalization is not None:
            return self._normalization
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.normalization in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.normalization]
            else:
                return None

    @normalization.setter
    def normalization(self, x: str):
        self.__erasing_protection(attr="normalization", val_old=self._normalization, val_new=x)
        self._value_protection(attr="normalization", allowed=self._ontology_container_sfaira.ontology_normalization,
                               attempted=x)
        self._normalization = x

    @property
    def age_obs_key(self) -> str:
        return self._age_obs_key

    @age_obs_key.setter
    def age_obs_key(self, x: str):
        self.__erasing_protection(attr="age_obs_key", val_old=self._age_obs_key, val_new=x)
        self._age_obs_key = x

    @property
    def assay_obs_key(self) -> str:
        return self._assay_obs_key

    @assay_obs_key.setter
    def assay_obs_key(self, x: str):
        self.__erasing_protection(attr="assay_obs_key", val_old=self._assay_obs_key, val_new=x)
        self._assay_obs_key = x

    @property
    def bio_sample_obs_key(self) -> str:
        return self._bio_sample_obs_key

    @bio_sample_obs_key.setter
    def bio_sample_obs_key(self, x: str):
        self.__erasing_protection(attr="bio_sample_obs_key", val_old=self._bio_sample_obs_key, val_new=x)
        self._bio_sample_obs_key = x

    @property
    def cellontology_id_obs_key(self) -> str:
        return self._cellontology_id_obs_key

    @cellontology_id_obs_key.setter
    def cellontology_id_obs_key(self, x: str):
        self.__erasing_protection(attr="cellontology_id_obs_key", val_old=self._cellontology_id_obs_key, val_new=x)
        self._cellontology_id_obs_key = x

    @property
    def cellontology_original_obs_key(self) -> str:
        return self._cellontology_original_obs_key

    @cellontology_original_obs_key.setter
    def cellontology_original_obs_key(self, x: str):
        self.__erasing_protection(attr="cellontology_original_obs_key", val_old=self._cellontology_original_obs_key,
                                  val_new=x)
        self._cellontology_original_obs_key = x

    @property
    def development_stage_obs_key(self) -> str:
        return self._development_stage_obs_key

    @development_stage_obs_key.setter
    def development_stage_obs_key(self, x: str):
        self.__erasing_protection(attr="dev_stage_obs_key", val_old=self._development_stage_obs_key, val_new=x)
        self._development_stage_obs_key = x

    @property
    def ethnicity_obs_key(self) -> str:
        return self._ethnicity_obs_key

    @ethnicity_obs_key.setter
    def ethnicity_obs_key(self, x: str):
        self.__erasing_protection(attr="ethnicity_obs_key", val_old=self._ethnicity_obs_key, val_new=x)
        self._ethnicity_obs_key = x

    @property
    def healthy_obs_key(self) -> str:
        return self._healthy_obs_key

    @healthy_obs_key.setter
    def healthy_obs_key(self, x: str):
        self.__erasing_protection(attr="healthy_obs_key", val_old=self._healthy_obs_key, val_new=x)
        self._healthy_obs_key = x

    @property
    def individual_obs_key(self) -> str:
        return self._individual_obs_key

    @individual_obs_key.setter
    def individual_obs_key(self, x: str):
        self.__erasing_protection(attr="individual_obs_key", val_old=self._individual_obs_key, val_new=x)
        self._individual_obs_key = x

    @property
    def organ_obs_key(self) -> str:
        return self._organ_obs_key

    @organ_obs_key.setter
    def organ_obs_key(self, x: str):
        self.__erasing_protection(attr="organ_obs_key", val_old=self._organ_obs_key, val_new=x)
        self._organ_obs_key = x

    @property
    def organism_obs_key(self) -> str:
        return self._organism_obs_key

    @organism_obs_key.setter
    def organism_obs_key(self, x: str):
        self.__erasing_protection(attr="organism_obs_key", val_old=self._organism_obs_key, val_new=x)
        self._organism_obs_key = x

    @property
    def sex_obs_key(self) -> str:
        return self._sex_obs_key

    @sex_obs_key.setter
    def sex_obs_key(self, x: str):
        self.__erasing_protection(attr="sex_obs_key", val_old=self._sex_obs_key, val_new=x)
        self._sex_obs_key = x

    @property
    def state_exact_obs_key(self) -> str:
        return self._state_exact_obs_key

    @state_exact_obs_key.setter
    def state_exact_obs_key(self, x: str):
        self.__erasing_protection(attr="state_exact_obs_key", val_old=self._state_exact_obs_key, val_new=x)
        self._state_exact_obs_key = x

    @property
    def tech_sample_obs_key(self) -> str:
        return self._tech_sample_obs_key

    @tech_sample_obs_key.setter
    def tech_sample_obs_key(self, x: str):
        self.__erasing_protection(attr="tech_sample_obs_key", val_old=self._tech_sample_obs_key, val_new=x)
        self._tech_sample_obs_key = x

    @property
    def organ(self) -> Union[None, str]:
        if self._organ is not None:
            return self._organ
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.organ in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.organ]
            else:
                return None

    @organ.setter
    def organ(self, x: str):
        self.__erasing_protection(attr="organ", val_old=self._organ, val_new=x)
        self._value_protection(attr="organ", allowed=self._ontology_container_sfaira.ontology_organ, attempted=x)
        self._organ = x

    @property
    def organism(self) -> Union[None, str]:
        if self._organism is not None:
            return self._organism
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.organism in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.organism]
            else:
                return None

    @organism.setter
    def organism(self, x: str):
        self.__erasing_protection(attr="organism", val_old=self._organism, val_new=x)
        self._value_protection(attr="organism", allowed=self._ontology_container_sfaira.ontology_organism, attempted=x)
        self._organism = x

    @property
    def sex(self) -> Union[None, str]:
        if self._sex is not None:
            return self._sex
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.sex in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.sex]
            else:
                return None

    @sex.setter
    def sex(self, x: str):
        self.__erasing_protection(attr="sex", val_old=self._sex, val_new=x)
        self._value_protection(attr="sex", allowed=self._ontology_container_sfaira.ontology_sex, attempted=x)
        self._sex = x

    @property
    def source(self) -> str:
        return self._source

    @source.setter
    def source(self, x: Union[str, None]):
        self.__erasing_protection(attr="source", val_old=self._source, val_new=x)
        self._source = x

    @property
    def state_exact(self) -> Union[None, str]:
        if self._state_exact is not None:
            return self._state_exact
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.state_exact in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.state_exact]
            else:
                return None

    @state_exact.setter
    def state_exact(self, x: str):
        self.__erasing_protection(attr="state_exact", val_old=self._state_exact, val_new=x)
        self._state_exact = x

    @property
    def tech_sample(self) -> Union[None, str]:
        if self._tech_sample is not None:
            return self._tech_sample
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.tech_sample in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.tech_sample]
            else:
                return None

    @tech_sample.setter
    def tech_sample(self, x: str):
        self.__erasing_protection(attr="tech_sample", val_old=self._tech_sample, val_new=x)
        self._tech_sample = x

    @property
    def var_ensembl_col(self) -> str:
        return self._var_ensembl_col

    @var_ensembl_col.setter
    def var_ensembl_col(self, x: str):
        self.__erasing_protection(attr="var_ensembl_col", val_old=self._var_ensembl_col, val_new=x)
        self._var_ensembl_col = x

    @property
    def var_symbol_col(self) -> str:
        return self._var_symbol_col

    @var_symbol_col.setter
    def var_symbol_col(self, x: str):
        self.__erasing_protection(attr="var_symbol_col", val_old=self._var_symbol_col, val_new=x)
        self._var_symbol_col = x

    @property
    def year(self) -> Union[None, int]:
        if self._year is not None:
            return self._year
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.year in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.year]
            else:
                return None

    @year.setter
    def year(self, x: int):
        self.__erasing_protection(attr="year", val_old=self._year, val_new=x)
        self._value_protection(attr="year", allowed=self._ontology_container_sfaira.ontology_year, attempted=x)
        self._year = x

    @property
    def ontology_celltypes(self):
        return self._ontology_container_sfaira.ontology_cell_types

    @property
    def ontology_organ(self):
        return self._ontology_container_sfaira.ontology_organ

    @property
    def celltypes_universe(self):
        if self._celltype_universe:
            self._celltype_universe = CelltypeUniverse(
                cl=self.ontology_celltypes,
                uberon=self._ontology_container_sfaira.ontology_organ,
                organism=self.organism,
            )
        return self._celltype_universe

    @property
    def cell_ontology_map(self) -> dict:
        return self._ontology_class_map

    @cell_ontology_map.setter
    def cell_ontology_map(self, x: pd.DataFrame):
        self.__erasing_protection(attr="ontology_class_map", val_old=self._ontology_class_map, val_new=x)
        assert x.shape[1] in [2, 3], f"{x.shape} in {self.id}"
        assert x.columns[0] == self._adata_ids_sfaira.classmap_source_key
        assert x.columns[1] == self._adata_ids_sfaira.classmap_target_key
        # Transform data frame into a mapping dictionary:
        self._ontology_class_map = dict(list(zip(
            x[self._adata_ids_sfaira.classmap_source_key].values.tolist(),
            x[self._adata_ids_sfaira.classmap_target_key].values.tolist()
        )))

    # Private methods:

    def __erasing_protection(self, attr, val_old, val_new):
        """
        This is called when a erasing protected attribute is set to check whether it was set before.

        :param attr: Attribute to be set.
        :param val_old: Old value for attribute to be set.
        :param val_new: New value for attribute to be set.
        """
        if val_old is not None:
            raise ValueError(f"attempted to set erasing protected attribute {attr}: "
                             f"previously was {str(val_old)}, attempted to set {str(val_new)}")

    def _value_protection(
            self,
            attr: str,
            allowed: Union[Ontology, bool, int, float, str, List[bool], List[int], List[float], List[str]],
            attempted
    ):
        """
        Check whether value is from set of allowed values.

        Does not check if allowed is None.

        :param attr: Attribut to set.
        :param allowed: Constraint for values of `attr`.
            Either ontology instance used to constrain entries, or list of allowed values.
        :param attempted: Value(s) to attempt to set in `attr`.
        :return:
        """
        if allowed is not None:
            if not isinstance(attempted, list) and not isinstance(attempted, tuple):
                attempted = [attempted]
            if isinstance(allowed, Ontology):
                for x in attempted:
                    allowed.validate_node(x)
            else:
                for x in attempted:
                    if x not in allowed:
                        raise ValueError(f"{x} is not a valid entry for {attr}, choose from: {str(allowed)}")

    def subset_cells(self, key, values):
        """
        Subset list of adata objects based on cell-wise properties.

        These keys are properties that are not available in lazy model and require loading first because the
        subsetting works on the cell-level: .adata are maintained but reduced to matches.

        :param key: Property to subset by. Options:

            - "age" points to self.obs_key_age
            - "cell_ontology_class" points to self.obs_key_cellontology_original
            - "dev_stage" points to self.obs_key_dev_stage
            - "ethnicity" points to self.obs_key_ethnicity
            - "healthy" points to self.obs_key_healthy
            - "organ" points to self.obs_key_organ
            - "organism" points to self.obs_key_organism
            - "protocol" points to self.obs_key_protocol
            - "sex" points to self.obs_key_sex
            - "state_exact" points to self.obs_key_state_exact
        :param values: Classes to overlap to.
        :return:
        """
        if not isinstance(values, list):
            values = [values]

        def get_subset_idx(samplewise_key, cellwise_key):
            obs_key = getattr(self, cellwise_key)
            sample_attr = getattr(self, samplewise_key)
            if sample_attr is not None and obs_key is None:
                if not isinstance(sample_attr, list):
                    sample_attr = [sample_attr]
                if np.any([x in values for x in sample_attr]):
                    idx = np.arange(1, self.ncells)
                else:
                    idx = np.array([])
            elif sample_attr is None and obs_key is not None:
                assert self.adata is not None, "adata was not yet loaded"
                values_found = self.adata.obs[obs_key].values
                idx = np.where([x in values for x in values_found])
            elif sample_attr is not None and obs_key is not None:
                assert False, f"both cell-wise and sample-wise attribute {samplewise_key} given"
            else:
                assert False, "no subset chosen"
            return idx

        idx_keep = get_subset_idx(samplewise_key="obs_key_" + key, cellwise_key=key)
        self.adata = self.adata[idx_keep, :].copy()


class DatasetGroup:
    """
    Container class that co-manages multiple data sets, removing need to call Dataset() methods directly through
    wrapping them.

    Example:

    #query loaders lung
    #from sfaira.dev.data.loaders.lung import DatasetGroupLung as DatasetGroup
    #dsg_humanlung = DatasetGroupHuman(path='path/to/data')
    #dsg_humanlung.load_all(match_to_reference='Homo_sapiens_GRCh38_97')
    #dsg_humanlung[some_id]
    #dsg_humanlung.adata
    """
    datasets: Dict[str, DatasetBase]

    def __init__(self, datasets: dict):
        self._adata_ids_sfaira = AdataIdsSfaira()
        self.datasets = datasets

    @property
    def _unknown_celltype_identifiers(self):
        return np.unqiue(np.concatenate([v._unknown_celltype_identifiers for _, v in self.datasets.items()]))

    def load(
            self,
            annotated_only: bool = False,
            remove_gene_version: bool = True,
            match_to_reference: Union[str, bool, None] = None,
            load_raw: bool = False,
            allow_caching: bool = True,
            processes: int = 1,
            func=None,
            kwargs_func: Union[None, dict] = None,
    ):
        """
        Load all datasets in group (option for temporary loading).

        Note: This method automatically subsets to the group to the data sets for which input files were found.

        This method also allows temporarily loading data sets to execute function on loaded data sets (supply func).
        In this setting, datasets are removed from memory after the function has been executed.

        :param annotated_only:
        :param processes: Processes to parallelise loading over. Uses python multiprocessing if > 1, for loop otherwise.
        :param func: Function to run on loaded datasets. map_fun should only take one argument, which is a Dataset
            instance. The return can be empty:

                def func(dataset, **kwargs_func):
                    # code manipulating dataset and generating output x.
                    return x
        :param kwargs_func: Kwargs of func.
        """
        args = [
            remove_gene_version,
            match_to_reference,
            load_raw,
            allow_caching,
            func,
            kwargs_func
        ]

        if processes > 1 and len(self.datasets.items()) > 1:  # multiprocessing parallelisation
            print(f"using python multiprocessing (processes={processes}), "
                  f"for easier debugging revert to sequential execution (processes=1)")
            with multiprocessing.Pool(processes=processes) as pool:
                res = pool.starmap(map_fn, [
                    (tuple([v] + args),)
                    for k, v in self.datasets.items() if v.annotated or not annotated_only
                ])
            # Clear data sets that were not successfully loaded because of missing data:
            for x in res:
                if x is not None:
                    print(x[1])
                    del self.datasets[x[0]]
        else:  # for loop
            datasets_to_remove = []
            for k, v in self.datasets.items():
                print(f"loading {k}")
                x = map_fn(tuple([v] + args))
                # Clear data sets that were not successfully loaded because of missing data:
                if x is not None:
                    warnings.warn(f"data set {k} not loaded")
                    datasets_to_remove.append(k)
            for k in datasets_to_remove:
                del self.datasets[k]

    load.__doc__ += load_doc

    def streamline(self, format: str = "sfaira", clean: bool = False):
        """
        Streamline the adata instance in each data set to output format.

        Output format are saved in ADATA_FIELDS* classes.

        :param format: Export format.

            - "sfaira"
            - "cellxgene"
        :param clean: Whether to delete non-streamlined fields.
        :return:
        """
        for x in self.ids:
            self.datasets[x].streamline(format=format, clean=clean)

    def fragment(self) -> Dict[str, anndata.AnnData]:
        """
        Fragment data sets into largest consistent parititions based on meta data.

        ToDo return this as a DatasetGroup again.
          the streamlined Datasets are similar to anndata instances here, worth considering whether to use anndata
          instead because it can be indexed.

        :return:
        """
        # TODO: assert that data is streamlined.
        datasets_new = {}
        for k, v in self.datasets.items():
            # Define fragments and fragment names.
            # Because the data is streamlined, fragments are partitions of the .obs space, excluding the cell-wise
            # annotation columns:
            #       - cellontology_class
            #       - cellontology_id
            #       - cellontology_original
            cols_exclude = ["cellontology_class", "cellontology_id", "cellontology_original"]
            tab = v.adata.obs.loc[:, [x not in cols_exclude for x in tab.columns]]
            tab_unique = tab.drop_duplicates()
            idx_sets = [
                np.where([np.all(tab_unique.iloc[i, :] == tab.iloc[j, :])[0] for j in range(tab.shape[0])])
                for i in range(tab_unique.shape[0])
            ]
            for i, x in enumerate(idx_sets):
                datasets_new[k + "_fragment" + str(i)] = v.adata[x, :]
        return datasets_new

    def load_tobacked(
            self,
            adata_backed: anndata.AnnData,
            genome: str,
            idx: List[np.ndarray],
            annotated_only: bool = False,
            load_raw: bool = False,
            allow_caching: bool = True,
    ):
        """
        Loads data set group into slice of backed anndata object.

        Subsets self.datasets to the data sets that were found. Note that feature space is automatically formatted as
        this is necessary for concatenation.

        :param adata_backed: Anndata instance to load into.
        :param genome: Genome container target genomes loaded.
        :param idx: Indices in adata_backed to write observations to. This can be used to immediately create a
            shuffled object. This has to be a list of the length of self.data, one index array for each dataset.
        :param annotated_only:
        :param load_raw: See .load().
        :param allow_caching: See .load().
        :return: New row index for next element to be written into backed anndata.
        """
        i = 0
        for x in self.ids:
            # if this is for celltype prediction, only load the data with have celltype annotation
            try:
                if self.datasets[x].annotated or not annotated_only:
                    self.datasets[x].load_tobacked(
                        adata_backed=adata_backed,
                        genome=genome,
                        idx=idx[i],
                        load_raw=load_raw,
                        allow_caching=allow_caching
                    )
                    i += 1
            except FileNotFoundError:
                del self.datasets[x]

    def write_ontology_class_map(
            self,
            fn,
            protected_writing: bool = True,
            **kwargs
    ):
        """
        Write cell type maps of free text cell types to ontology classes.

        :param fn: File name of csv to load class maps from.
        :param protected_writing: Only write if file was not already found.
        """
        tab = []
        for k, v in self.datasets.items():
            if v.annotated:
                labels_original = np.sort(np.unique(np.concatenate([
                    v.adata.obs[self._adata_ids_sfaira.cell_types_original].values
                ])))
                tab.append(v.celltypes_universe.prepare_celltype_map_tab(
                    source=labels_original,
                    match_only=False,
                    anatomical_constraint=v.organ,
                    include_synonyms=True,
                    omit_list=v._unknown_celltype_identifiers,
                    **kwargs
                ))
        if len(tab) == 0:
            warnings.warn("attempted to write ontology classmaps for group without annotated data sets")
        else:
            tab = pandas.concat(tab, axis=0)
            # Take out columns with the same source:
            tab = tab.loc[[x not in tab.iloc[:i, 0].values for i, x in enumerate(tab.iloc[:, 0].values)], :].copy()
            tab = tab.sort_values(self._adata_ids_sfaira.classmap_source_key)
            if not os.path.exists(fn) or not protected_writing:
                tab.to_csv(fn, index=False, sep="\t")

    def download(self, **kwargs):
        for _, v in self.datasets.items():
            v.download(**kwargs)

    @property
    def ids(self):
        return list(self.datasets.keys())

    @property
    def adata_ls(self):
        adata_ls = []
        for i in self.ids:
            if self.datasets[i] is not None:
                if self.datasets[i].adata is not None:
                    adata_ls.append(self.datasets[i].adata)
        return adata_ls

    @property
    def adata(self):
        if not self.adata_ls:
            return None
        self.streamline(format="sfaira", clean=True)
        adata_ls = self.adata_ls

        # .var entries are renamed and copied upon concatenation.
        # To preserve gene names in .var, the target gene names are copied into var_names and are then copied
        # back into .var.
        for adata in adata_ls:
            adata.var.index = adata.var[self._adata_ids_sfaira.gene_id_ensembl].tolist()
        if len(adata_ls) > 1:
            # TODO: need to keep this? -> yes, still catching errors here (March 2020)
            # Fix for loading bug: sometime concatenating sparse matrices fails the first time but works on second try.
            try:
                adata_concat = adata_ls[0].concatenate(
                    *adata_ls[1:],
                    join="outer",
                    batch_key=self._adata_ids_sfaira.dataset,
                    batch_categories=[i for i in self.ids if self.datasets[i].adata is not None]
                )
            except ValueError:
                adata_concat = adata_ls[0].concatenate(
                    *adata_ls[1:],
                    join="outer",
                    batch_key=self._adata_ids_sfaira.dataset,
                    batch_categories=[i for i in self.ids if self.datasets[i].adata is not None]
                )

            adata_concat.var[self._adata_ids_sfaira.gene_id_ensembl] = adata_concat.var.index

            if len(set([a.uns[self._adata_ids_sfaira.mapped_features] for a in adata_ls])) == 1:
                adata_concat.uns[self._adata_ids_sfaira.mapped_features] = \
                    adata_ls[0].uns[self._adata_ids_sfaira.mapped_features]
            else:
                adata_concat.uns[self._adata_ids_sfaira.mapped_features] = False
        else:
            adata_concat = adata_ls[0]
            adata_concat.obs[self._adata_ids_sfaira.dataset] = self.ids[0]

        adata_concat.var_names_make_unique()
        return adata_concat

    def obs_concat(self, keys: Union[list, None] = None):
        """
        Returns concatenation of all .obs.

        Uses union of all keys if keys is not provided.

        :param keys:
        :return:
        """
        if keys is None:
            keys = np.unique(np.concatenate([list(x.obs.columns) for x in self.adata_ls]))
        obs_concat = pandas.concat([pandas.DataFrame(dict(
            [
                (k, self.datasets[x].adata.obs[k]) if k in self.datasets[x].adata.obs.columns
                else (k, ["nan" for _ in range(self.datasets[x].adata.obs.shape[0])])
                for k in keys
            ] + [(self._adata_ids_sfaira.dataset, [x for _ in range(self.datasets[x].adata.obs.shape[0])])]
        )) for x in self.ids if self.datasets[x].adata is not None])
        return obs_concat

    def ncells_bydataset(self, annotated_only: bool = False) -> np.ndarray:
        cells = []
        for x in self.ids:
            # if this is for celltype prediction, only load the data with have celltype annotation
            try:
                if self.datasets[x].annotated or not annotated_only:
                    cells.append(self.datasets[x].ncells)
            except FileNotFoundError:
                del self.datasets[x]
        return np.asarray(cells)

    def ncells(self, annotated_only: bool = False):
        cells = self.ncells_bydataset(annotated_only=annotated_only)
        return np.sum(cells)

    @property
    def ontology_celltypes(self):
        organism = np.unique([v.organism for _, v in self.datasets.items()])
        if len(organism) > 1:
            # ToDo: think about whether this should be handled differently.
            warnings.warn("found more than one organism in group, this could cause problems with using a joined cell "
                          "type ontology. Using only the ontology of the first data set in the group.")
        return self.datasets[self.ids[0]].ontology_celltypes

    def project_celltypes_to_ontology(self):
        """
        Project free text cell type names to ontology based on mapping table.
        :return:
        """
        for _, v in self.datasets.items():
            v.project_celltypes_to_ontology()

    def subset(self, key, values):
        """
        Subset list of adata objects based on sample-wise properties.

        These keys are properties that are available in lazy model.
        Subsetting happens on .datasets.

        :param key: Property to subset by.
        :param values: Classes to overlap to.
        :return:
        """
        ids_del = []
        if not isinstance(values, list):
            values = [values]
        for x in self.ids:
            try:
                values_found = getattr(self.datasets[x], key)
                if not isinstance(values_found, list):
                    values_found = [values_found]
                if not np.any([xx in values for xx in values_found]):
                    ids_del.append(x)
            except AttributeError:
                raise ValueError(f"{key} not a valid property of data set object")
        for x in ids_del:
            del self.datasets[x]

    def subset_cells(self, key, values: Union[str, List[str]]):
        """
        Subset list of adata objects based on cell-wise properties.

        These keys are properties that are not available in lazy model and require loading first because the
        subsetting works on the cell-level: .adata are maintained but reduced to matches.

        :param key: Property to subset by. Options:

            - "age" points to self.obs_key_age
            - "cell_ontology_class" points to self.obs_key_cellontology_original
            - "dev_stage" points to self.obs_key_dev_stage
            - "ethnicity" points to self.obs_key_ethnicity
            - "healthy" points to self.obs_key_healthy
            - "organ" points to self.obs_key_organ
            - "organism" points to self.obs_key_organism
            - "protocol" points to self.obs_key_protocol
            - "sex" points to self.obs_key_sex
            - "state_exact" points to self.obs_key_state_exact
        :param values: Classes to overlap to.
        :return:
        """
        for x in self.ids:
            self.datasets[x].subset_cells(key=key, values=values)
            if self.datasets[x].ncells == 0:  # none left
                del self.datasets[x]


class DatasetGroupDirectoryOriented(DatasetGroup):

    _cwd: os.PathLike

    def __init__(
            self,
            file_base: str,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
    ):
        """
        Automatically collects Datasets from python files in directory.

        Uses a pre-built DatasetGroup if this is defined in a group.py file, otherwise, the DatasetGroup is initialised
        here.

        :param file_base:
        :param data_path:
        :param meta_path:
        :param cache_path:
        """
        # Collect all data loaders from files in directory:
        datasets = []
        self._cwd = os.path.dirname(file_base)
        dataset_module = str(self._cwd.split("/")[-1])
        loader_pydoc_path = "sfaira.data.dataloaders.loaders." if str(self._cwd.split("/")[-5]) == "sfaira" else \
            "sfaira_extension.data.dataloaders.loaders."
        if "group.py" in os.listdir(self._cwd):
            DatasetGroupFound = pydoc.locate(loader_pydoc_path + dataset_module + ".group.DatasetGroup")
            dsg = DatasetGroupFound(data_path=data_path, meta_path=meta_path, cache_path=cache_path)
            datasets.extend(list(dsg.datasets.values))
        else:
            for f in os.listdir(self._cwd):
                if os.path.isfile(os.path.join(self._cwd, f)):  # only files
                    # Narrow down to data set files:
                    if f.split(".")[-1] == "py" and f.split(".")[0] not in ["__init__", "base", "group"]:
                        datasets_f = []
                        file_module = ".".join(f.split(".")[:-1])
                        DatasetFound = pydoc.locate(loader_pydoc_path + dataset_module + "." + file_module + ".Dataset")
                        # Load objects from name space:
                        # - load(): Loading function that return anndata instance.
                        # - SAMPLE_FNS: File name list for DatasetBaseGroupLoadingManyFiles
                        load_func = pydoc.locate(loader_pydoc_path + dataset_module + "." + file_module + ".load")
                        sample_fns = pydoc.locate(loader_pydoc_path + dataset_module + "." + file_module +
                                                  ".SAMPLE_FNS")
                        fn_yaml = os.path.join(self._cwd, file_module + ".yaml")
                        fn_yaml = fn_yaml if os.path.exists(fn_yaml) else None
                        # Check for sample_fns and sample_ids in yaml:
                        if fn_yaml is not None:
                            assert os.path.exists(fn_yaml), f"did not find yaml {fn_yaml}"
                            yaml_vals = read_yaml(fn=fn_yaml)
                            if sample_fns is None and yaml_vals["meta"]["sample_fns"] is not None:
                                sample_fns = yaml_vals["meta"]["sample_fns"]
                        if sample_fns is None:
                            sample_fns = [None]
                        # Here we distinguish between class that are already defined and those that are not.
                        # The latter case arises if meta data are defined in YAMLs and _load is given as a function.
                        if DatasetFound is None:
                            for x in sample_fns:
                                datasets_f.append(
                                    DatasetBase(
                                        data_path=data_path,
                                        meta_path=meta_path,
                                        cache_path=cache_path,
                                        load_func=load_func,
                                        sample_fn=x,
                                        sample_fns=sample_fns if sample_fns != [None] else None,
                                        yaml_path=fn_yaml,
                                    )
                                )
                        else:
                            for x in sample_fns:
                                datasets_f.append(
                                    DatasetFound(
                                        data_path=data_path,
                                        meta_path=meta_path,
                                        cache_path=cache_path,
                                        load_func=load_func,
                                        sample_fn=x,
                                        sample_fns=sample_fns if sample_fns != [None] else None,
                                        yaml_path=fn_yaml,
                                    )
                                )
                        # Load cell type maps:
                        for x in datasets_f:
                            x.load_ontology_class_map(fn=os.path.join(self._cwd, file_module + ".tsv"))
                        datasets.extend(datasets_f)

        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)))

    def clean_ontology_class_map(self):
        """
        Finalises processed class maps of free text cell types to ontology classes.

        Checks that the assigned ontology class names appear in the ontology.
        Adds a third column with the corresponding ontology IDs into the file.

        :return:
        """
        for f in os.listdir(self._cwd):
            if os.path.isfile(os.path.join(self._cwd, f)):  # only files
                # Narrow down to data set files:
                if f.split(".")[-1] == "py" and f.split(".")[0] not in ["__init__", "base", "group"]:
                    file_module = ".".join(f.split(".")[:-1])
                    fn_map = os.path.join(self._cwd, file_module + ".tsv")
                    if os.path.exists(fn_map):
                        # Access reading and value protection mechanisms from first data set loaded in group.
                        tab = list(self.datasets.values())[0]._read_class_map(fn=fn_map)
                        # Checks that the assigned ontology class names appear in the ontology.
                        list(self.datasets.values())[0]._value_protection(
                            attr="celltypes",
                            allowed=self.ontology_celltypes,
                            attempted=np.unique(tab[self._adata_ids_sfaira.classmap_target_key].values).tolist()
                        )
                        # Adds a third column with the corresponding ontology IDs into the file.
                        tab[self._adata_ids_sfaira.classmap_target_id_key] = [
                            self.ontology_celltypes.id_from_name(x) if x != self._adata_ids_sfaira.unknown_celltype_name
                            else self._adata_ids_sfaira.unknown_celltype_name
                            for x in tab[self._adata_ids_sfaira.classmap_target_key].values
                        ]
                        list(self.datasets.values())[0]._write_class_map(fn=fn_map, tab=tab)


class DatasetSuperGroup:
    """
    Container for multiple DatasetGroup instances.

    Used to manipulate structured dataset collections. Primarly designed for this manipulation, convert to DatasetGroup
    via flatten() for more functionalities.
    """
    adata: Union[None, anndata.AnnData]
    fn_backed: Union[None, PathLike]
    dataset_groups: Union[list, List[DatasetGroup], List[DatasetSuperGroup]]

    def __init__(self, dataset_groups: Union[None, List[DatasetGroup], List[DatasetSuperGroup]]):
        self.adata = None
        self.fn_backed = None
        self.set_dataset_groups(dataset_groups=dataset_groups)

        self._adata_ids_sfaira = AdataIdsSfaira()

    def set_dataset_groups(self, dataset_groups: Union[DatasetGroup, DatasetSuperGroup, List[DatasetGroup],
                                                       List[DatasetSuperGroup]]):
        if isinstance(dataset_groups, DatasetGroup) or isinstance(dataset_groups, DatasetSuperGroup):
            dataset_groups = [dataset_groups]
        if len(dataset_groups) > 0:
            if isinstance(dataset_groups[0], DatasetGroup):
                self.dataset_groups = dataset_groups
            elif isinstance(dataset_groups[0], DatasetSuperGroup):
                # Decompose super groups first
                dataset_groups_proc = []
                for x in dataset_groups:
                    dataset_groups_proc.extend(x.dataset_groups)
                self.dataset_groups = dataset_groups_proc
            else:
                assert False
        else:
            self.dataset_groups = []

    def extend_dataset_groups(self, dataset_groups: Union[List[DatasetGroup], List[DatasetSuperGroup]]):
        if isinstance(dataset_groups[0], DatasetGroup):
            self.dataset_groups.extend(dataset_groups)
        elif isinstance(dataset_groups[0], DatasetSuperGroup):
            # Decompose super groups first
            dataset_groups_proc = []
            for x in dataset_groups:
                dataset_groups_proc.extend(x.datasets)
            self.dataset_groups.extend(dataset_groups_proc)
        else:
            assert False

    def get_gc(
            self,
            genome: str = None
    ):
        if genome.lower().startswith("homo_sapiens"):
            g = SuperGenomeContainer(
                organism="human",
                genome=genome
            )
        elif genome.lower().startswith("mus_musculus"):
            g = SuperGenomeContainer(
                organism="mouse",
                genome=genome
            )
        else:
            raise ValueError(f"Genome {genome} not recognised. Needs to start with 'Mus_Musculus' or 'Homo_Sapiens'.")
        return g

    def ncells_bydataset(self, annotated_only: bool = False):
        """
        List of list of length of all data sets by data set group.
        :return:
        """
        return [x.ncells_bydataset(annotated_only=annotated_only) for x in self.dataset_groups]

    def ncells_bydataset_flat(self, annotated_only: bool = False):
        """
        Flattened list of length of all data sets.
        :return:
        """
        return [xx for x in self.ncells_bydataset(annotated_only=annotated_only) for xx in x]

    def ncells(self, annotated_only: bool = False):
        return np.sum(self.ncells_bydataset_flat(annotated_only=annotated_only))

    def flatten(self) -> DatasetGroup:
        """
        Returns DatasetGroup (rather than self = DatasetSuperGroup) containing all listed data sets.

        :return:
        """
        ds = {}
        for x in self.dataset_groups:
            for k, v in x.datasets.items():
                assert k not in ds.keys(), f"{k} was duplicated in super group, purge duplicates before flattening"
                ds[k] = v
        return DatasetGroup(datasets=ds)

    def download(self, **kwargs):
        for x in self.dataset_groups:
            x.download(**kwargs)

    def load_all(
            self,
            annotated_only: bool = False,
            match_to_reference: Union[str, bool, None] = None,
            remove_gene_version: bool = True,
            load_raw: bool = False,
            allow_caching: bool = True,
            processes: int = 1,
    ):
        """
        Loads data set human into anndata object.

        :param annotated_only:
        :param match_to_reference: See .load().
        :param remove_gene_version: See .load().
        :param load_raw: See .load().
        :param allow_caching: See .load().
        :param processes: Processes to parallelise loading over. Uses python multiprocessing if > 1, for loop otherwise.
            Note: parallelises loading of each dataset group, but not across groups.
        :return:
        """
        for x in self.dataset_groups:
            x.load(
                annotated_only=annotated_only,
                remove_gene_version=remove_gene_version,
                match_to_reference=match_to_reference,
                load_raw=load_raw,
                allow_caching=allow_caching,
                processes=processes,
            )
        # Make sure that concatenate is not used on a None adata object:
        adatas = [x.adata for x in self.dataset_groups if x.adata is not None]
        if len(adatas) > 1:
            self.adata = adatas[0].adata.concatenate(
                *adatas[1:],
                join="outer",
                batch_key=self._adata_ids_sfaira.dataset_group
            )
        elif len(adatas) == 1:
            self.adata = adatas[0]
        else:
            warnings.warn("no anndata instances to concatenate")

    def load_all_tobacked(
            self,
            fn_backed: PathLike,
            genome: str,
            shuffled: bool = False,
            as_dense: bool = False,
            annotated_only: bool = False,
            load_raw: bool = False,
            allow_caching: bool = True,
    ):
        """
        Loads data set human into backed anndata object.

        Example usage:

            ds = DatasetSuperGroup([...])
            ds.load_all_tobacked(
                fn_backed="...",
                target_genome="...",
                annotated_only=False
            )
            adata_backed = anndata.read(ds.fn_backed, backed='r')
            adata_slice = ad_full[idx]

        :param fn_backed: File name to save backed anndata to temporarily.
        :param genome: ID of target genomes.
        :param shuffled: Whether to shuffle data when writing to backed.
        :param as_dense: Whether to load into dense count matrix.
        :param annotated_only:
        :param load_raw: See .load().
        :param allow_caching: See .load().
        """
        if shuffled and not as_dense:
            raise ValueError("cannot write backed shuffled and sparse")
        scatter_update = shuffled or as_dense
        self.fn_backed = fn_backed
        n_cells = self.ncells(annotated_only=annotated_only)
        gc = self.get_gc(genome=genome)
        n_genes = gc.ngenes
        if scatter_update:
            self.adata = anndata.AnnData(
                scipy.sparse.csr_matrix((n_cells, n_genes), dtype=np.float32)
            )  # creates an empty anndata object with correct dimensions that can be filled with cells from data sets
        else:
            self.adata = anndata.AnnData(
                scipy.sparse.csr_matrix((0, n_genes), dtype=np.float32)
            )
        self.adata.filename = fn_backed  # setting this attribute switches this anndata to a backed object
        # Note that setting .filename automatically redefines .X as dense, so we have to redefine it as sparse:
        if not as_dense:
            X = scipy.sparse.csr_matrix(self.adata.X)  # redefines this backed anndata as sparse
            X.indices = X.indices.astype(np.int64)
            X.indptr = X.indptr.astype(np.int64)
            self.adata.X = X
        keys = [
            self._adata_ids_sfaira.annotated,
            self._adata_ids_sfaira.assay,
            self._adata_ids_sfaira.author,
            self._adata_ids_sfaira.dataset,
            self._adata_ids_sfaira.cell_ontology_class,
            self._adata_ids_sfaira.development_stage,
            self._adata_ids_sfaira.normalization,
            self._adata_ids_sfaira.organ,
            self._adata_ids_sfaira.state_exact,
            self._adata_ids_sfaira.year,
        ]
        if scatter_update:
            self.adata.obs = pandas.DataFrame({
                k: ["nan" for _ in range(n_cells)] for k in keys
            })
        else:
            for k in keys:
                self.adata.obs[k] = []
        # Define index vectors to write to:
        idx_vector = np.arange(0, n_cells)
        if shuffled:
            np.random.shuffle(idx_vector)
        idx_ls = []
        row = 0
        ncells = self.ncells_bydataset(annotated_only=annotated_only)
        if np.all([len(x) == 0 for x in ncells]):
            raise ValueError("no datasets found")
        for x in ncells:
            temp_ls = []
            for y in x:
                temp_ls.append(idx_vector[row:(row + y)])
                row += y
            idx_ls.append(temp_ls)
        print("checking expected and received data set sizes, rerun meta data generation if mismatch is found:")
        print(self.ncells_bydataset(annotated_only=annotated_only))
        print([[len(x) for x in xx] for xx in idx_ls])
        for i, x in enumerate(self.dataset_groups):
            x.load_tobacked(
                adata_backed=self.adata,
                genome=genome,
                idx=idx_ls[i],
                annotated_only=annotated_only,
                load_raw=load_raw,
                allow_caching=allow_caching,
            )
        # If the sparse non-shuffled approach is used, make sure that self.adata.obs.index is unique() before saving
        if not scatter_update:
            self.adata.obs.index = pd.RangeIndex(0, len(self.adata.obs.index))
        # Explicitly write backed file to disk again to make sure that obs are included and that n_obs is set correctly
        self.adata.write()
        # Saving obs separately below is therefore no longer required (hence commented out)
        # fn_backed_obs = ".".join(self.fn_backed.split(".")[:-1]) + "_obs.csv"
        # self.adata.obs.to_csv(fn_backed_obs)

    def delete_backed(self):
        del self.adata
        self.adata = None
        os.remove(str(self.fn_backed))

    def load_cached_backed(self, fn: PathLike):
        self.adata = anndata.read(fn, backed='r')

    def streamline(self, format: str = "sfaira", clean: bool = False):
        """
        Streamline the adata instance in each group and each data set to output format.

        Output format are saved in ADATA_FIELDS* classes.

        :param format: Export format.

            - "sfaira"
            - "cellxgene"
        :param clean: Whether to delete non-streamlined fields.
        :return:
        """
        for x in self.dataset_groups:
            for xx in x.ids:
                x.datasets[xx].streamline(format=format, clean=clean)

    def subset(self, key, values):
        """
        Subset list of adata objects based on match to values in key property.

        These keys are properties that are available in lazy model.
        Subsetting happens on .datasets.

        :param key: Property to subset by.
        :param values: Classes to overlap to.
        :return:
        """
        for x in self.dataset_groups:
            x.subset(key=key, values=values)
        self.dataset_groups = [x for x in self.dataset_groups if x.datasets]  # Delete empty DatasetGroups

    def subset_cells(self, key, values: Union[str, List[str]]):
        """
        Subset list of adata objects based on cell-wise properties.

        These keys are properties that are not available in lazy model and require loading first because the
        subsetting works on the cell-level: .adata are maintained but reduced to matches.

        :param key: Property to subset by. Options:

            - "age" points to self.obs_key_age
            - "cell_ontology_class" points to self.obs_key_cellontology_original
            - "dev_stage" points to self.obs_key_dev_stage
            - "ethnicity" points to self.obs_key_ethnicity
            - "healthy" points to self.obs_key_healthy
            - "organ" points to self.obs_key_organ
            - "organism" points to self.obs_key_organism
            - "protocol" points to self.obs_key_protocol
            - "sex" points to self.obs_key_sex
            - "state_exact" points to self.obs_key_state_exact
        :param values: Classes to overlap to.
        :return:
        """
        for i in range(len(self.dataset_groups)):
            self.dataset_groups[i].subset_cells(key=key, values=values)

    def project_celltypes_to_ontology(self):
        """
        Project free text cell type names to ontology based on mapping table.
        :return:
        """
        for _, v in self.dataset_groups:
            v.project_celltypes_to_ontology()
