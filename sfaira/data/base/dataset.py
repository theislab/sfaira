from __future__ import annotations

import abc
import anndata
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

from sfaira.versions.genome_versions import SuperGenomeContainer
from sfaira.versions.metadata import Ontology, CelltypeUniverse
from sfaira.consts import AdataIds, AdataIdsSfaira, META_DATA_FIELDS, OCS
from sfaira.data.utils import collapse_matrix, read_yaml

UNS_STRING_META_IN_OBS = "__obs__"


load_doc = \
    """
    :param remove_gene_version: Remove gene version string from ENSEMBL ID so that different versions in different data sets are superimposed.
    :param match_to_reference: Reference genomes name or False to keep original feature space.
    :param load_raw: Loads unprocessed version of data if available in data loader.
    :param allow_caching: Whether to allow method to cache adata object for faster re-loading.
    """


def is_child(
        query,
        ontology: Union[Ontology, bool, int, float, str, List[bool], List[int], List[float], List[str]],
        ontology_parent=None,
) -> True:
    """
    Check whether value is from set of allowed values using ontology.

    :param query: Value to attempt to set, only yield a single value and not a list.
    :param ontology: Constraint for values.
        Either ontology instance used to constrain entries, or list of allowed values.
    :param ontology_parent: If ontology is a DAG, not only check if node is a DAG node but also whether it is a child
        of this parent node.
    :return: Whether attempted term is sub-term of allowed term in ontology
    """
    if ontology_parent is None and ontology is None:
        return True
    else:
        if isinstance(ontology, Ontology):
            if ontology_parent is None:
                return ontology.is_node(query)
            else:
                return ontology.is_a(query=query, reference=ontology_parent)
        elif ontology is None:
            return query == ontology_parent
        else:
            raise ValueError(f"did not recognize ontology type {type(ontology)}")


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
    _assay_sc: Union[None, str]
    _assay_differentiation: Union[None, str]
    _assay_type_differentiation: Union[None, str]
    _author: Union[None, str]
    _bio_sample: Union[None, str]
    _cell_line: Union[None, str]
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
    _sample_source: Union[None, str]
    _state_exact: Union[None, str]
    _bio_sample: Union[None, str]
    _year: Union[None, int]

    _age_obs_key: Union[None, str]
    _assay_sc_obs_key: Union[None, str]
    _assay_differentiation_obs_key: Union[None, str]
    _assay_type_differentiation_obs_key: Union[None, str]
    _assay_cell_line_obs_key: Union[None, str]
    _cellontology_class_obs_key: Union[None, str]
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
    _sample_source_obs_key: Union[None, str]
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
        self.ontology_container_sfaira = OCS  # Using a pre-instantiated version of this yields drastic speed-ups.

        self.adata = None
        self.meta = None
        self.genome = None
        self.data_dir_base = data_path
        self.meta_path = meta_path
        self.cache_path = cache_path

        self._age = None
        self._author = None
        self._assay_sc = None
        self._assay_differentiation = None
        self._assay_type_differentiation = None
        self._bio_sample = None
        self._cell_line = None
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
        self._sample_source = None
        self._sex = None
        self._source = None
        self._state_exact = None
        self._tech_sample = None
        self._year = None

        self._age_obs_key = None
        self._assay_sc_obs_key = None
        self._assay_differentiation_obs_key = None
        self._assay_type_differentiation_obs_key = None
        self._bio_sample_obs_key = None
        self._cell_line_obs_key = None
        self._cellontology_class_obs_key = None
        self._cellontology_id_obs_key = None
        self._cellontology_original_obs_key = None
        self._development_stage_obs_key = None
        self._ethnicity_obs_key = None
        self._healthy_obs_key = None
        self._individual_obs_key = None
        self._organ_obs_key = None
        self._organism_obs_key = None
        self._sample_source_obs_key = None
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
                if v is not None and k not in ["sample_fns", "dataset_index"]:
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
                    self.adata = self.load_func(data_dir=self.data_dir, sample_fn=self.sample_fn)
            else:
                self.adata = self.load_func(data_dir=self.data_dir, sample_fn=self.sample_fn)

        def _cached_writing(filename):
            if filename is not None:
                dir_cache = os.path.dirname(filename)
                if not os.path.exists(dir_cache):
                    os.makedirs(dir_cache)
                self.adata.write_h5ad(filename)

        if load_raw and allow_caching:
            self.adata = self.load_func(data_dir=self.data_dir, sample_fn=self.sample_fn)
            _cached_writing(self.cache_fn)
        elif load_raw and not allow_caching:
            self.adata = self.load_func(data_dir=self.data_dir, sample_fn=self.sample_fn)
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
        self._collapse_genes(remove_gene_version=remove_gene_version)
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
        if not ensembl_col and not (isinstance(match_to_reference, bool) and not match_to_reference):
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

        if not symbol_col and not (isinstance(match_to_reference, bool) and not match_to_reference):
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

    def _collapse_genes(self, remove_gene_version):
        """
        Remove version tag on ensembl gene ID so that different versions are superimposed downstream.

        :param remove_gene_version:
        :return:
        """
        if remove_gene_version:
            self.adata.var_names = [
                x.split(".")[0] for x in self.adata.var[self._adata_ids_sfaira.gene_id_index].values
            ]
        # Collapse if necessary:
        self.adata = collapse_matrix(adata=self.adata)

        self.adata.var[self._adata_ids_sfaira.gene_id_index] = self.adata.var_names
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

    def _set_metadata_in_adata(self, adata_ids: AdataIds):
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
            [self.age, adata_ids.age, self.age_obs_key, self.ontology_container_sfaira.age],
            [self.assay_sc, adata_ids.assay_sc, self.assay_sc_obs_key, self.ontology_container_sfaira.assay_sc],
            [self.assay_differentiation, adata_ids.assay_differentiation, self.assay_differentiation_obs_key,
             self.ontology_container_sfaira.assay_differentiation],
            [self.assay_type_differentiation, adata_ids.assay_type_differentiation,
             self.assay_type_differentiation_obs_key, self.ontology_container_sfaira.assay_type_differentiation],
            [self.bio_sample, adata_ids.bio_sample, self.bio_sample_obs_key, None],
            [self.cell_line, adata_ids.cell_line, self.cell_line_obs_key, adata_ids.cell_line],
            [self.development_stage, adata_ids.development_stage, self.development_stage_obs_key,
             self.ontology_container_sfaira.developmental_stage],
            [self.ethnicity, adata_ids.ethnicity, self.ethnicity_obs_key,
             self.ontology_container_sfaira.ethnicity],
            [self.individual, adata_ids.individual, self.individual_obs_key, None],
            [self.organ, adata_ids.organ, self.organ_obs_key, self.ontology_container_sfaira.organ],
            [self.organism, adata_ids.organism, self.organism_obs_key,
             self.ontology_container_sfaira.organism],
            [self.sample_source, adata_ids.sample_source, self.sample_source_obs_key,
             self.ontology_container_sfaira.sample_source],
            [self.sex, adata_ids.sex, self.sex_obs_key, self.ontology_container_sfaira.sex],
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
                        attr=y, allowed=v, attempted=np.unique(self.adata.obs[z].values).tolist())
                    self.adata.obs[y] = self.adata.obs[z].values.tolist()
            else:
                assert False, "switch option should not occur"
        # Load boolean labels:
        for x, y, z, v, w in (
            [self.healthy, adata_ids.healthy, self.healthy_obs_key, self.ontology_container_sfaira.healthy,
             self.healthy_state_healthy],
        ):
            if x is None and z is None:
                self.adata.uns[y] = None
            elif x is not None and z is None:
                # Attribute supplied per data set: Write into .uns.
                if w is None:
                    self.adata.uns[y] = x
                else:
                    self.adata.uns[y] = x == w
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
                    label_y = self.adata.obs[z].values
                    # Use reference string to establish equality if available:
                    if w is not None:
                        label_y = label_y == w
                    self._value_protection(
                        attr=y, allowed=v, attempted=np.unique(label_y).tolist())
                    self.adata.obs[y] = label_y.tolist()
            else:
                assert False, "switch option should not occur"
        # Set cell-wise attributes (.obs):
        # None so far other than celltypes.
        # Set cell types:
        # Map cell type names from raw IDs to ontology maintained ones:
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
        elif format == "cellxgene":
            from sfaira.consts import AdataIdsCellxgene
            adata_fields = AdataIdsCellxgene()
        else:
            raise ValueError(f"did not recognize format {format}")
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
            self.adata.uns = dict([
                (k, v) for k, v in self.adata.uns.items() if k in [
                    adata_fields.annotated,
                    adata_fields.author,
                    adata_fields.doi,
                    adata_fields.download_url_data,
                    adata_fields.download_url_meta,
                    adata_fields.id,
                    adata_fields.normalization,
                    adata_fields.year,
                ]
            ])
            # Only retain target elements in adata.var:
            self.adata.var = self.adata.var.reindex(columns=[
                adata_fields.gene_id_names,
                adata_fields.gene_id_ensembl,
            ])
            # Only retain target columns in adata.obs:
            self.adata.obs = self.adata.obs.reindex(columns=[
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
            ])

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
        labels_original = self.adata.obs[self.cellontology_original_obs_key].values
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
            attempted=[
                x for x in np.unique(labels_mapped).tolist()
                if x != self._adata_ids_sfaira.unknown_celltype_identifier and
                x != self._adata_ids_sfaira.not_a_cell_celltype_identifier
            ]
        )
        self.adata.obs[self._adata_ids_sfaira.cell_ontology_class] = labels_mapped
        self.cellontology_class_obs_key = self._adata_ids_sfaira.cell_ontology_class
        self.adata.obs[self._adata_ids_sfaira.cell_types_original] = labels_original
        # Add cell type IDs into object:
        # The IDs are not read from a source file but inferred based on the class name.
        # TODO this could be changed in the future, this allows this function to be used both on cell type name mapping
        #  files with and without the ID in the third column.
        ids_mapped = [
            self.ontology_container_sfaira.cellontology_class.id_from_name(x)
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
            self._adata_ids_sfaira.assay_sc,
            self._adata_ids_sfaira.assay_differentiation,
            self._adata_ids_sfaira.assay_type_differentiation,
            self._adata_ids_sfaira.bio_sample,
            self._adata_ids_sfaira.cell_line,
            self._adata_ids_sfaira.development_stage,
            self._adata_ids_sfaira.ethnicity,
            self._adata_ids_sfaira.healthy,
            self._adata_ids_sfaira.individual,
            self._adata_ids_sfaira.organ,
            self._adata_ids_sfaira.organism,
            self._adata_ids_sfaira.sample_source,
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
                  f"{clean(self._assay_sc)}_" \
                  f"{clean(author)}_" \
                  f"{idx}_" \
                  f"{self.doi_main}"

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
        self._value_protection(attr="age", allowed=self.ontology_container_sfaira.age, attempted=x)
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
    def assay_sc(self) -> Union[None, str]:
        if self._assay_sc is not None:
            return self._assay_sc
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.assay_sc in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.assay_sc]
            else:
                return None

    @assay_sc.setter
    def assay_sc(self, x: str):
        self.__erasing_protection(attr="assay_sc", val_old=self._assay_sc, val_new=x)
        self._value_protection(attr="assay_sc", allowed=self.ontology_container_sfaira.assay_sc, attempted=x)
        self._assay_sc = x

    @property
    def assay_differentiation(self) -> Union[None, str]:
        if self._assay_differentiation is not None:
            return self._assay_differentiation
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.assay_differentiation in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.assay_differentiation]
            else:
                return None

    @assay_differentiation.setter
    def assay_differentiation(self, x: str):
        self.__erasing_protection(attr="assay_differentiation", val_old=self._assay_differentiation, val_new=x)
        self._value_protection(attr="assay_differentiation",
                               allowed=self.ontology_container_sfaira.assay_differentiation, attempted=x)
        self._assay_differentiation = x

    @property
    def assay_type_differentiation(self) -> Union[None, str]:
        if self._assay_type_differentiation is not None:
            return self._assay_type_differentiation
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.assay_type_differentiation in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.assay_type_differentiation]
            else:
                return None

    @assay_type_differentiation.setter
    def assay_type_differentiation(self, x: str):
        self.__erasing_protection(attr="assay_type_differentiation", val_old=self._assay_type_differentiation, val_new=x)
        self._value_protection(attr="assay_type_differentiation",
                               allowed=self.ontology_container_sfaira.assay_type_differentiation, attempted=x)
        self._assay_type_differentiation = x

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
    def cell_line(self) -> Union[None, str]:
        if self._cell_line is not None:
            return self._cell_line
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.cell_line in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.cell_line]
            else:
                return None

    @cell_line.setter
    def cell_line(self, x: str):
        self.__erasing_protection(attr="cell_line", val_old=self._cell_line, val_new=x)
        self._cell_line = x

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
        self._value_protection(attr="dev_stage", allowed=self.ontology_container_sfaira.developmental_stage,
                               attempted=x)
        self._development_stage = x

    @property
    def doi(self) -> Union[str, List[str]]:
        if self._doi is not None:
            return self._doi
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is None or self._adata_ids_sfaira.doi not in self.meta.columns:
                raise ValueError("doi must be set but was neither set in constructor nor in meta data")
            return self.meta[self._adata_ids_sfaira.doi]

    @doi.setter
    def doi(self, x: Union[str, List[str]]):
        self.__erasing_protection(attr="doi", val_old=self._doi, val_new=x)
        self._doi = x

    @property
    def doi_main(self) -> str:
        """
        Yields the main DOI associated with the study, defined as the DOI that comes first in alphabetical order.
        """
        return self.doi if isinstance(self.doi, str) else np.sort(self.doi)[0]

    @property
    def directory_formatted_doi(self) -> str:
        # Chose first doi in list.
        return "d" + "_".join("_".join("_".join(self.doi_main.split("/")).split(".")).split("-"))

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
            raise AttributeError(f"Dataset ID was not set in dataloader in {self.doi_main}, please ensure the "
                                 f"dataloader constructor of this dataset contains a call to self.set_dataset_id()")

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
        self._value_protection(attr="normalization", allowed=self.ontology_container_sfaira.normalization,
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
    def assay_sc_obs_key(self) -> str:
        return self._assay_sc_obs_key

    @assay_sc_obs_key.setter
    def assay_sc_obs_key(self, x: str):
        self.__erasing_protection(attr="assay_sc_obs_key", val_old=self._assay_sc_obs_key, val_new=x)
        self._assay_sc_obs_key = x

    @property
    def assay_differentiation_obs_key(self) -> str:
        return self._assay_differentiation_obs_key

    @assay_differentiation_obs_key.setter
    def assay_differentiation_obs_key(self, x: str):
        self.__erasing_protection(attr="assay_differentiation_obs_key", val_old=self._assay_differentiation_obs_key, val_new=x)
        self._assay_differentiation_obs_key = x

    @property
    def assay_type_differentiation_obs_key(self) -> str:
        return self._assay_type_differentiation_obs_key

    @assay_type_differentiation_obs_key.setter
    def assay_type_differentiation_obs_key(self, x: str):
        self.__erasing_protection(attr="assay_type_differentiation_otype_bs_key", val_old=self._assay_differentiation_obs_key, val_new=x)
        self._assay_type_differentiation_obs_key = x

    @property
    def bio_sample_obs_key(self) -> str:
        return self._bio_sample_obs_key

    @bio_sample_obs_key.setter
    def bio_sample_obs_key(self, x: str):
        self.__erasing_protection(attr="bio_sample_obs_key", val_old=self._bio_sample_obs_key, val_new=x)
        self._bio_sample_obs_key = x

    @property
    def cell_line_obs_key(self) -> str:
        return self._cell_line_obs_key

    @cell_line_obs_key.setter
    def cell_line_obs_key(self, x: str):
        self.__erasing_protection(attr="cell_line_obs_key", val_old=self._cell_line_obs_key, val_new=x)
        self._cell_line_obs_key = x

    @property
    def cellontology_class_obs_key(self) -> str:
        return self._cellontology_class_obs_key

    @cellontology_class_obs_key.setter
    def cellontology_class_obs_key(self, x: str):
        self._cellontology_class_obs_key = x

    @property
    def cellontology_id_obs_key(self) -> str:
        return self._cellontology_id_obs_key

    @cellontology_id_obs_key.setter
    def cellontology_id_obs_key(self, x: str):
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
    def sample_source_obs_key(self) -> str:
        return self._sample_source_obs_key

    @sample_source_obs_key.setter
    def sample_source_obs_key(self, x: str):
        self.__erasing_protection(attr="sample_source_obs_key", val_old=self._sample_source_obs_key, val_new=x)
        self._sample_source_obs_key = x

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
        self._value_protection(attr="organ", allowed=self.ontology_container_sfaira.organ, attempted=x)
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
        self._value_protection(attr="organism", allowed=self.ontology_container_sfaira.organism, attempted=x)
        self._organism = x

    @property
    def sample_source(self) -> Union[None, str]:
        if self._sample_source is not None:
            return self._sample_source
        else:
            if self.meta is None:
                self.load_meta(fn=None)
            if self.meta is not None and self._adata_ids_sfaira.sample_source in self.meta.columns:
                return self.meta[self._adata_ids_sfaira.sample_source]
            else:
                return None

    @sample_source.setter
    def sample_source(self, x: str):
        self.__erasing_protection(attr="sample_source", val_old=self._sample_source, val_new=x)
        self._value_protection(attr="sample_source", allowed=self.ontology_container_sfaira.sample_source, attempted=x)
        self._sample_source = x

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
        self._value_protection(attr="sex", allowed=self.ontology_container_sfaira.sex, attempted=x)
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
        self._value_protection(attr="year", allowed=self.ontology_container_sfaira.year, attempted=x)
        self._year = x

    @property
    def ontology_celltypes(self):
        return self.ontology_container_sfaira.cellontology_class

    @property
    def ontology_organ(self):
        return self.ontology_container_sfaira.organ

    @property
    def celltypes_universe(self):
        if self._celltype_universe:
            self._celltype_universe = CelltypeUniverse(
                cl=self.ontology_celltypes,
                uberon=self.ontology_container_sfaira.organ,
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

        :param attr: Attribute to set.
        :param allowed: Constraint for values of `attr`.
            Either ontology instance used to constrain entries, or list of allowed values.
        :param attempted: Value(s) to attempt to set in `attr`.
        :return:
        """
        if isinstance(attempted, np.ndarray):
            attempted = attempted.tolist()
        if isinstance(attempted, tuple):
            attempted = list(attempted)
        if not isinstance(attempted, list):
            attempted = [attempted]
        for x in attempted:
            if not is_child(query=x, ontology=allowed):
                if isinstance(allowed, Ontology):
                    # use node names instead of ontology object to produce a readable error message
                    allowed = allowed.node_names
                raise ValueError(f"{x} is not a valid entry for {attr}, choose from: {allowed}")

    def subset_cells(self, key, values):
        """
        Subset list of adata objects based on cell-wise properties.

        These keys are properties that are not available in lazy model and require loading first because the
        subsetting works on the cell-level: .adata are maintained but reduced to matches.

        :param key: Property to subset by. Options:

            - "age" points to self.age_obs_key
            - "assay_sc" points to self.assay_sc_obs_key
            - "assay_differentiation" points to self.assay_differentiation_obs_key
            - "assay_type_differentiation" points to self.assay_type_differentiation_obs_key
            - "cell_line" points to self.cell_line
            - "cellontology_class" points to self.cellontology_class_obs_key
            - "developmental_stage" points to self.developmental_stage_obs_key
            - "ethnicity" points to self.ethnicity_obs_key
            - "healthy" points to self.healthy_obs_key
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
            except AttributeError:
                sample_attr = None
            obs_key = getattr(self, cellwise_key)
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
            elif sample_attr is not None and obs_key is not None:
                assert False, f"both cell-wise and sample-wise attribute {samplewise_key} given"
            else:
                assert False, "no subset chosen"
            return idx

        idx_keep = get_subset_idx(samplewise_key=key, cellwise_key=key + "_obs_key")
        self.adata = self.adata[idx_keep, :].copy()  # if len(idx_keep) > 0 else None
