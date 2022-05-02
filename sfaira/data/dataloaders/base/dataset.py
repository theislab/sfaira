from __future__ import annotations

import anndata
import numpy as np
import pandas as pd
import os
from os import PathLike
import pandas
import scipy.sparse
from typing import Dict, List, Union
import warnings

from sfaira.versions.genomes import GenomeContainer
from sfaira.versions.metadata import CelltypeUniverse
from sfaira.consts import AdataIdsCellxgene_v2_0_0, AdataIdsSfaira, META_DATA_FIELDS, OCS
from sfaira.data.dataloaders.base.utils import is_child, get_directory_formatted_doi, identify_tsv
from sfaira.data.dataloaders.crossref import crossref_query
from sfaira.data.dataloaders.download_utils import download
from sfaira.data.dataloaders.obs_utils import streamline_obs_uns
from sfaira.data.store.io.io_dao import write_dao
from sfaira.data.dataloaders.base.annotation_container import AnnotationContainer
from sfaira.data.dataloaders.var_utils import streamline_var, collapse_x_var_by_feature
from sfaira.consts.utils import clean_doi, clean_id_str
from sfaira.versions.metadata.maps import prepare_ontology_map_tab


class DatasetBase(AnnotationContainer):
    adata: Union[None, anndata.AnnData] = None
    _meta: Union[None, pandas.DataFrame] = None
    data_dir_base: Union[None, str]
    meta_path: Union[None, str]
    cache_path: Union[None, str]
    genome_container: Union[None, GenomeContainer] = None
    supplier: str

    _celltype_universe: Union[None, CelltypeUniverse] = None
    _ontology_class_maps: Union[dict]

    load_raw: Union[None, bool] = None
    mapped_features: Union[None, str, bool] = None
    remove_gene_version: Union[None, bool] = None
    subset_gene_type: Union[None, str] = None
    streamlined_meta: bool = False

    sample_fn: Union[None, str]
    _sample_fns: Union[None, List[str]]

    _title: Union[None, str] = None

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

        self.data_dir_base = data_path
        self.meta_path = meta_path
        self.cache_path = cache_path

        self._ontology_class_maps = dict([(k, None) for k in self._adata_ids.ontology_constrained])

        self.sample_fn = sample_fn
        self._sample_fns = sample_fns

        # Check if YAML files exists, read meta data from there if available:
        if yaml_path is not None:
            self.read_from_yaml(yaml_path=yaml_path, sample_fn=self.sample_fn)
            self.set_dataset_id(idx=self.dataset_index)

        self.load_func = load_func
        self.dict_load_func_annotation = dict_load_func_annotation
        self._additional_annotation_key = additional_annotation_key

        self.supplier = "sfaira"

    def clear(self):
        """
        Remove loaded .adata to reduce memory footprint.

        :return:
        """
        import gc
        self.adata = None
        gc.collect()

    def download(self, **kwargs):
        assert self.download_url_data is not None, f"The `download_url_data` attribute of dataset {self.id_base} " \
                                                   f"is not set, cannot download dataset."
        assert self.data_dir_base is not None, "No path was provided when instantiating the dataset container, " \
                                               "cannot download datasets."

        if not os.path.exists(os.path.join(self.data_dir_base, self.directory_formatted_doi)):
            os.makedirs(os.path.join(self.data_dir_base, self.directory_formatted_doi))

        urls = self.download_url_data[0] + self.download_url_meta[0]
        download(urls=urls, data_dir=self.data_dir, directory_formatted_doi=self.directory_formatted_doi,
                 dataset_id=self.id_base, **kwargs)

    @property
    def cache_fn(self):
        if self.cache_path is None:
            cache = self.data_dir
        else:
            cache = os.path.join(self.cache_path, self.directory_formatted_doi)
        return os.path.join(cache, "cache", self.id + ".h5ad")

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
            raise ValueError(f"adata of {self.id_base} already loaded.")
        if self.data_dir is None:
            raise ValueError("No sfaira data repo path provided in constructor.")

        def _error_buffered_reading(**load_kwargs):
            self.adata = self.load_func(data_dir=self.data_dir, sample_fn=self.sample_fn, **load_kwargs)

        # Run data set-specific loading script:
        def _assembly_wrapper():
            if self.load_func is None:
                raise ValueError(f"Tried to access load_func for {self.id_base} but did not find any.")
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

    def collapse_counts(self):
        """
        Collapse count matrix along duplicated index.
        """
        if len(np.unique(self.adata.var.index)) < self.adata.var.shape[0]:
            obs_names = self.adata.obs_names
            x, var = collapse_x_var_by_feature(x=self.adata.X, var=self.adata.var, var_column="index")
            self.adata = anndata.AnnData(
                X=x,
                obs=self.adata.obs,
                obsm=self.adata.obsm,
                var=var,
                uns=self.adata.uns
            )
            self.adata.obs_names = obs_names

    def streamline_var(
            self,
            match_to_release: Union[str, Dict[str, str], None] = None,
            clean_var: bool = True,
            remove_gene_version: bool = True,
            schema: str = "sfaira",
            subset_genes_to_type: Union[None, str, List[str]] = None,
            verbose: int = 1,
            **kwargs
    ):
        """
        Subset and sort genes to genes defined in an assembly or genes of a particular type, such as protein coding.
        This also adds missing ensid or gene symbol columns if match_to_reference is not set to False and removes all
        adata.var columns that are not defined as gene_id_ensembl_var_key or gene_id_symbol_var_key in the dataloader.

        :param clean_var: Whether to delete non-streamlined fields in .var, .varm and .varp.
        :param match_to_release: Which ensembl genome annotation release to map the feature space to.
            Uses default set in schema if not set. Can be:
                - str: Provide the name of the release (eg. "104").
                - dict: Mapping of organism to name of the release (see str format). Chooses release for each
                    data set based on organism annotation.
        :param remove_gene_version: Whether to remove the version number after the colon sometimes found in ensembl
            gene ids.
            Uses default set in schema if not set.
        :param schema: Export format.
            - "sfaira"
            - "cellxgene"
        :param subset_genes_to_type: Type(s) to subset to. Can be a single type or a list of types or None.
            Uses default set in schema if not set.
            Types can be:
                - None: Keep the subset of input gene set that can be mapped to assembly.
                - "all": All genes in assembly.
                - "protein_coding": All protein coding genes in assembly.
        :param verbose: Report feature transformation statistics.
        """
        self.__assert_loaded()
        if schema.startswith("sfaira"):
            adata_target_ids = AdataIdsSfaira()
        elif schema.startswith("cellxgene"):
            adata_target_ids = AdataIdsCellxgene_v2_0_0()
        else:
            raise ValueError(f"did not recognize schema {schema}")
        if isinstance(match_to_release, dict):
            match_to_release = match_to_release[self.organism]

        streamline_output = streamline_var(
            adata=self.adata,
            adata_target_ids=adata_target_ids,
            clean_var=clean_var,
            dataset_id=self.id,
            organism=self.organism,
            feature_id_var_key=self.feature_id_var_key,
            feature_symbol_var_key=self.feature_symbol_var_key,
            layer_counts=self.layer_counts,
            layer_processed=self.layer_processed,
            match_to_release=match_to_release,
            remove_gene_version=remove_gene_version,
            subset_genes_to_type=subset_genes_to_type,
            verbose=verbose,
            **kwargs)
        self.adata = streamline_output["adata"]
        self.genome_container = streamline_output["genome_container"]
        self.layer_counts = streamline_output["layer_counts"]
        self.layer_processed = streamline_output["layer_processed"]
        self.mapped_features = self.genome_container.release
        self.remove_gene_version = remove_gene_version
        self.subset_gene_type = subset_genes_to_type
        self.feature_id_var_key = adata_target_ids.feature_id
        self.feature_symbol_var_key = adata_target_ids.feature_symbol

    def streamline_obs_uns(
            self,
            schema: str = "sfaira",
            clean_obs: Union[bool, list, np.ndarray, tuple] = True,
            clean_uns: bool = True,
            clean_obs_names: bool = True,
            keep_orginal_obs: Union[bool, list, np.ndarray, tuple] = False,
            keep_symbol_obs: bool = True,
            keep_id_obs: bool = True,
            **kwargs
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
             Alternatively, list of .obs fields to remove.
        :param clean_uns: Whether to delete non-streamlined fields in .uns.
        :param clean_obs_names: Whether to replace obs_names with a string comprised of dataset id and an increasing
            integer.
        :param keep_orginal_obs: For ontology-constrained .obs columns, whether to keep a column with original
            annotation. Alternatively, list of original fields to keep, others are removed.
        :param keep_symbol_obs: For ontology-constrained .obs columns, whether to keep a column with ontology symbol
            annotation.
        :param keep_id_obs: For ontology-constrained .obs columns, whether to keep a column with ontology ID annotation.
        :return:
        """
        self.__assert_loaded()
        if isinstance(keep_orginal_obs, tuple):
            keep_orginal_obs = list(keep_orginal_obs)
        elif isinstance(keep_orginal_obs, np.ndarray):
            keep_orginal_obs = keep_orginal_obs.tolist()

        # Set schema as provided by the user
        if schema.startswith("sfaira"):
            adata_target_ids = AdataIdsSfaira()
        elif schema.startswith("cellxgene"):
            adata_target_ids = AdataIdsCellxgene_v2_0_0()
        else:
            raise ValueError(f"did not recognize schema {schema}")

        self.adata = streamline_obs_uns(adata=self.adata,
                                        adata_input_ids=self._adata_ids,
                                        adata_target_ids=adata_target_ids,
                                        annotation_container=self,
                                        organism=self.organism,
                                        schema=schema,
                                        clean_obs=clean_obs,
                                        clean_uns=clean_uns,
                                        clean_obs_names=clean_obs_names,
                                        dataset_id=self.id,
                                        keep_orginal_obs=keep_orginal_obs,
                                        keep_symbol_obs=keep_symbol_obs,
                                        keep_id_obs=keep_id_obs,
                                        ontology_class_maps=self.ontology_class_maps,
                                        **kwargs)
        self._adata_ids = adata_target_ids
        self.streamlined_meta = True

    def write_distributed_store(
            self,
            dir_cache: Union[str, os.PathLike],
            store_format: str = "dao",
            dense: bool = False,
            compression_kwargs: Union[dict, None] = None,
            chunks: Union[int, None] = None,
            shuffle_data: bool = False
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
        :param shuffle_data: If True -> shuffle ordering of cells in datasets before writing store
        """
        if compression_kwargs is None:
            compression_kwargs = {}
        self.__assert_loaded()
        if store_format == "h5ad":
            if not isinstance(self.adata.X, scipy.sparse.csr_matrix):
                print(f"WARNING: high-perfomances caches based on .h5ad work better with .csr formatted expression "
                      f"data, found {type(self.adata.X)}")
            fn = os.path.join(dir_cache, self.id + ".h5ad")
            as_dense = ("X",) if dense else ()
            print(f"writing {self.adata.shape} into {fn}")
            self.adata.write_h5ad(filename=fn, as_dense=as_dense, **compression_kwargs)
        elif store_format == "dao":
            # Convert data object to sparse / dense as required:
            if not dense:
                raise ValueError("WARNING: sparse zarr array performance is not be optimal and not supported yet, "
                                 "consider writing as dense and consider that zarr arrays are compressed on disk!")
            fn = os.path.join(dir_cache, self.id)
            chunks = (chunks, self.adata.X.shape[1]) if chunks is not None else True
            write_dao(store=fn, adata=self.adata, chunks=chunks, compression_kwargs=compression_kwargs,
                      shuffle_data=shuffle_data)
        else:
            raise ValueError()

    def _write_ontology_class_map(self, fn, tab: pd.DataFrame):
        """
        Write class map to file.

        Helper to allow direct interaction with written table instead of using table from instance.

        :param fn: File name of csv to write class maps to.
        :param tab: Class map table.
        """
        tab.to_csv(fn, index=False, sep="\t")

    def write_ontology_class_maps(
            self,
            fn,
            attrs: List[str],
            protected_writing: bool = True,
            **kwargs
    ):
        """
        Load class maps of ontology-controlled field to ontology classes.

        TODO: deprecate and only keep DatasetGroup writing?

        :param fn: File name of tsv to write class maps to.
        :param attrs: Attributes to create a tsv for. Must correspond to *_obs_key in yaml.
        :param protected_writing: Only write if file was not already found.
        :return:
        """
        for x in attrs:
            k = getattr(self, x + "_obs_key")
            if k is None:
                warnings.warn(f"attempted to write ontology class maps for data set {self.id_base} without annotation "
                              f"of meta data {x}")
            elif k not in self.adata.obs.columns:
                warnings.warn(f"attempted to write ontology class maps for data set {self.id_base} but did not find "
                              f"column {k} in .obs which should correspond to {x}")
            else:
                fn_x = fn + "_" + x + ".tsv"
                labels_original = np.unique(self.adata.obs[k].values)
                tab = prepare_ontology_map_tab(
                    source=labels_original,
                    onto=x,
                    organism=self.organism,
                    match_only=False,
                    anatomical_constraint=self.organ,
                    include_synonyms=True,
                    omit_list=[self._adata_ids.unknown_metadata_identifier],
                    **kwargs
                )
                if not os.path.exists(fn_x) or not protected_writing:
                    self._write_ontology_class_map(fn=fn_x, tab=tab)

    def _read_ontology_class_map(self, fn) -> pd.DataFrame:
        """
        Read class map.

        Helper to allow direct interaction with resulting table instead of loading into instance.

        :param fn: File name of csv to load class maps from.
        :return:
        """
        try:
            # Need dtype="str" to force numeric cell type identifiers, e.g. cluster numbers to be in string format.
            tab = pd.read_csv(fn, header=0, index_col=None, sep="\t", dtype="str").astype(str)
        except pandas.errors.ParserError as e:
            print(f"{self.id_base}")
            raise pandas.errors.ParserError(e)
        return tab

    def read_ontology_class_maps(self, fns: List[str]):
        """
        Load class maps of free text class labels to ontology classes.

        :param fns: File names of tsv to load class maps from.
        :return:
        """
        assert isinstance(fns, list)
        ontology_class_maps = {}
        for fn in fns:
            k = identify_tsv(fn=fn, ontology_names=self._adata_ids.ontology_constrained)
            if os.path.exists(fn):
                ontology_class_maps[k] = self._read_ontology_class_map(fn=fn)
        self.ontology_class_maps = ontology_class_maps

    @property
    def citation(self):
        """
        Return all information necessary to cite data set.

        :return:
        """
        return [self.author, self.year, self.doi_journal]

    @property
    def meta_fn(self):
        if self.meta_path is None:
            meta = self.data_dir
        else:
            meta = os.path.join(self.meta_path, self.directory_formatted_doi)
        if meta is None:
            return None
        else:
            return os.path.join(meta, "meta", self.id + "_meta.csv")

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
            fn_meta = os.path.join(dir_out, self.id + "_meta.csv")
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
        self.id_base = f"{clean_id_str(self._organism)}_" \
                       f"{clean_id_str(self._organ)}_" \
                       f"{self._year}_" \
                       f"{clean_id_str(self._assay_sc)}_" \
                       f"{clean_id_str(author)}_" \
                       f"{idx}"

    @property
    def additional_annotation_key(self) -> Union[None, str]:
        return self._additional_annotation_key

    @additional_annotation_key.setter
    def additional_annotation_key(self, x: str):
        self._additional_annotation_key = x

    @property
    def annotated(self) -> Union[bool, None]:
        return self.cell_type_obs_key is not None

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
    def id(self) -> str:
        """
        Extends base ID by directory formatted (cleaned) DOI.

        In contrast to the based ID, this ID is guaranteed to be unique across studies and is also suitable as a file
        name.
        """
        return self.id_base + "_" + clean_doi(self.doi_main)

    @property
    def id_base(self) -> str:
        """
        Base ID of a dataset, unique identifies dataset within study.

        See also .id().
        """
        if self._id_base is not None:
            return self._id_base
        else:
            raise AttributeError(f"Dataset ID was not set in dataloader in {self.doi_main}, please ensure the "
                                 f"dataloader constructor of this dataset contains a call to self.set_dataset_id()")

    @id_base.setter
    def id_base(self, x: str):
        self._id_base = x

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
    def ncells(self) -> Union[None, int]:
        # ToDo cache this if it was loaded from meta?
        if self.adata is not None:
            return self.adata.n_obs
        else:
            return None

    @property
    def title(self):
        if self._title is None:
            x = crossref_query(doi=self.doi_main, k="title")
            # This is a list of length 1 if a title is found:
            if isinstance(x, list):
                x = x[0]
            self._title = x
        return self._title

    @property
    def celltypes_universe(self):
        if self._celltype_universe is None:
            self._celltype_universe = CelltypeUniverse(
                cl=getattr(self.ontology_container_sfaira, "cell_type"),
                uberon=getattr(self.ontology_container_sfaira, "organ"),
            )
        return self._celltype_universe

    @property
    def ontology_class_maps(self) -> Dict[str, pd.DataFrame]:
        return self._ontology_class_maps

    @ontology_class_maps.setter
    def ontology_class_maps(self, x: Dict[str, pd.DataFrame]):
        for k, v in x.items():
            assert v.shape[1] in [2, 3], f"{v.shape} in {self.id_base}"
            assert v.columns[0] == self._adata_ids.classmap_source_key
            assert v.columns[1] == self._adata_ids.classmap_target_key
            # Check for weird entries:
            # nan arises if columns was empty in that row.
            nan_vals = np.where([
                False if isinstance(x, str) else (np.isnan(x) or x is None)
                for x in v[self._adata_ids.classmap_target_key].values.tolist()
            ])[0]
            assert len(nan_vals) == 0, \
                f"Found nan target values in {self.id_base} for {x[self._adata_ids.classmap_target_key].values[nan_vals]}," \
                f" check if all entries in cell type .tsv file are non-empty, for example." \
                f"This bug may arise if a tab separator of columns is missing in one or multiple rows, for example."
            # Transform data frame into a mapping dictionary:
            self._ontology_class_maps[k] = dict(list(zip(
                v[self._adata_ids.classmap_source_key].values.tolist(),
                v[self._adata_ids.classmap_target_key].values.tolist()
            )))

    def subset(self, key, values, allow_missing_annotation: bool = False, allow_partial_match: bool = True) -> bool:
        """
        Subset list of adata objects based on cell-wise properties.

        These keys are properties that are not available in lazy model and require loading first because the
        subsetting works on the cell-level: .adata are maintained but reduced to matches.

        :param key: Property to subset by. Choose from sfaira controlled meta data:

                - "assay_sc"
                - "assay_differentiation"
                - "assay_type_differentiation
                - "bio_sample"
                - "cell_line"
                - "cell_type"
                - "development_stage"
                - "disease"
                - "ethnicity"
                - "id"
                - "individual"
                - "organ"
                - "organism"
                - "sex"
                - "state_exact"
                - "sample_source"
                - "tech_sample"

            The name key is independent of the schema used for streamlining uns and obs annotation.
        :param values: Target labels. All observations that are annotated with one of these labels or their descendants
            (for ontology-based meta data) are kept in the subset.
        :param allow_missing_annotation: Whether to add cells and data sets with missing annotation for queried key into
            the selection of the subset.
        :param allow_partial_match: Whether to allow all or no cells in a dataset that is not yet loaded into memory but
            that has (partially) matching cell-wise annotation for the given subset filter.
        :return: If the dataset should be kept for downstream operations.
            If the dataset was already loaded, this says if it was empty after sub-setting.
            If it was not loaded yet but fully satisfied the sub-setting, this is True.
            If it was not loaded yet but partially satisfied the sub-setting, this is allow_partial_match. Note in this
            case, we know that some cells in the data set fit this criterion, but because the data set is not yet loaded
            we do not know which ones. If you want to avoid this ambiguity, load the data before subset().
            If it was not loaded yet and did not satisfy the sub-setting, this is False.
        """
        if not hasattr(self._adata_ids, key):
            raise ValueError(f"attempted to subset based on non-supported key={key}, "
                             "choose from sfaira controlled meta data")
        if not (isinstance(values, list) or isinstance(values, tuple) or isinstance(values, np.ndarray)):
            values = [values]

        def get_subset_idx(samplewise_key, cellwise_key):
            sample_attr = getattr(self, samplewise_key)
            if sample_attr is not None and not (isinstance(sample_attr, list) or isinstance(sample_attr, tuple) or
                                                isinstance(sample_attr, np.ndarray)):
                sample_attr = [sample_attr]
            if sample_attr is not None:
                if len(sample_attr) == 1:
                    # Only use sample-wise subsetting if the sample-wise attribute is unique (not mixed).
                    if np.any([x in values for x in sample_attr]):
                        idx_ = np.arange(1, self.ncells)
                        keep_dataset_ = True
                    else:
                        idx_ = np.array([])
                        keep_dataset_ = False
                else:
                    # No per cell annotation and ambiguous sample annotation: pass entire data set if some keys match.
                    raise NotImplementedError(f"{self.id}: {(samplewise_key, cellwise_key)}")
            elif hasattr(self, cellwise_key) and getattr(self, cellwise_key) is not None:
                if self.loaded:
                    # Look for the streamlined ID column in .obs.
                    # Alternatively look for the streamlined symbol column in .obs-
                    col_ids = getattr(self._adata_ids, samplewise_key) + self._adata_ids.onto_id_suffix
                    col_symbols = getattr(self._adata_ids, samplewise_key)
                    # Check for availability:
                    col_ids = col_ids if col_ids in self.adata.obs.columns else None
                    col_symbols = col_symbols if col_symbols in self.adata.obs.columns else None
                    if col_ids in self.adata.obs.columns:
                        col_lookup = col_ids
                    elif col_symbols in self.adata.obs.columns:
                        col_lookup = col_symbols
                    else:
                        # TODO is this ever useful?
                        col_lookup = cellwise_key
                    values_found = self.adata.obs[col_lookup].values
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
                    idx_ = np.where([x in values_found_unique_matched for x in values_found])[0]
                    keep_dataset_ = len(idx_) > 0
                else:
                    if allow_partial_match:
                        # TODO find keys from tsvs here
                        idx_ = None
                        keep_dataset_ = False
                    else:
                        idx_ = None
                        keep_dataset_ = False
            else:
                if allow_missing_annotation:
                    # Pass all cells in object:
                    idx_ = np.arange(0, self.adata.n_obs)
                    keep_dataset_ = True
                else:
                    # Pass none of the cells in object:
                    idx_ = np.array([])
                    keep_dataset_ = False
            return idx_, keep_dataset_

        idx, keep_dataset = get_subset_idx(samplewise_key=key, cellwise_key=key + "_obs_key")
        if self.loaded:
            self.adata = self.adata[idx, :].copy() if len(idx) > 0 else None
        return keep_dataset

    def show_summary(self):
        print(f"{(self.supplier, self.organism, self.organ, self.assay_sc, self.disease)}")

    def __assert_loaded(self):
        if self.adata is None:
            raise ValueError("adata was not loaded, this is necessary for this operation")
