from __future__ import annotations

import anndata
import multiprocessing
import numpy as np
import pandas as pd
import os
from os import PathLike
import pandas
import pydoc
import scipy.sparse
from typing import Dict, List, Union
import warnings

from sfaira.data.dataloaders.base.dataset import DatasetBase
from sfaira.data.dataloaders.base.utils import is_child
from sfaira.versions.genomes.genomes import GenomeContainer
from sfaira.consts import AdataIds, AdataIdsSfaira
from sfaira.data.utils import read_yaml

UNS_STRING_META_IN_OBS = "__obs__"


def map_fn(inputs):
    """
    Functional to load data set with predefined additional actions.

    :param inputs:
    :return: None if function ran, error report otherwise
    """
    ds, load_raw, allow_caching, kwargs, func, kwargs_func = \
        inputs
    try:
        ds.load(load_raw=load_raw, allow_caching=allow_caching, **kwargs)
        if func is not None:
            x = func(ds, **kwargs_func)
            ds.clear()
            return x
        else:
            return None
    except (FileNotFoundError, OSError) as e:
        return ds.id, e,


def merge_uns_from_list(adata_ls):
    """
    Merge .uns from list of adata objects.

    Merges values for innert join of keys across all objects. This will retain uns streamlining.
    Keeps shared uns values for a given key across data sets as single value (not a list of 1 unique value).
    Other values are represented as a list of all unique values found.
    """
    uns_keys = [list(x.uns.keys()) for x in adata_ls]
    uns_keys_shared = set(uns_keys[0])
    for x in uns_keys[1:]:
        uns_keys_shared = uns_keys_shared.intersection(set(x))
    uns_keys_shared = list(uns_keys_shared)
    uns = {}
    for k in uns_keys_shared:
        uns_k = []
        for y in adata_ls:
            x = y.uns[k]
            if isinstance(x, list):
                pass
            elif isinstance(x, tuple):
                x = list(x)
            elif isinstance(x, np.ndarray):
                x = x.tolist()
            else:
                x = [x]
            uns_k.extend(x)
        uns_k = np.sort(np.unique(uns_k)).tolist()
        if len(uns_k) == 1:
            uns_k = uns_k[0]
        uns[k] = uns_k
    return uns


load_doc = \
    """
    :param remove_gene_version: Remove gene version string from ENSEMBL ID so that different versions in different data sets are superimposed.
    :param match_to_reference: Reference genomes name or False to keep original feature space.
    :param load_raw: Loads unprocessed version of data if available in data loader.
    :param allow_caching: Whether to allow method to cache adata object for faster re-loading.
    """


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
    _collection_id: str

    def __init__(self, datasets: dict, collection_id: str = "default"):
        self._adata_ids = AdataIdsSfaira()
        self.datasets = datasets
        self._collection_id = collection_id

    @property
    def _unknown_celltype_identifiers(self):
        return np.unqiue(np.concatenate([v._unknown_celltype_identifiers for _, v in self.datasets.items()]))

    def load(
            self,
            annotated_only: bool = False,
            load_raw: bool = False,
            allow_caching: bool = True,
            processes: int = 1,
            func=None,
            kwargs_func: Union[None, dict] = None,
            verbose: int = 0,
            **kwargs
    ):
        """
        Load all datasets in group (option for temporary loading).

        Note: This method automatically subsets to the group to the data sets for which input files were found.

        This method also allows temporarily loading data sets to execute function on loaded data sets (supply func).
        In this setting, datasets are removed from memory after the function has been executed.

        :param annotated_only:
        :param load_raw:
        :param allow_caching:
        :param processes: Processes to parallelise loading over. Uses python multiprocessing if > 1, for loop otherwise.
        :param func: Function to run on loaded datasets. map_fun should only take one argument, which is a Dataset
            instance. The return can be empty:

                def func(dataset, **kwargs_func):
                    # code manipulating dataset and generating output x.
                    return x
        :param kwargs_func: Kwargs of func.
        :param verbose: Verbosity of description of loading failure.

                - 0: no indication of failure
                - 1: indication of which data set failed in warning
                - 2: 1 with error report in warning
                - 3: reportin as in 2 but aborts with OSError
        """
        args = [
            load_raw,
            allow_caching,
            kwargs,
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
                    if verbose == 1:
                        warnings.warn(f"data set {k} not loaded")
                    if verbose == 2:
                        warnings.warn(f"data set {k} not loaded\nin data set {x[0]}: {x[1]}")
                    if verbose == 3:
                        raise OSError(f"data set {k} not loaded\nin data set {x[0]}: {x[1]}")
                    datasets_to_remove.append(k)
            for k in datasets_to_remove:
                del self.datasets[k]

    load.__doc__ += load_doc

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
        Streamline the adata instance in each data set to output format.
        Output format are saved in ADATA_FIELDS* classes.

        :param schema: Export format.
            - "sfaira"
            - "cellxgene"
        :param clean_obs: Whether to delete non-streamlined fields in .obs, .obsm and .obsp.
        :param clean_var: Whether to delete non-streamlined fields in .var, .varm and .varp.
        :param clean_uns: Whether to delete non-streamlined fields in .uns.
        :param clean_obs_names: Whether to replace obs_names with a string comprised of dataset id and an increasing integer.
        :param clean_obs_names: Whether to replace obs_names with a string comprised of dataset id and an increasing
            integer.
        :param keep_orginal_obs: For ontology-constrained .obs columns, whether to keep a column with original
            annotation.
        :param keep_symbol_obs: For ontology-constrained .obs columns, whether to keep a column with ontology symbol
            annotation.
        :param keep_id_obs: For ontology-constrained .obs columns, whether to keep a column with ontology ID annotation.
        :return:
        """
        for x in self.ids:
            self.datasets[x].streamline_metadata(
                schema=schema,
                clean_obs=clean_obs,
                clean_var=clean_var,
                clean_uns=clean_uns,
                clean_obs_names=clean_obs_names,
                keep_orginal_obs=keep_orginal_obs,
                keep_symbol_obs=keep_symbol_obs,
                keep_id_obs=keep_id_obs,
            )

    def streamline_features(
            self,
            match_to_reference: Union[str, Dict[str, str], None],
            remove_gene_version: bool = True,
            subset_genes_to_type: Union[None, str, List[str]] = None,
    ):
        """
        Subset and sort genes to genes defined in an assembly or genes of a particular type, such as protein coding.
        :param match_to_reference: Which annotation to map the feature space to. Can be:
            - str: Provide the name of the annotation in the format Organism.Assembly.Release
            - dict: Mapping of organism to name of the annotation (see str format). Chooses annotation for each data set
                based on organism annotation.
        :param remove_gene_version: Whether to remove the version number after the colon sometimes found in ensembl gene ids.
        :param subset_genes_to_type: Type(s) to subset to. Can be a single type or a list of types or None. Types can be:
            - None: All genes in assembly.
            - "protein_coding": All protein coding genes in assembly.
        """
        for x in self.ids:
            self.datasets[x].streamline_features(
                match_to_reference=match_to_reference,
                remove_gene_version=remove_gene_version,
                subset_genes_to_type=subset_genes_to_type,
            )

    def collapse_counts(self):
        """
        Collapse count matrix along duplicated index.
        """
        for x in self.ids:
            self.datasets[x].collapse_counts()

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
        This method writes a separate file for each data set in this object.

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
        for _, v in self.datasets.items():
            v.write_distributed_store(dir_cache=dir_cache, store_format=store_format, dense=dense,
                                      compression_kwargs=compression_kwargs, chunks=chunks)

    def write_backed(
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
                    self.datasets[x].write_backed(
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

        :param fn: File name of tsv to write class maps to.
        :param protected_writing: Only write if file was not already found.
        """
        tab = []
        for k, v in self.datasets.items():
            if v.annotated:
                labels_original = np.sort(np.unique(np.concatenate([
                    v.adata.obs[v.cell_type_original_obs_key].values
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
            tab = tab.sort_values(self._adata_ids.classmap_source_key)
            if not os.path.exists(fn) or not protected_writing:
                tab.to_csv(fn, index=False, sep="\t")

    def download(self, **kwargs):
        for _, v in self.datasets.items():
            v.download(**kwargs)

    @property
    def ids(self):
        return list(self.datasets.keys())

    @property
    def collection_id(self):
        return self._collection_id

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
        adata_ls = self.adata_ls
        if not adata_ls:
            return None
        if len(adata_ls) == 1:
            for i in self.ids:
                if self.datasets[i] is not None:
                    if self.datasets[i].adata is not None:
                        ds_id = i
            adata_concat = adata_ls[0]
            adata_concat.obs[self._adata_ids.dataset] = ds_id
        else:
            # Check that all individual adata objects in linked Dataset instances have identicall streamlined features and metadata
            match_ref_list = []
            rm_gene_ver_list = []
            gene_type_list = []
            for d_id in self.ids:
                if self.datasets[d_id].adata is not None:
                    assert self.datasets[d_id].mapped_features, f"Dataset {d_id} does not seem to have a streamlined " \
                                                                f"featurespace. To obtain an adata object from this " \
                                                                f"DatasetGroup, all contained Datasets need to have a " \
                                                                f"streamlined featurespace. Run .streamline_features()" \
                                                                f" first."
                    assert self.datasets[d_id].streamlined_meta, f"Dataset {d_id} does not seem to have streamlined " \
                                                                 f"metadata. To obtain an adata object from this " \
                                                                 f"DatasetGroup, all contained Datasets need to have " \
                                                                 f"streamlined metadata. Run .streamline_metadata() first."
                    match_ref_list.append(self.datasets[d_id].mapped_features)
                    rm_gene_ver_list.append(self.datasets[d_id].remove_gene_version)
                    gene_type_list.append(self.datasets[d_id].subset_gene_type)
            assert len(set(match_ref_list)) == 1, \
                "Not all datasets in this group had their features matched to the same reference (argument " \
                "'match_to_reference' of method .streamline_features())." \
                "This is however a prerequisite for creating a combined adata object."
            assert len(set(rm_gene_ver_list)) == 1, \
                "Not all datasets in this group have had their gene version removed (argument 'remove_gene_version' of " \
                "method .streamline_features()). This is however a prerequisite for creating a combined adata object."
            assert len(set(gene_type_list)) == 1, \
                "Not all datasets in this group had their featurespace subsetted to the same gene type (argument " \
                "'subset_gene_type' of method .streamline_features()). This is however a prerequisite for creating a " \
                "combined adata object."

            var_original = adata_ls[0].var.copy()
            for a in adata_ls:
                a.var_names_make_unique()
            # TODO: need to keep this? -> yes, still catching errors here (March 2020)
            # Fix for loading bug: sometime concatenating sparse matrices fails the first time but works on second try.
            try:
                adata_concat = adata_ls[0].concatenate(
                    *adata_ls[1:],
                    join="outer",
                    batch_key=self._adata_ids.dataset,
                    batch_categories=[i for i in self.ids if self.datasets[i].adata is not None],
                    index_unique=None
                )
            except ValueError:
                adata_concat = adata_ls[0].concatenate(
                    *adata_ls[1:],
                    join="outer",
                    batch_key=self._adata_ids.dataset,
                    batch_categories=[i for i in self.ids if self.datasets[i].adata is not None],
                    index_unique=None
                )
            adata_concat.var = var_original
            adata_concat.uns = merge_uns_from_list(adata_ls)
            adata_concat.uns[self._adata_ids.mapped_features] = match_ref_list[0]

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
            ] + [(self._adata_ids.dataset, [x for _ in range(self.datasets[x].adata.obs.shape[0])])]
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

    def project_celltypes_to_ontology(self, adata_fields: Union[AdataIds, None] = None, copy=False):
        """
        Project free text cell type names to ontology based on mapping table.
        :return:
        """
        for _, v in self.datasets.items():
            v.project_free_to_ontology(adata_fields=adata_fields, copy=copy)

    def subset(self, key, values: Union[list, tuple, np.ndarray]):
        """
        Subset list of adata objects based on sample-wise properties.

        These keys are properties that are available in lazy model.
        Subsetting happens on .datasets.

        :param key: Property to subset by.
        :param values: Classes to overlap to. Return if elements match any of these classes.
        :return:
        """
        ids_del = []
        if isinstance(values, np.ndarray):
            values = values.tolist()
        if isinstance(values, tuple):
            values = list(values)
        if not isinstance(values, list):
            values = [values]
        for x in self.ids:
            try:
                values_found = getattr(self.datasets[x], key)
            except AttributeError:
                raise ValueError(f"{key} not a valid property of data set object")
            try:
                ontology = getattr(self.datasets[x].ontology_container_sfaira, key)
            except AttributeError:
                raise ValueError(f"{key} not a valid property of ontology_container object")
            if values_found is None:
                # Delete entries which do not have this meta data item annotated.
                ids_del.append(x)
            else:
                if not isinstance(values_found, list):
                    values_found = [values_found]
                if not np.any([
                    np.any([
                        is_child(query=y, ontology=ontology, ontology_parent=z)
                        for z in values
                    ]) for y in values_found
                ]):
                    # Delete entries which a non-matching meta data value associated with this item.
                    ids_del.append(x)
        for x in ids_del:
            del self.datasets[x]

    def subset_cells(self, key, values: Union[str, List[str]]):
        """
        Subset list of adata objects based on cell-wise properties.

        These keys are properties that are not available in lazy model and require loading first because the
        subsetting works on the cell-level: .adata are maintained but reduced to matches.

        :param key: Property to subset by. Options:

            - "assay_differentiation" points to self.assay_differentiation_obs_key
            - "assay_sc" points to self.assay_sc_obs_key
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
        for x in self.ids:
            self.datasets[x].subset_cells(key=key, values=values)
            if self.datasets[x].ncells == 0:  # No observations (cells) left.
                del self.datasets[x]

    @property
    def additional_annotation_key(self) -> Dict[str, Union[None, str]]:
        """"
        Return dictionary of additional_annotation_key for each data set with ids as keys.
        """
        return dict([
            (k, self.datasets[k].additional_annotation_key)
            for k, v in self.datasets.items()
        ])

    @additional_annotation_key.setter
    def additional_annotation_key(self, x: Dict[str, Union[None, str]]):
        """
        Allows setting of additional_annotation_key in a subset of datasets identifed by keys in x.

        :param x: Dictionary with data set ids in keys and new _additional_annotation_key values to be setted in values.
            Note that you can either add  or change secondary annotation by setting a value to a string or remove it
            by setting a value to None.
        :return:
        """
        for k, v in x.items():
            self.datasets[k].additional_annotation_key = v

    @property
    def doi(self) -> List[str]:
        """
        Propagates DOI annotation from contained datasets.
        """
        dois = []
        for _, v in self.datasets.items():
            vdoi = v.doi_journal
            if isinstance(vdoi, str):
                vdoi = [vdoi]
            dois.extend(vdoi)
        return np.sort(np.unique(dois)).tolist()

    @property
    def supplier(self) -> List[str]:
        """
        Propagates supplier annotation from contained datasets.
        """
        supplier = [v.supplier for _, v in self.datasets.items()]
        return np.sort(np.unique(supplier)).tolist()

    def show_summary(self):
        for k, v in self.datasets.items():
            print(k)
            print(f"\t {(v.supplier, v.organism, v.organ, v.assay_sc, v.disease)}")


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
        try:
            collection_id = str(self._cwd).split(os.sep)[-1]
            package_source = str(self._cwd).split(os.sep)[-5]
            if package_source == "sfaira":
                pass
            elif package_source == "sfaira_extension":
                package_source = "sfairae"
            else:
                raise ValueError(f"invalid package source {package_source} for {self._cwd}, {self.collection_id}")
        except IndexError as e:
            raise IndexError(f"{e} for {self._cwd}")
        loader_pydoc_path_sfaira = "sfaira.data.dataloaders.loaders."
        loader_pydoc_path_sfairae = "sfaira_extension.data.dataloaders.loaders."
        loader_pydoc_path = loader_pydoc_path_sfaira if package_source == "sfaira" else loader_pydoc_path_sfairae
        if "group.py" in os.listdir(self._cwd):
            DatasetGroupFound = pydoc.locate(loader_pydoc_path + collection_id + ".group.DatasetGroup")
            dsg = DatasetGroupFound(data_path=data_path, meta_path=meta_path, cache_path=cache_path)
            datasets.extend(list(dsg.datasets.values))
        else:
            for f in os.listdir(self._cwd):
                if os.path.isfile(os.path.join(self._cwd, f)):  # only files
                    # Narrow down to data set files:
                    if f.split(".")[-1] == "py" and f.split(".")[0] not in ["__init__", "base", "group"]:
                        datasets_f = []
                        file_module = ".".join(f.split(".")[:-1])
                        DatasetFound = pydoc.locate(loader_pydoc_path + collection_id + "." + file_module + ".Dataset")
                        # Load objects from name space:
                        # - load(): Loading function that return anndata instance.
                        # - SAMPLE_FNS: File name list for DatasetBaseGroupLoadingManyFiles
                        load_func = pydoc.locate(loader_pydoc_path + collection_id + "." + file_module + ".load")
                        load_func_annotation = \
                            pydoc.locate(loader_pydoc_path + collection_id + "." + file_module + ".LOAD_ANNOTATION")
                        # Also check sfaira_extension for additional load_func_annotation:
                        if package_source != "sfairae":
                            load_func_annotation_sfairae = pydoc.locate(loader_pydoc_path_sfairae + collection_id +
                                                                        "." + file_module + ".LOAD_ANNOTATION")
                            # LOAD_ANNOTATION is a dictionary so we can use update to extend it.
                            if load_func_annotation_sfairae is not None and load_func_annotation is not None:
                                load_func_annotation.update(load_func_annotation_sfairae)
                            elif load_func_annotation_sfairae is not None and load_func_annotation is None:
                                load_func_annotation = load_func_annotation_sfairae
                        sample_fns = pydoc.locate(loader_pydoc_path + collection_id + "." + file_module +
                                                  ".SAMPLE_FNS")
                        fn_yaml = os.path.join(self._cwd, file_module + ".yaml")
                        fn_yaml = fn_yaml if os.path.exists(fn_yaml) else None
                        # Check for sample_fns in yaml:
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
                                        dict_load_func_annotation=load_func_annotation,
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
                                        load_func_annotation=load_func_annotation,
                                        sample_fn=x,
                                        sample_fns=sample_fns if sample_fns != [None] else None,
                                        yaml_path=fn_yaml,
                                    )
                                )
                        # Load cell type maps:
                        for x in datasets_f:
                            x.read_ontology_class_map(fn=os.path.join(self._cwd, file_module + ".tsv"))
                        datasets.extend(datasets_f)

        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)), collection_id=collection_id)

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
                        tab = list(self.datasets.values())[0]._read_ontology_class_map(fn=fn_map)
                        # Checks that the assigned ontology class names appear in the ontology.
                        list(self.datasets.values())[0]._value_protection(
                            attr="cell_type",
                            allowed=self.ontology_celltypes,
                            attempted=[
                                x for x in np.unique(tab[self._adata_ids.classmap_target_key].values)
                                if x not in [
                                    self._adata_ids.unknown_metadata_identifier,
                                    self._adata_ids.not_a_cell_celltype_identifier
                                ]
                            ]
                        )
                        # Adds a third column with the corresponding ontology IDs into the file.
                        tab[self._adata_ids.classmap_target_id_key] = [
                            self.ontology_celltypes.convert_to_id(x)
                            if (x != self._adata_ids.unknown_metadata_identifier and
                                x != self._adata_ids.not_a_cell_celltype_identifier)
                            else self._adata_ids.unknown_metadata_identifier
                            for x in tab[self._adata_ids.classmap_target_key].values
                        ]
                        # Get writing function from any (first) data set instance:
                        k = list(self.datasets.keys())[0]
                        self.datasets[k]._write_ontology_class_map(fn=fn_map, tab=tab)


class DatasetSuperGroup:
    """
    Container for multiple DatasetGroup instances.

    Used to manipulate structured dataset collections. Primarly designed for this manipulation, convert to DatasetGroup
    via flatten() for more functionalities.
    """
    _adata: Union[None, anndata.AnnData]
    fn_backed: Union[None, PathLike]
    dataset_groups: Union[list, List[DatasetGroup], List[DatasetSuperGroup]]

    def __init__(self, dataset_groups: Union[None, List[DatasetGroup], List[DatasetSuperGroup]]):
        self._adata = None
        self.fn_backed = None
        self.set_dataset_groups(dataset_groups=dataset_groups)

        self._adata_ids = AdataIdsSfaira()

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

    @property
    def ids(self):
        ids = []
        for x in self.dataset_groups:
            ids.extend(x.ids)
        return ids

    def get_gc(
            self,
            genome: str = None
    ):
        if genome.lower().startswith("homo_sapiens") or genome.lower().startswith("mus_musculus"):
            g = GenomeContainer(
                assembly=genome
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

    @property
    def datasets(self) -> Dict[str, DatasetBase]:
        """
        Returns DatasetGroup (rather than self = DatasetSuperGroup) containing all listed data sets.

        :return:
        """
        ds = {}
        for x in self.dataset_groups:
            for k, v in x.datasets.items():
                assert k not in ds.keys(), f"{k} was duplicated in super group, remove duplicates before flattening"
                ds[k] = v
        return ds

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

    def load(
            self,
            annotated_only: bool = False,
            load_raw: bool = False,
            allow_caching: bool = True,
            processes: int = 1,
            **kwargs
    ):
        """
        Loads data set human into anndata object.

        :param annotated_only:
        :param load_raw: See .load().
        :param allow_caching: See .load().
        :param processes: Processes to parallelise loading over. Uses python multiprocessing if > 1, for loop otherwise.
            Note: parallelises loading of each dataset group, but not across groups.
        :return:
        """
        for x in self.dataset_groups:
            x.load(
                annotated_only=annotated_only,
                load_raw=load_raw,
                allow_caching=allow_caching,
                processes=processes,
                **kwargs
            )

    def streamline_features(
            self,
            match_to_reference: Union[str, Dict[str, str], None],
            remove_gene_version: bool = True,
            subset_genes_to_type: Union[None, str, List[str]] = None,
    ):
        """
        Subset and sort genes to genes defined in an assembly or genes of a particular type, such as protein coding.
        :param match_to_reference: Which annotation to map the feature space to. Can be:
            - str: Provide the name of the annotation in the format Organism.Assembly.Release
            - dict: Mapping of organism to name of the annotation (see str format). Chooses annotation for each data set
                based on organism annotation.
        :param remove_gene_version: Whether to remove the version number after the colon sometimes found in ensembl gene ids.
        :param subset_genes_to_type: Type(s) to subset to. Can be a single type or a list of types or None. Types can be:
            - None: All genes in assembly.
            - "protein_coding": All protein coding genes in assembly.
        """
        for x in self.dataset_groups:
            x.streamline_features(
                match_to_reference=match_to_reference,
                remove_gene_version=remove_gene_version,
                subset_genes_to_type=subset_genes_to_type,
            )

    def collapse_counts(self):
        """
        Collapse count matrix along duplicated index.
        """
        for x in self.dataset_groups:
            x.collapse_counts()

    @property
    def adata_ls(self):
        adata_ls = []
        for k, v in self.flatten().datasets.items():
            if v.adata is not None:
                adata_ls.append(v.adata)
        return adata_ls

    @property
    def adata(self):
        adata_ls = self.adata_ls
        if not adata_ls:
            return None
        if len(adata_ls) == 1:
            for i in self.ids:
                if self.datasets[i] is not None:
                    if self.datasets[i].adata is not None:
                        ds_id = i
            adata_concat = adata_ls[0]
            adata_concat.obs[self._adata_ids.dataset] = ds_id
        else:
            # Check that all individual adata objects in linked Dataset instances have identicall streamlined features and metadata
            match_ref_list = []
            rm_gene_ver_list = []
            gene_type_list = []
            for d_id in self.flatten().ids:
                if self.flatten().datasets[d_id].adata is not None:
                    assert self.flatten().datasets[d_id].mapped_features, f"Dataset {d_id} does not seem to have a streamlined " \
                                                                          f"featurespace. To obtain an adata object from this " \
                                                                          f"DatasetGroup, all contained Datasets need to have a " \
                                                                          f"streamlined featurespace. Run .streamline_features()" \
                                                                          f" first."
                    assert self.flatten().datasets[d_id].streamlined_meta, f"Dataset {d_id} does not seem to have streamlined " \
                                                                           f"metadata. To obtain an adata object from this " \
                                                                           f"DatasetGroup, all contained Datasets need to have " \
                                                                           f"streamlined metadata. Run .streamline_metadata() first."
                    match_ref_list.append(self.flatten().datasets[d_id].mapped_features)
                    rm_gene_ver_list.append(self.flatten().datasets[d_id].remove_gene_version)
                    gene_type_list.append(self.flatten().datasets[d_id].subset_gene_type)
            assert len(set(match_ref_list)) == 1, \
                "Not all datasets in this group had their features matched to the same reference (argument " \
                "'match_to_reference' of method .streamline_features()). This is however a prerequisite for creating a " \
                "combined adata object."
            assert len(set(rm_gene_ver_list)) == 1, \
                "Not all datasets in this group have had their gene version removed (argument 'remove_gene_version' of " \
                "method .streamline_features()). This is however a prerequisite for creating a combined adata object."
            assert len(set(gene_type_list)) == 1, \
                "Not all datasets in this group had their featurespace subsetted to the same gene type (argument " \
                "'subset_gene_type' of method .streamline_features()). This is however a prerequisite for creating a " \
                "combined adata object."

            var_original = adata_ls[0].var.copy()
            for a in adata_ls:
                a.var_names_make_unique()
            # TODO: need to keep this? -> yes, still catching errors here (March 2020)
            # Fix for loading bug: sometime concatenating sparse matrices fails the first time but works on second try.
            try:
                adata_concat = adata_ls[0].concatenate(
                    *adata_ls[1:],
                    join="outer",
                    batch_key=self._adata_ids.dataset,
                    batch_categories=[i for i in self.ids if self.flatten().datasets[i].adata is not None],
                    index_unique=None
                )
            except ValueError:
                adata_concat = adata_ls[0].concatenate(
                    *adata_ls[1:],
                    join="outer",
                    batch_key=self._adata_ids.dataset,
                    batch_categories=[i for i in self.ids if self.flatten().datasets[i].adata is not None],
                    index_unique=None
                )
            adata_concat.var = var_original
            adata_concat.uns = merge_uns_from_list(adata_ls)
            adata_concat.uns[self._adata_ids.mapped_features] = match_ref_list[0]

        return adata_concat

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
        This method writes a separate file for each data set in this object.

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
        for x in self.dataset_groups:
            x.write_distributed_store(dir_cache=dir_cache, store_format=store_format, dense=dense,
                                      compression_kwargs=compression_kwargs, chunks=chunks)

    def write_backed(
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

        TODO replace streamlining in here by required call to .streamline() before.

        Example usage:

            ds = DatasetSuperGroup([...])
            ds.write_backed(
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
        n_genes = gc.n_var
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
            self._adata_ids.annotated,
            self._adata_ids.assay_sc,
            self._adata_ids.assay_differentiation,
            self._adata_ids.assay_type_differentiation,
            self._adata_ids.author,
            self._adata_ids.cell_line,
            self._adata_ids.dataset,
            self._adata_ids.cell_type,
            self._adata_ids.development_stage,
            self._adata_ids.normalization,
            self._adata_ids.organ,
            self._adata_ids.bio_sample,
            self._adata_ids.state_exact,
            self._adata_ids.year,
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
            x.write_backed(
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
        Streamline the adata instance in each group and each data set to output format.
        Output format are saved in ADATA_FIELDS* classes.

        :param schema: Export format.
            - "sfaira"
            - "cellxgene"
        :param clean_obs: Whether to delete non-streamlined fields in .obs, .obsm and .obsp.
        :param clean_var: Whether to delete non-streamlined fields in .var, .varm and .varp.
        :param clean_uns: Whether to delete non-streamlined fields in .uns.
        :param clean_obs_names: Whether to replace obs_names with a string comprised of dataset id and an increasing integer.
        :param clean_obs_names: Whether to replace obs_names with a string comprised of dataset id and an increasing
            integer.
        :param keep_orginal_obs: For ontology-constrained .obs columns, whether to keep a column with original
            annotation.
        :param keep_symbol_obs: For ontology-constrained .obs columns, whether to keep a column with ontology symbol
            annotation.
        :param keep_id_obs: For ontology-constrained .obs columns, whether to keep a column with ontology ID annotation.
        :return:
        """
        for x in self.dataset_groups:
            for xx in x.ids:
                x.datasets[xx].streamline_metadata(
                    schema=schema,
                    clean_obs=clean_obs,
                    clean_var=clean_var,
                    clean_uns=clean_uns,
                    clean_obs_names=clean_obs_names,
                    keep_orginal_obs=keep_orginal_obs,
                    keep_symbol_obs=keep_symbol_obs,
                    keep_id_obs=keep_id_obs,
                )

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

    def remove_duplicates(
            self,
            supplier_hierarchy: str = "cellxgene,sfaira"
    ):
        """
        Remove duplicate data loaders from super group, e.g. loaders that map to the same DOI.

        Any DOI match is removed (pre-print or journal publication).
        Data sets without DOI are removed, too.
        Loaders are kept in the hierarchy indicated in supplier_hierarchy.
        Requires a super group with homogenous suppliers across DatasetGroups, throws an error otherwise.
        This is given for sfaira maintained libraries but may not be the case if custom assembled DatasetGroups are
        used.

        :param supplier_hierarchy: Hierarchy to resolve duplications by.
            Comma separated string that indicates which data provider takes priority.
            Choose "cellxgene,sfaira" to prioritise use of data sets downloaded from cellxgene.
            Choose "sfaira,cellxgene" to prioritise use of raw data processing pipelines locally.

                - cellxgene: cellxgene downloads
                - sfaira: local raw file processing
        :return:
        """
        # Build a pairing of provider and DOI:
        report_list = []
        idx_tokeep = []
        supplier_hierarchy = supplier_hierarchy.split(",")
        for i, (x, y) in enumerate([(xx.supplier, xx.doi) for xx in self.dataset_groups]):
            if len(x) > 1:
                raise ValueError(f"found more than one supplier for DOI {str(y)}")
            else:
                x = x[0]
                if x not in supplier_hierarchy:
                    raise ValueError(f"could not associate supplier {x} with hierarchy {supplier_hierarchy} in "
                                     f"data set {y}")
                if len(report_list) > 0:
                    matched_idx = np.where([
                        np.any([
                            zz in y
                            for zz in z[1]
                        ])
                        for z in report_list
                    ])[0]
                    assert len(matched_idx) < 1, f"more matches than expected for {(x, y)}"
                else:
                    matched_idx = []
                if len(matched_idx) > 0:
                    # Establish which entry takes priority
                    supplier_old = report_list[matched_idx[0]][0]
                    priority = supplier_hierarchy.index(supplier_old) > supplier_hierarchy.index(x)
                    print(f"removing duplicate data set {y} from supplier: {supplier_old if priority else x}")
                    if priority:
                        idx_tokeep.append(i)
                        del idx_tokeep[matched_idx[0]]
                else:
                    report_list.append([x, y])
                    idx_tokeep.append(i)
        self.dataset_groups = [self.dataset_groups[i] for i in idx_tokeep]

    def subset_cells(self, key, values: Union[str, List[str]]):
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
        for i in range(len(self.dataset_groups)):
            self.dataset_groups[i].subset_cells(key=key, values=values)

    def project_celltypes_to_ontology(self, adata_fields: Union[AdataIds, None] = None, copy=False):
        """
        Project free text cell type names to ontology based on mapping table.
        :return:
        """
        for _, v in self.dataset_groups:
            v.project_free_to_ontology(adata_fields=adata_fields, copy=copy)

    def write_config(self, fn: Union[str, os.PathLike]):
        """
        Writes a config file that describes the current data sub-setting.

        This config file can be loaded later to recreate a sub-setting.

        :param fn: Output file.
        """
        pd.DataFrame({"id": np.sort(self.ids)}).to_csv(fn, index=False, sep="\t")

    def load_config(self, fn: Union[str, os.PathLike]):
        """
        Load a config file and recreates a data sub-setting.

        :param fn: Output file.
        """
        tab = pd.read_csv(fn, header=0, index_col=None, sep="\t")
        ids_keep = tab["id"].values
        self.subset(key="id", values=ids_keep)

    @property
    def additional_annotation_key(self) -> List[Dict[str, Union[None, str]]]:
        """"
        Return list (by data set group) of dictionaries of additional_annotation_key for each data set with ids as keys.
        """
        return [
            dict([
                (k, x.datasets[k].additional_annotation_key)
                for k, v in x.datasets.items()
            ]) for x in self.dataset_groups
        ]

    @additional_annotation_key.setter
    def additional_annotation_key(self, x: Dict[str, Union[None, str]]):
        """
        Allows setting of additional_annotation_key in a subset of datasets identifed by keys in x.

        The input is not structured by DatasetGroups but only by ID, all groups are checked for matching IDs.

        :param x: Dictionary with data set ids in keys and new _additional_annotation_key values to be setted in values.
            Note that you can either add  or change secondary annotation by setting a value to a string or remove it
            by setting a value to None.
        :return:
        """
        for k, v in x.items():
            counter = 0
            for x in self.dataset_groups:
                if k in x.ids:
                    x.datasets[k].additional_annotation_key = v
                    counter += 1
            if counter == 0:
                warnings.warn(f"did not data set matching ID {k}")
            elif counter > 1:
                warnings.warn(f"found more than one ({counter}) data set matching ID {k}")

    def show_summary(self):
        for k, v in self.datasets.items():
            print(k)
            print(f"\t {(v.supplier, v.organism, v.organ, v.assay_sc, v.disease)}")
