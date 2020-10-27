import abc
import anndata
import h5py
import numpy as np
import pandas as pd
import os
from os import PathLike
import pandas
import scipy.sparse
from typing import Dict, List, Union
import warnings

from .external import SuperGenomeContainer
from .external import ADATA_IDS


class DatasetBase(abc.ABC):

    adata: Union[None, anndata.AnnData]
    class_maps: dict
    meta: Union[None, pandas.DataFrame]
    download_website: Union[None, str]
    download_website_meta: Union[None, str]
    path: Union[None, str]
    id: Union[None, str]
    download_website: Union[None, str]
    organ: Union[None, str]
    sub_tissue: Union[None, str]
    has_celltypes: Union[None, bool]
    species: Union[None, str]
    genome: Union[None, str]

    def __init__(
            self,
            path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            **kwargs
    ):
        self.species = None
        self.adata = None
        self.download_website_meta = None
        self.id = None
        self.download_website = None
        self.organ = None
        self.sub_tissue = None
        self.has_celltypes = None
        self.meta = None
        self.genome = None
        self.path = path
        self.meta_path = meta_path
        self._load_raw = None


    @abc.abstractmethod
    def _load(self, fn):
        pass

    def load(
            self,
            celltype_version: Union[str, None] = None,
            fn: Union[str, None] = None,
            remove_gene_version: bool = True,
            match_to_reference: Union[str, None] = None,
            load_raw: bool = False
    ):
        """

        :param celltype_version: Version of cell type ontology to use. Uses most recent if None.
        :param fn: Optional target file name, otherwise infers from defined directory structure.
        :param remove_gene_version: Remove gene version string from ENSEMBL ID so that different versions in different
            data sets are superimposed.
        :param match_to_reference: Reference genomes name.
        :param load_raw: Loads unprocessed version of data if available in data loader.
        :return:
        """

        self._load_raw = load_raw

        if match_to_reference and not remove_gene_version:
            warnings.warn("it is not recommended to enable matching the feature space to a genomes reference"
                          "while not removing gene versions. this can lead to very poor matching performance")

        # set default genomes if none provided
        if match_to_reference:
            genome = match_to_reference
            self._set_genome(genome=genome)
        elif self.species == "human":
            genome = "Homo_sapiens_GRCh38_97"
            warnings.warn(f"using default genomes {genome}")
            self._set_genome(genome=genome)
        elif self.species == "mouse":
            genome = "Mus_musculus_GRCm38_97"
            warnings.warn(f"using default genomes {genome}")
            self._set_genome(genome=genome)

        self._load(fn=fn)

        if ADATA_IDS.cell_ontology_id not in self.adata.obs.columns:
            self.adata.obs[ADATA_IDS.cell_ontology_id] = None

        # Map cell type names from raw IDs to ontology maintained ones::
        self.adata.obs[ADATA_IDS.cell_ontology_class] = self.map_ontology_class(
            raw_ids=self.adata.obs[ADATA_IDS.cell_ontology_class].values,
            celltype_version=celltype_version
        )

        # Remove version tag on ensembl gene ID so that different versions are superimposed downstream:
        if remove_gene_version:
            new_index = [x.split(".")[0] for x in self.adata.var_names.tolist()]
            # Collapse if necessary:
            new_index_collapsed = list(np.unique(new_index))
            if len(new_index_collapsed) < self.adata.n_vars:
                raise ValueError("duplicate features detected after removing gene versions."
                                 "the code to collapse these features is implemented but not tested.")
                idx_map = np.array([new_index_collapsed.index(x) for x in new_index])
                # Need reverse sorting to find index of last element in sorted list to split array using list index().
                idx_map_sorted_rev = np.argsort(idx_map)[::-1]
                n_genes = len(idx_map_sorted_rev)
                # 1. Sort array in non-reversed order: idx_map_sorted_rev[::-1]
                # 2. Split into chunks based on blocks of identical entries in idx_map, using the occurrence of the
                # last element of each block as block boundaries:
                # n_genes - 1 - idx_map_sorted_rev.index(x)
                # Note that the blocks are named as positive integers starting at 1, without gaps.
                counts = np.concatenate([
                    np.sum(x, axis=1, keepdims=True)
                    for x in np.split(
                        self.adata[:, idx_map_sorted_rev[::-1]].X,  # forward ordered data
                        indices_or_sections=[
                            n_genes - 1 - idx_map_sorted_rev.index(x)  # last occurrence of element in forward order
                            for x in np.arange(0, len(new_index_collapsed)-1)  # -1: do not need end of last partition
                        ],
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
            self.adata.var[ADATA_IDS.gene_id_ensembl] = new_index
            self.adata.var.index = self.adata.var[ADATA_IDS.gene_id_ensembl].values

        # Match feature space to a genomes provided with sfaira
        if match_to_reference:
            # Convert data matrix to csc matrix
            if isinstance(self.adata.X, np.ndarray):
                # Change NaN to zero. This occurs for example in concatenation of anndata instances.
                if np.any(np.isnan(self.adata.X)):
                    self.adata.X[np.isnan(self.adata.X)] = 0
                x = scipy.sparse.csc_matrix(self.adata.X)
            elif isinstance(self.adata.X, scipy.sparse.spmatrix):
                x = self.adata.X.tocsc()
            else:
                raise ValueError("data type %s not recognized" % type(self.adata.X))

            # Compute indices of genes to keep
            data_ids = self.adata.var[ADATA_IDS.gene_id_ensembl].values
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
                                           ADATA_IDS.gene_id_ensembl: self.genome_container.ensembl},
                                     index=self.genome_container.ensembl),
                    uns=self.adata.uns
            )

        self.adata.uns['mapped_features'] = match_to_reference

    def _convert_and_set_var_names(
            self,
            symbol_col: str = None,
            ensembl_col: str = None,
            new_index: str = ADATA_IDS.gene_id_ensembl
    ):
        if symbol_col and ensembl_col:
            if symbol_col == 'index':
                self.adata.var.index.name = 'index'
                self.adata.var = self.adata.var.reset_index().rename(
                    {'index': ADATA_IDS.gene_id_names},
                    axis='columns'
                )
            else:
                self.adata.var = self.adata.var.rename(
                    {symbol_col:  ADATA_IDS.gene_id_names},
                    axis='columns'
                )

            if ensembl_col == 'index':
                self.adata.var.index.name = 'index'
                self.adata.var = self.adata.var.reset_index().rename(
                    {'index': ADATA_IDS.gene_id_ensembl},
                    axis='columns'
                )
            else:
                self.adata.var = self.adata.var.rename(
                    {ensembl_col: ADATA_IDS.gene_id_ensembl},
                    axis='columns'
                )

        elif symbol_col:
            id_dict = self.genome_container.names_to_id_dict
            id_strip_dict = self.genome_container.strippednames_to_id_dict
            if symbol_col == 'index':
                self.adata.var.index.name = 'index'
                self.adata.var = self.adata.var.reset_index().rename(
                    {'index':  ADATA_IDS.gene_id_names},
                    axis='columns'
                )
            else:
                self.adata.var = self.adata.var.rename(
                    {symbol_col:  ADATA_IDS.gene_id_names},
                    axis='columns'
                )

            # Matching gene names to ensembl ids in the following way: if the gene is present in the ensembl dictionary,
            # match it straight away, if it is not in there we try to match everything in front of the first period in
            # the gene name with a dictionary that was modified in the same way, if there is still no match we append na
            ensids = []
            for n in self.adata.var[ADATA_IDS.gene_id_names]:
                if n in id_dict.keys():
                    ensids.append(id_dict[n])
                elif n.split(".")[0] in id_strip_dict.keys():
                    ensids.append(id_strip_dict[n.split(".")[0]])
                else:
                    ensids.append('n/a')
            self.adata.var[ADATA_IDS.gene_id_ensembl] = ensids

        elif ensembl_col:
            id_dict = self.genome_container.id_to_names_dict
            if ensembl_col == 'index':
                self.adata.var.index.name = 'index'
                self.adata.var = self.adata.var.reset_index().rename(
                    {'index': ADATA_IDS.gene_id_ensembl}, 
                    axis='columns'
                )
            else:
                self.adata.var = self.adata.var.rename(
                    {ensembl_col: ADATA_IDS.gene_id_names}, 
                    axis='columns'
                )

            self.adata.var[ADATA_IDS.gene_id_names] = [
                id_dict[n.split(".")[0]] if n.split(".")[0] in id_dict.keys() else 'n/a'
                for n in self.adata.var[ADATA_IDS.gene_id_ensembl]
            ]

        else:
            raise ValueError('Please provide the name of at least the name of the var column containing ensembl ids or'
                             'the name of the var column containing gene symbols')

        self.adata.var.index = self.adata.var[new_index].tolist()
        self.adata.var_names_make_unique()

    def subset_organs(self, subset: Union[None, List]):
        if self.organ == "mixed":
            self.organsubset = subset
        else:
            raise ValueError("Only data that contain multiple organs can be subset.")
        if self.adata is not None:
            warnings.warn("You are trying to subset organs after loading the dataset."
                          "This will have no effect unless the dataset is loaded again.")

    def load_tobacked(
            self,
            adata_backed: anndata.AnnData,
            genome: str,
            idx: np.ndarray,
            keys: List[str] = [],
            fn: Union[None, str] = None,
            celltype_version: Union[str, None] = None,
    ):
        """
        Loads data set into slice of backed anndata object.

        Note: scatter updates to backed sparse arrays are not yet supported by anndata. Accordingly, we need to work
        around below using .append() of the backed matrix.

        :param adata_backed:
        :param genome: Genome name to use as refernce.
        :param idx: Indices in adata_backed to write observations to. This can be used to immediately create a
            shuffled object.
        :param keys:
        :param fn:
        :param celltype_version: Version of cell type ontology to use. Uses most recent if None.
        :return: New row index for next element to be written into backed anndata.
        """
        self.load(
            fn=fn,
            celltype_version=celltype_version,
            remove_gene_version=True,
            match_to_reference=genome
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
                if k == ADATA_IDS.dataset:
                    adata_backed.obs.loc[np.sort(idx), ADATA_IDS.dataset] = [self.id for i in range(len(idx))]
                elif k in self.adata.obs.columns:
                    adata_backed.obs.loc[np.sort(idx), k] = self.adata.obs[k].values[np.argsort(idx)]
                elif k in list(self.adata.uns.keys()):
                    adata_backed.obs.loc[np.sort(idx), k] = [self.adata.uns[k] for i in range(len(idx))]
                else:
                    # Need to fill this instead of throwing an exception as this condition can trigger for one element
                    # within a loop over multiple data sets (ie in data set groups).
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
                    (k, [self.id for i in range(len(idx))]) if k == ADATA_IDS.dataset
                    else (k, self.adata.obs[k].values[np.argsort(idx)]) if k in self.adata.obs.columns
                    else (k, [self.adata.uns[k] for i in range(len(idx))]) if k in list(self.adata.uns.keys())
                    else (k, ["key_not_found" for i in range(len(idx))])
                    for k in adata_backed.obs.columns
                ]))
            )
        else:
             raise ValueError("did not reccognize backed anndata.X format %s" % type(adata_backed.X))

    def set_unkown_class_id(self, ids: list):
        """
        Sets list of custom identifiers of unkown cell types in adata.obs["cell_ontology_class"] to the target one.

        :param ids: IDs in adata.obs["cell_ontology_class"] to replace.
        :return:
        """
        target_id = "unknown"
        ontology_classes = [
            x if x not in ids else target_id
            for x in self.adata.obs[ADATA_IDS.cell_ontology_class].tolist()
        ]
        self.adata.obs[ADATA_IDS.cell_ontology_class] = ontology_classes

    def _set_genome(self,
                    genome: str
    ):

        if genome.lower().startswith("homo_sapiens"):
            g = SuperGenomeContainer(
                species="human",
                genome=genome
            )
        elif genome.lower().startswith("mus_musculus"):
            g = SuperGenomeContainer(
                species="mouse",
                genome=genome
            )
        else:
            raise ValueError("genomes %s not recognised. please provide valid genomes." % genome)

        self.genome_container = g

    @property
    def doi_cleaned_id(self):
        return "_".join(self.id.split("_")[:-1])

    def load_meta(self, fn: Union[PathLike, str]):
        if fn is None:
            if self.meta_path is None:
                raise ValueError("provide either fn in load or path in constructor")
            fn = os.path.join(self.meta_path, self.doi_cleaned_id + "_meta.csv")
        else:
            if isinstance(fn, str):
                fn = os.path.normpath(fn)
        self.meta = pandas.read_csv(fn)

    @property
    def ncells(self):
        if self.meta is None:
            self.load_meta(fn=None)
        return int(self.meta["ncells"])

    def write_meta(
            self,
            fn_meta: Union[None, str] = None,
            fn_data: Union[None, str] = None,
            dir_out: Union[None, str] = None,
    ):
        if fn_meta is None:
            if self.path is None and dir_out is None:
                raise ValueError("provide either fn in load or path in constructor")
            if dir_out is None:
                dir_out = self.meta_path
            fn_meta = os.path.join(dir_out, self.doi_cleaned_id + "_meta.csv")
        if self.adata is None:
            self.load(fn=fn_data, remove_gene_version=False, match_to_reference=None)
        meta = pandas.DataFrame({
            "ncells": self.adata.n_obs,
            "animal": self.adata.uns[ADATA_IDS.animal],
            "organ": self.adata.uns[ADATA_IDS.organ],
            "subtissue": self.adata.uns[ADATA_IDS.subtissue],
            "id": self.adata.uns[ADATA_IDS.id],
            "lab": self.adata.uns[ADATA_IDS.lab],
            "year": self.adata.uns[ADATA_IDS.year],
            "protocol": self.adata.uns[ADATA_IDS.protocol],
            "counts": self.adata.uns[ADATA_IDS.normalization] if ADATA_IDS.normalization in self.adata.uns.keys() else None,
            "has_celltypes": self.has_celltypes
        }, index=range(1))
        meta.to_csv(fn_meta)

    @property
    def available_type_versions(self):
        return np.array(list(self.class_maps.keys()))

    def set_default_type_version(self):
        """
        Choose most recent version.

        :return: Version key corresponding to most recent version.
        """
        return self.available_type_versions[np.argmax([int(x) for x in self.available_type_versions])]

    def assert_celltype_version_key(
            self,
            celltype_version
    ):
        if celltype_version not in self.available_type_versions:
            raise ValueError(
                "required celltype version %s not found. available are: %s" %
                (celltype_version, str(self.available_type_versions))
            )

    def map_ontology_class(
            self,
            raw_ids,
            celltype_version
    ):
        """

        :param raw_ids:
        :param class_maps:
        :param celltype_version: Version of cell type ontology to use. Uses most recent if None.
        :return:
        """
        if celltype_version is None:
            celltype_version = self.set_default_type_version()
        self.assert_celltype_version_key(celltype_version=celltype_version)
        return [
            self.class_maps[celltype_version][x] if x in self.class_maps[celltype_version].keys() else x
            for x in raw_ids
        ]


class DatasetGroupBase(abc.ABC):
    """

    Example:

    #query human lung
    #from sfaira.dev.data.human.lung import DatasetGroupLung as DatasetGroup
    #dsg_humanlung = DatasetGroupHuman(path='path/to/data')
    #dsg_humanlung.load_all(match_to_reference='Homo_sapiens_GRCh38_97')
    #dsg_humanlung[some_id]
    #dsg_humanlung.adata
    """
    datasets: Dict


    def subset_organs(self, subset: Union[None, List]):
        for i in self.ids:
            if self.datasets[i].organ == "mixed":
                self.datasets[i].subset_organs(subset)
            else:
                raise ValueError("Only data that contain multiple organs can be subset.")

    def load_all(
            self,
            celltype_version: Union[str, None] = None,
            annotated_only: bool = False,
            remove_gene_version: bool = True,
            match_to_reference: Union[str, None] = None,
            load_raw: bool = False
    ):
        """

        :param celltype_version: Version of cell type ontology to use. Uses most recent if None.
        :param annotated_only:
        :param remove_gene_version:
        :param match_to_reference:
        :param load_raw: Loads unprocessed version of data if available in data loader.
        :return:
        """
        for i in self.ids:
            if self.datasets[i].has_celltypes or not annotated_only:
                self.datasets[i].load(
                    celltype_version=self.format_type_version(celltype_version),
                    remove_gene_version=remove_gene_version,
                    match_to_reference=match_to_reference,
                    load_raw=load_raw
                )

    def load_all_tobacked(
            self,
            adata_backed: anndata.AnnData,
            genome: str,
            idx: List[np.ndarray],
            keys: List[str] = [],
            annotated_only: bool = False,
            celltype_version: Union[str, None] = None,
    ):
        """
        Loads data set group into slice of backed anndata object.

        :param adata_backed:
        :param genome: Genome container target genomes loaded.
        :param idx: Indices in adata_backed to write observations to. This can be used to immediately create a
            shuffled object. This has to be a list of the length of self.data, one index array for each dataset.
        :param keys:
        :param annotated_only:
        :param celltype_version: Version of cell type ontology to use. Uses most recent if None.
        :return: New row index for next element to be written into backed anndata.
        """
        keys_to_always_load = ["organ"]
        for x in keys_to_always_load:
            if x not in keys:
                keys.append(x)
        for i, id in enumerate(self.ids):
            # if this is for celltype prediction, only load the data with have celltype annotation
            if self.datasets[id].has_celltypes or not annotated_only:
                self.datasets[id].load_tobacked(
                    adata_backed=adata_backed,
                    genome=genome,
                    idx=idx[i],
                    keys=keys,
                    celltype_version=self.format_type_version(celltype_version)
                )

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
        adata_ls = self.adata_ls
        # Save uns attributes that are fixed for entire data set to .obs to retain during concatenation:
        for adata in adata_ls:
            adata.obs[ADATA_IDS.lab] = adata.uns[ADATA_IDS.lab]
            adata.obs[ADATA_IDS.year] = adata.uns[ADATA_IDS.year]
            adata.obs[ADATA_IDS.protocol] = adata.uns[ADATA_IDS.protocol]
            adata.obs[ADATA_IDS.subtissue] = adata.uns[ADATA_IDS.subtissue]
            if ADATA_IDS.normalization in adata.uns.keys():
                adata.obs[ADATA_IDS.normalization] = adata.uns[ADATA_IDS.normalization]
            if ADATA_IDS.dev_stage in adata.obs.columns:
                adata.obs[ADATA_IDS.dev_stage] = adata.uns[ADATA_IDS.dev_stage]
            adata.obs[ADATA_IDS.has_celltypes] = adata.uns[ADATA_IDS.has_celltypes]
        # Workaround related to anndata bugs:  # TODO remove this in future.
        for adata in adata_ls:
            # Fix 1:
            if adata.raw is not None:
                adata.raw._varm = None
            # Fix 2:
            if adata.uns is not None:
                keys_to_keep = [
                    'neighbors',
                    ADATA_IDS.lab,
                    ADATA_IDS.year,
                    ADATA_IDS.protocol,
                    ADATA_IDS.subtissue,
                    ADATA_IDS.normalization,
                    ADATA_IDS.dev_stage,
                    ADATA_IDS.has_celltypes,
                    "mapped_features"
                ]
                for k in list(adata.uns.keys()):
                    if k not in keys_to_keep:
                        del adata.uns[k]
            # Fix 3:
            if not isinstance(adata.X, scipy.sparse.csr_matrix):
                adata.X = scipy.sparse.csr_matrix(adata.X)
        # .var entries are renamed and copied upon concatenation.
        # To preserve gene names in .var, the target gene names are copied into var_names and are then copied
        # back into .var.
        for adata in adata_ls:
            adata.var.index = adata.var[ADATA_IDS.gene_id_ensembl].tolist()
        if len(adata_ls) > 1:
            # TODO: need to keep this? -> yes, still catching errors here (March 2020)
            # Fix for loading bug: sometime concatenating sparse matrices fails the first time but works on second try.
            try:
                adata_concat = adata_ls[0].concatenate(
                    *adata_ls[1:],
                    join="outer",
                    batch_key=ADATA_IDS.dataset,
                    batch_categories=[i for i in self.ids if self.datasets[i].adata is not None]
                )
            except ValueError:
                adata_concat = adata_ls[0].concatenate(
                    *adata_ls[1:],
                    join="outer",
                    batch_key=ADATA_IDS.dataset,
                    batch_categories=[i for i in self.ids if self.datasets[i].adata is not None]
                )

            adata_concat.var[ADATA_IDS.gene_id_ensembl] = adata_concat.var.index

            if len(set([a.uns['mapped_features'] for a in adata_ls])) == 1:
                adata_concat.uns['mapped_features'] = adata_ls[0].uns['mapped_features']
            else:
                adata_concat.uns['mapped_features'] = False
        else:
            adata_concat = adata_ls[0]
            adata_concat.obs[ADATA_IDS.dataset] = self.ids[0]

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
                else (k, ["nan" for i in range(self.datasets[x].adata.obs.shape[0])])
                for k in keys
            ] + [(ADATA_IDS.dataset, [x for i in range(self.datasets[x].adata.obs.shape[0])])]
        )) for x in self.ids if self.datasets[x].adata is not None])
        return obs_concat

    @property
    def ncells(self):
        return sum([self.datasets[i].ncells for i in self.ids])

    @property
    def ncells_bydataset(self):
        return [self.datasets[i].ncells for i in self.ids]

    def assert_celltype_version_key(
            self,
            celltype_version
    ):
        """
        Assert that version key exists in each data set.
        :param celltype_version:
        :return:
        """
        for x in self.ids:
            if celltype_version not in self.datasets[x].available_type_versions:
                raise ValueError(
                    "required celltype version %s not found in data set %s. available are: %s" %
                    (celltype_version, x, str(self.datasets[x].available_type_versions))
                )

    def format_type_version(self, version):
        """
        Choose most recent version available in each dataset if None, otherwise return input version after checking.

        :return: Version key corresponding to default version.
        """
        if version is None:
            versions = set(self.datasets[self.ids[0]].available_type_versions)
            for x in self.ids[1:]:
                versions = versions.intersection(set(self.datasets[x].available_type_versions))
            versions = np.array(list(versions))
            return versions[np.argmax([int(x) for x in versions])]
        else:
            self.assert_celltype_version_key()
            return version


class DatasetSuperGroup:
    """
    Container for multiple DatasetGroup instances.

    Can be used to grid_searches models across organs. Supports backed anndata objects.
    """
    adata: Union[None, anndata.AnnData]
    fn_backed: Union[None, PathLike]
    dataset_groups: List[DatasetGroupBase]

    def __init__(self, dataset_groups: Union[None, List[DatasetGroupBase]]):
        self.adata = None
        self.fn_backed = None
        self.set_dataset_groups(dataset_groups=dataset_groups)

    def get_gc(
            self,
            genome: str = None
    ):
        if genome.lower().startswith("homo_sapiens"):
            g = SuperGenomeContainer(
                species="human",
                genome=genome
            )
        elif genome.lower().startswith("mus_musculus"):
            g = SuperGenomeContainer(
                species="mouse",
                genome=genome
            )
        else:
            raise ValueError("genomes %s not recognised. please provide valid genomes." % genome)
        return g


    @property
    def ncells(self):
        return sum([x.ncells for x in self.dataset_groups])

    @property
    def ncells_bydataset(self):
        """
        List of list of length of all data sets by data set group.
        :return:
        """
        return [x.ncells_bydataset for x in self.dataset_groups]

    @property
    def ncells_bydataset_flat(self):
        """
        Flattened list of length of all data sets.
        :return:
        """
        return [xx for x in self.dataset_groups for xx in x.ncells_bydataset]

    def set_dataset_groups(self, dataset_groups: List[DatasetGroupBase]):
        self.dataset_groups = dataset_groups

    def subset_organs(self, subset: Union[None, List]):
        for x in self.dataset_groups:
            if x.datasets[0].organ == "mixed":
                x.subset_organs(subset)

    def load_all(
            self,
            celltype_version: Union[str, None] = None,
            match_to_reference: Union[str, None] = None,
            remove_gene_version: bool = True,
            annotated_only: bool = False,
            load_raw: bool = False
    ):
        """
        Loads data set groups into anndata object.

        :param celltype_version: Version of cell type ontology to use.
            Uses most recent within each DatasetGroup if None.
        """
        for x in self.dataset_groups:
            x.load_all(
                annotated_only=annotated_only,
                remove_gene_version=remove_gene_version,
                match_to_reference=match_to_reference,
                celltype_version=celltype_version,
                load_raw=load_raw
            )
        # making sure that concatenate is not used on a None adata object resulting from organ filtering
        for i in range(len(self.dataset_groups)):
            if self.dataset_groups[i].adata is not None:
                break
        self.adata = self.dataset_groups[i].adata.concatenate(
            *[x.adata for x in self.dataset_groups[1:] if x is not None],
            join="outer",
            batch_key=ADATA_IDS.dataset_group
        )

    def load_all_tobacked(
            self,
            fn_backed: PathLike,
            genome: str,
            shuffled: bool = False,
            as_dense: bool = False,
            annotated_only: bool = False,
            celltype_version: Union[str, None] = None,
    ):
        """
        Loads data set groups into backed anndata object.

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
        :param annotated_only:
        :param celltype_version: Version of cell type ontology to use. Uses most recent if None.
        """
        if shuffled and not as_dense:
            raise ValueError("cannot write backed shuffled and sparse")
        scatter_update = shuffled
        self.fn_backed = fn_backed
        ncells = self.ncells
        gc = self.get_gc(genome=genome)
        ngenes = gc.ngenes
        if scatter_update:
            self.adata = anndata.AnnData(
                scipy.sparse.csr_matrix((ncells, ngenes), dtype=np.float32)
            )  # creates an empty anndata object with correct dimensions that can be filled with cells from data sets
        else:
            self.adata = anndata.AnnData(
                scipy.sparse.csr_matrix((0, ngenes), dtype=np.float32)
            )
        self.adata.filename = fn_backed  # setting this attribute switches this anndata to a backed object
        # Note that setting .filename automatically redefines .X as dense, so we have to redefine it as sparse:
        if not as_dense:
            self.adata.X = scipy.sparse.csr_matrix(self.adata.X)  # redefines this backed anndata as sparse
        keys = [
            ADATA_IDS.lab,
            ADATA_IDS.year,
            ADATA_IDS.protocol,
            ADATA_IDS.organ,
            ADATA_IDS.subtissue,
            ADATA_IDS.cell_ontology_class,
            ADATA_IDS.state_exact,
            ADATA_IDS.normalization,
            ADATA_IDS.dev_stage,
            ADATA_IDS.has_celltypes,
            ADATA_IDS.dataset
        ]
        if scatter_update:
            self.adata.obs = pandas.DataFrame({
                k: ["nan" for x in range(ncells)] for k in keys
            })
        else:
            for k in keys:
                self.adata.obs[k] = []
        # Define index vectors to write to:
        idx_vector = np.arange(0, ncells)
        if shuffled:
            np.random.shuffle(idx_vector)
        idx_ls = []
        row = 0
        for x in self.ncells_bydataset:
            temp_ls = []
            for y in x:
                temp_ls.append(idx_vector[row:(row+y)])
                row += y
            idx_ls.append(temp_ls)
        print("checking expected and received data set sizes, rerun meta data generation if mismatch is found:")
        print(self.ncells_bydataset)
        print([[len(x) for x in xx] for xx in idx_ls])
        for i, x in enumerate(self.dataset_groups):
            x.load_all_tobacked(
                adata_backed=self.adata,
                genome=genome,
                idx=idx_ls[i],
                keys=keys,
                annotated_only=annotated_only,
                celltype_version=celltype_version
            )
        # Save obs separately as this is not included in backed h5ad.
        fn_backed_obs = ".".join(self.fn_backed.split(".")[:-1]) + "_obs.csv"
        self.adata.obs.to_csv(fn_backed_obs)

    def delete_backed(self):
        del self.adata
        self.adata = None
        os.remove(str(self.fn_backed))

    def load_cached_backed(self, fn: PathLike):
        self.adata = anndata.read(fn, backed='r')
