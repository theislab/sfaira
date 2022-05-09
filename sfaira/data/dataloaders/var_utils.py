from typing import Dict, List, Union

import anndata
import numpy as np
import pandas as pd
import scipy.sparse

from sfaira.consts import AdataIds
from sfaira.versions.genomes import GenomeContainer


def streamline_var(adata: anndata.AnnData,
                   adata_target_ids: AdataIds,
                   organism: str,
                   dataset_id: str = "",
                   clean_var: bool = False,
                   feature_id_var_key: Union[str, None] = None,
                   feature_symbol_var_key: Union[str, None] = None,
                   layer_counts: Union[str, None] = None,
                   layer_processed: Union[str, None] = None,
                   match_to_release: Union[str, Dict[str, str], None] = None,
                   matrix_format: str = "csr",
                   remove_gene_version: bool = True,
                   subset_genes_to_type: Union[None, str, List[str]] = None,
                   verbose: int = 1):
    """
    Subset and sort genes to genes defined in an assembly or genes of a particular type, such as protein coding.
    This also adds missing ENSEMBL ID or gene symbol columns if match_to_reference is not set to False and removes all
    adata.var columns that are not defined as feature_id_var_key or feature_symbol_var_key in the data loader.


    :param adata: anndata.AnnData instance to transform.
    :param adata_target_ids: AdataIds instance for output (transformed) adata. See also sfaira.consts.adata_fields for
        examples.
    :param clean_var: Whether to remove non-controlled fields in .var.
    :param dataset_id: Identifier of dataset used in logging.
    :param feature_id_var_key: Column name of ENSEMBL IDs in adata.var, use "index" if they are in the index and None if
        they are not included in adata. Supply either this feature_id_var_key feature_symbol_var_key.
    :param feature_symbol_var_key: Column name of gene symbols in adata.var, use "index" if they are in the index and
        None if they are not included in adata. Supply either this feature_id_var_key feature_symbol_var_key.
    :param layer_counts: Key of layer with count data, use "X" if these are in adata.X and "raw" if they are in
        adata.raw.X. Need to supply at least one of layer_counts and layer_processed.
    :param layer_processed: Key of layer with processed data, use "X" if these are in adata.X and "raw" if they are in
        adata.raw.X. Need to supply at least one of layer_counts and layer_processed.
    :param match_to_release: Which ENSEMBL genome annotation release to map the feature space to.
        Uses default set in schema if not set. Can be:
            - str: Provide the name of the release (eg. "104").
            - dict: Mapping of organism to name of the release (see str format). Chooses release for each
                data set based on organism annotation.
    :param matrix_format: Format to convert data matrices to.
    :param organism: Main organism in data that defines the genome used, use NCBI taxon elements.
    :param remove_gene_version: Whether to remove the version number after the colon sometimes found in ensembl
        gene ids.
        Uses default set in schema if not set.
    :param subset_genes_to_type: Type(s) to subset to. Can be a single type or a list of types or None.
        Uses default set in schema if not set.
        Types can be:
            - None: Keep the subset of input gene set that can be mapped to assembly.
            - "all": All genes in assembly.
            - "protein_coding": All protein coding genes in assembly.
    :param verbose: Report transformation statistics.
    :return:
    """
    if match_to_release is None:
        match_to_release = adata_target_ids.feature_kwargs["match_to_release"]
    if remove_gene_version is None:
        remove_gene_version = adata_target_ids.feature_kwargs["remove_gene_version"]
    if subset_genes_to_type is None:
        subset_genes_to_type = adata_target_ids.feature_kwargs["subset_genes_to_type"]

    # Define target feature space.
    # Set genome container if mapping of gene labels is requested
    if isinstance(match_to_release, dict):
        match_to_release = match_to_release[organism]
    genome_container = GenomeContainer(organism=organism, release=match_to_release)
    # Define target gene sets based on defined genome.
    if subset_genes_to_type is None:
        target_ids = None
    elif isinstance(subset_genes_to_type, str) and subset_genes_to_type == "all":
        target_ids = genome_container.ensembl
    else:
        if isinstance(subset_genes_to_type, str):
            subset_genes_to_type = [subset_genes_to_type]
        keys = np.unique(genome_container.biotype)
        if subset_genes_to_type not in keys:
            raise ValueError(f"subset type {subset_genes_to_type} not available in list {keys}")
        target_ids = np.asarray([
            x.upper() for x, y in zip(genome_container.ensembl, genome_container.biotype)
            if y in subset_genes_to_type
        ])
    allowed_ids = np.asarray([x.upper() for x in genome_container.ensembl]) if target_ids is None else None

    # Extract layers from loaded object in order to manipulate feature dimension of data matrices.
    # Sfaira controls a processed and a counts layer, both are extracted here.
    # Any change to the size of the feature dimension of the data matrix also requires changes to .var which is also
    # extracted here.
    assert layer_processed or layer_counts
    if layer_processed is None:
        x_proc = None
        var_proc = None
    elif layer_processed == "X":
        x_proc = adata.X
        var_proc = adata.var
    elif layer_processed == "raw":
        x_proc = adata.raw.X
        var_proc = adata.raw.var
    else:
        if layer_processed not in adata.layers.keys():
            raise ValueError(f"layer {layer_processed} not in layers {list(adata.layers.keys())}")
        else:
            x_proc = adata.layers[layer_processed]
            var_proc = adata.var
    if layer_counts is None:
        x_counts = None
        var_counts = None
    elif layer_counts == "X":
        x_counts = adata.X
        var_counts = adata.var
    elif layer_counts == "raw":
        x_counts = adata.raw.X
        var_counts = adata.raw.var
    else:
        if layer_counts not in adata.layers.keys():
            raise ValueError(f"layer {layer_counts} not in layers {list(adata.layers.keys())}")
        else:
            x_counts = adata.layers[layer_counts]
            var_counts = adata.var

    # Create new adata using extracted layers.
    # First, we prepare moving data matrices to new, controlled location in .adata.
    # Note: one layer needs to be in X.
    if x_proc is not None and x_counts is not None:
        layer_counts = adata_target_ids.layer_counts
        layer_proc = adata_target_ids.layer_proc
        assert layer_counts == "X" or layer_proc == "X"
        if layer_counts == "X" and layer_proc != "X":
            x1 = x_counts
            var1 = var_counts
            x2 = x_proc
            var2 = var_proc
            layer2 = layer_proc
        elif layer_counts != "X" and layer_proc == "X":
            x1 = x_proc
            var1 = var_proc
            x2 = x_counts
            var2 = var_counts
            layer2 = layer_counts
        else:
            raise ValueError(f"one layer needs to be in X: layer_counts={layer_counts}, layer_proc={layer_proc}")
    elif x_proc is None and x_counts is not None:
        x1 = x_counts
        var1 = var_counts
        x2 = None
        var2 = None
        layer2 = None
        layer_counts = "X"
        layer_proc = None
    elif x_proc is not None and x_counts is None:
        x1 = x_proc
        var1 = var_proc
        x2 = None
        var2 = None
        layer2 = None
        layer_counts = None
        layer_proc = "X"
    else:
        raise ValueError("Neither layer_counts nor layer_proc are set in yaml. Aborting")
    # Second, we streamline the feature dimension of these data matrices to make sure all necessary feature
    # identifiers are provided:
    var1, feature_id_var_key, feature_symbol_var_key = format_var(
        adata_ids=adata_target_ids, clean_var=clean_var, feature_id_var_key=feature_id_var_key,
        feature_symbol_var_key=feature_symbol_var_key, gc=genome_container, var=var1)
    # Only need to process var2 if x2 has a separate feature space, ie if x2 is in .raw:
    if layer2 is not None and layer2 == "raw":
        # Note: Assume that IDs and var keys in .adata.raw.var are the same as in .adata.var.
        # This is tested here:
        if feature_id_var_key is not None:
            ids_counts = var2.index.values if feature_id_var_key == "index" \
                else var2[feature_id_var_key].values
            ids_proc = var2.index.values if feature_id_var_key == "index" \
                else var2[feature_id_var_key].values
        elif feature_symbol_var_key is not None:
            ids_counts = var2.index.values if feature_symbol_var_key == "index" \
                else var2[feature_symbol_var_key].values
            ids_proc = var2.index.values if feature_symbol_var_key == "index" \
                else var2[feature_symbol_var_key].values
        else:
            raise ValueError("Neither feature_id_var_key nor feature_symbol_var_key are set in yaml. Aborting")
        if set(ids_proc) - set(ids_counts):
            raise IndexError(f"Features of layer specified as `layer_processed` ('{layer_processed}') are "
                             f"not a subset of the features of layer specified as `layer_counts` "
                             f"('{layer_counts}'). This is not supported.")
        var2, _, _ = format_var(adata_ids=adata_target_ids, clean_var=clean_var, feature_id_var_key=feature_id_var_key,
                                feature_symbol_var_key=feature_symbol_var_key, gc=genome_container, var=var2)
    # Next, we impute and subset the feature dimension based on the target gene sets:
    x1_sum = np.log(x1.sum() + 1.) / np.log(10)
    x1_nonzero = (x1.sum(axis=0) > 0).sum()
    x1 = convert_matrix_format(x=x1, matrix_format=matrix_format)
    x1, var1 = reorder_adata_to_target_features(
        allowed_ids=allowed_ids,
        map_in_upper_case=True,
        map_without_version=True,
        target_ids=target_ids,
        var=var1,
        var_key=feature_id_var_key,
        x=x1)
    # Need to clean again to populate all gene identifier fields other than feature_id_var_key which were
    # added to the object, ie imputed features. There are NaN after reorder_adata_to_target_features in all fields
    # other than feature_id_var_key.
    var1, _, _ = format_var(adata_ids=adata_target_ids, feature_id_var_key="index", feature_symbol_var_key=None,
                            gc=genome_container, var=var1)
    x1_new_sum = np.log(x1.sum() + 1.) / np.log(10)
    x1_new_nonzero = (x1.sum(axis=0) > 0).sum()
    if layer2 is not None:
        x2_sum = np.log(x2.sum() + 1.) / np.log(10)
        x2_nonzero = (x2.sum(axis=0) > 0).sum()
        x2 = convert_matrix_format(x=x2, matrix_format=matrix_format)
        x2, var2 = reorder_adata_to_target_features(
            allowed_ids=allowed_ids,
            map_in_upper_case=True,
            map_without_version=True,
            target_ids=target_ids,
            var=var2,
            var_key=feature_id_var_key,
            x=x2)
        # Need to clean again to populate all gene identifier fields other than feature_id_var_key which were
        # added to the object, ie imputed features. There are NaN after reorder_adata_to_target_features in all
        # fields other than feature_id_var_key.
        var2, _, _ = format_var(adata_ids=adata_target_ids, feature_id_var_key="index", feature_symbol_var_key=None,
                                gc=genome_container, var=var2)
        x2_new_sum = np.log(x2.sum() + 1.) / np.log(10)
        x2_new_nonzero = (x2.sum(axis=0) > 0).sum()
    # Last, we build a new .adata instance from these manipulated data matrices.
    adata = anndata.AnnData(
        X=x1,
        obs=adata.obs,
        obsm=adata.obsm,
        obsp=adata.obsp,
        var=var1,
        uns=adata.uns,
        dtype="float32"
    )
    if layer2 is not None:
        if layer2 != "raw":
            adata.layers[layer2] = x2
        else:
            adata.raw = anndata.AnnData(x2, obs=pd.DataFrame({}, index=adata.obs_names), var=var2, dtype="float32")
    layer_counts = layer_counts
    layer_processed = layer_proc
    if hasattr(adata_target_ids, "mapped_features") and adata_target_ids.mapped_features is not None:
        adata.uns[adata_target_ids.mapped_features] = match_to_release
    # Reporting:
    if verbose > 0:
        report = f"transformed feature space on {dataset_id}: \n" \
                 f"log10 total counts {round(x1_sum, 2)} to {round(x1_new_sum, 2)}, " \
                 f"non-zero features {x1_nonzero} to {x1_new_nonzero}"
        if layer2 is not None:
            report += f"\n(secondary layer: log10 total total counts {round(x2_sum, 2)} to {round(x2_new_sum, 2)}, " \
                      f"non-zero features {x2_nonzero} to {x2_new_nonzero})"
        print(report)
    return {"adata": adata,
            "genome_container": genome_container,
            "layer_counts": layer_counts,
            "layer_processed": layer_processed}


def format_var(adata_ids: AdataIds,
               feature_id_var_key: Union[None, str],
               feature_symbol_var_key: Union[None, str],
               gc: GenomeContainer,
               var: pd.DataFrame,
               clean_var: bool = False):
    """
    Maps feature names in var table to feature names in ontology, imputing symbols or IDs respectively if necessary.

    Note on gene versions and case-sensitivity:
    In mapping of IDs to symbols and reverse, this is handled by methods of the GenomeContainer.
    In the output, both are defined by the assembly set on the GenomeContainer instance.

    :param adata_ids: Field container.
    :param clean_var: Whether to remove non-controlled fields in .var.
    :param feature_symbol_var_key: Location of feature symbol column in input var (column name or "index").
    :param feature_id_var_key:  Location of feature ID column in input var (column name or "index").
    :param gc: Genome container for feature identifier translation.
    :param var: Feature table to modify.
    :return: Tuple of
            - Modified table.
            - New location of IDs (column name or "index).
            - New location of symbols (column name or "index).
    """
    if feature_symbol_var_key is None and feature_id_var_key is None:
        raise ValueError("Either feature_symbol_var_key or feature_id_var_key needs to be provided in the "
                         "data loader")
    elif feature_id_var_key:
        ensids = var.index if feature_id_var_key == "index" else var[feature_id_var_key].values
        # Add new feature identifier:
        var[adata_ids.feature_symbol] = gc.translate_id_to_symbols(x=ensids)
        # Place existing identifier under streamlined key:
        if feature_id_var_key != adata_ids.feature_id:
            if adata_ids.feature_id == "index":
                var.index = ensids
            else:
                var[adata_ids.feature_id] = ensids
            # Delete old entry:
            if feature_id_var_key != "index":
                del var[feature_id_var_key]
    elif feature_symbol_var_key and feature_id_var_key is None:
        symbols = var.index if feature_symbol_var_key == "index" else var[feature_symbol_var_key]
        # Add new feature identifier:
        var[adata_ids.feature_id] = gc.translate_symbols_to_id(x=symbols)
        # Place existing identifier under streamlined key:
        if feature_symbol_var_key != adata_ids.feature_symbol:
            if adata_ids.feature_symbol == "index":
                var.index = symbols
            else:
                var[adata_ids.feature_symbol] = symbols
            # Delete old entry:
            if feature_symbol_var_key != "index":
                del var[feature_symbol_var_key]
    # Optionally remove all non-controlled fields:
    if clean_var:
        allowed_columns = [getattr(adata_ids, x) for x in adata_ids.var_keys]
        for k in list(var.columns):
            if k not in allowed_columns:
                del var[k]
    return var, adata_ids.feature_id, adata_ids.feature_symbol


def reorder_adata_to_target_features(
        map_in_upper_case: bool,
        map_without_version: bool,
        var: pd.DataFrame,
        var_key: str,
        x,
        allowed_ids: Union[np.ndarray, None] = None,
        target_ids: Union[np.ndarray, None] = None,
):
    """
    Reorders gene-dimensional x and var matrices to a target set of genes.

    If allowed_ids is used, the output matrix does not have a defined feature dimension but rather the length of the
    intersection between allowed_ids and the found IDs in var[var_key].
    if target_ids is used, the non-zero features in the output are the intersection between allowed_ids and the found
    IDs in var[var_key], but additionally all remaining target_ids are imputed as zero and the feature dimension is
    reordered to match target_ids.

    :param allowed_ids: Allowed identifiers, can use this to subset if target identifiers is not set.
        Note: These have to match the identifier class used in var[var_key], eg. gene symbols or ENSG IDs.
        Note: Original feature ordering is kept if allowed_ids is not None and target_ids is None.
    :param map_in_upper_case: Whether to map features in upper case, ie case-insensitive.
    :param map_without_version: Whether to map features in without accounting version (".* suffix of feature ID).
    :param target_ids: Target set of feature identifiers.
        Note: These have to match the identifier class used in var[var_key], eg. gene symbols or ENSG IDs.
    :param var: Feature meta data table (.var in anndata).
    :param var_key: Key of feature identifiers in var table.
    :param x: Expression matrix (.X or layer in anndata).
    :return: Tuple of processed x and var.
    """
    # Process gene annotations
    # Make features unique (to avoid na-matches in converted columns to be collapsed below.
    x, var = collapse_x_var_by_feature(sep_deduplication="-", var=var, var_column=var_key, x=x)
    input_ids = var.index if var_key == "index" else var[var_key].values
    # Set index of var to var_key if not already the case, this is useful to avoid index duplication events in
    # rearrangements.
    if var_key != "index":
        var.index = input_ids
    if allowed_ids is not None and target_ids is None:
        # Keeping input ordering:
        output_ids = np.array([z for z in input_ids if z in allowed_ids])
    elif allowed_ids is None and target_ids is not None:
        output_ids = target_ids
    else:
        raise ValueError(f"supply either allowed_ids or target_ids: {(allowed_ids, target_ids)}")
    map_mat = scipy.sparse.lil_matrix(np.zeros((len(input_ids), len(output_ids))))
    if map_in_upper_case:
        input_ids = np.array([z.upper() if isinstance(z, str) else "" for z in input_ids])
        output_ids = np.array([z.upper() if isinstance(z, str) else "" for z in output_ids])
    if map_without_version:
        input_ids = np.array([z.split(".")[0] for z in input_ids])
        output_ids = np.array([z.split(".")[0] for z in output_ids])
    for i, z in enumerate(output_ids):
        # TODO: this is an equality map of feature identifiers that could be replaced by ENSEMBL informed maps.
        map_mat[np.where(input_ids == z)[0], i] = 1.
    map_mat = map_mat.tocsr()
    x, var = reorder_x_var(map_mat=map_mat, var=var, var_index_target=output_ids, x=x)
    return x, var


def collapse_x_var_by_feature(x, var, var_column, sep_deduplication="-"):
    """
    Collapses (sum) features with the same var_name in a provided var column.

    Keeps .var column of first occurrence of duplicated variables.
    Does not modify data if index is already unique.
    Reverses potential previous deduplication of feature names via declared dedepulication suffixes.
    These would have been introduced by methods such as:
        https://anndata.readthedocs.io/en/refpaths/anndata.AnnData.var_names_make_unique.html.

    :param sep_deduplication: Separator for feature ID suffixes that may have been added to make feature names unique.
    :param var: Input var object.
    :param var_column: column name in .var that contains the duplicated features of interest
    :param x: Input data matrix.
    :return: Processed x and var.
    """
    old_index = var.index if var_column == "index" else var[var_column].values
    old_index = np.array([z.split(sep_deduplication)[0] for z in old_index])
    # Get unique elements maintaining original ordering:
    idx = np.unique(old_index, return_index=True)[1]
    new_index = old_index[np.sort(idx)]
    if len(new_index) < len(old_index):
        map_mat = scipy.sparse.lil_matrix(np.zeros((len(old_index), len(new_index))))
        for i, z in enumerate(new_index):
            map_mat[np.where(old_index == z)[0], i] = 1.
        map_mat = map_mat.tocsr()
        x, var = reorder_x_var(map_mat=map_mat, var=var, var_index_target=new_index, x=x)
    return x, var


def reorder_x_var(
        map_mat: Union[scipy.sparse.spmatrix, None],
        var: pd.DataFrame,
        var_index_target: np.ndarray,
        x: scipy.sparse.spmatrix,
):
    """
    Reorders gene-dimensional x and var matrices based on a mapping matrix.

    :param map_mat: Spare mapping matrix: (number of input features, number of target features).
        Setup:
            - Every entry is in [0, 1] and signifies the proportion of counts in the corresponding input gene (row) that
                is carried over to the designated target gene (column).
            - The sum over every row and column is 0 or 1 for 1:1, 0:1, and 1:0 input:output feature maps.
            - For 1:many and many:1 input:output feature maps, the sum over rows or columns may deviate.
            - Rows with all zero entries correspond to features that do not appear in the output feature space.
            - Columns with all zero entries correspond to features that do not appear in the input feature space.
            - map_mat is all-zero off-diagonal and has ones on the diagonal iff input == output features.
        Note on var manipulation:
            - 0:1 results in creation of a new entry in var with NaN entries.
            - many:1 keep var entry of first occurrence only. ToDo: This may not be desired in edge cases.
            - 1:1 and 1:many results in preservation of the corresponding input feature column in var.
            - *:0 results in deletion of the corresponding var columns.
            - ToDo many:many are not handled separately yet.
    :param var: Feature meta data table (.var in anndata).
    :param var_index_target: Index out output var, length of number of target features.
    :param x: Expression matrix (.X or layer in anndata).
    :return: Tuple of processed x and var.
    """
    assert len(var_index_target) == map_mat.shape[1]
    # Move data matrix to sparse:
    if isinstance(x, np.ndarray):
        # Change NaN to zero. This occurs for example in concatenation of anndata instances.
        if np.any(np.isnan(x)):
            x[np.isnan(x)] = 0
        x = scipy.sparse.csr_matrix(x)
    elif isinstance(x, scipy.sparse.spmatrix):
        pass
    else:
        raise ValueError(f"Data type {type(x)} not recognized.")

    row_idx, col_idx, _ = scipy.sparse.find(map_mat)
    is_same_size = map_mat.shape[0] == map_mat.shape[1]
    is_same_ordering = np.all(row_idx == np.arange(0, len(row_idx))) and np.all(col_idx == np.arange(0, len(col_idx)))
    # Only modify if map is not diagonal:
    if not is_same_size or not is_same_ordering:
        # Apply the transport map to x:
        x_new = x.dot(map_mat)
        # Modify var:
        idx = []
        for i in range(map_mat.shape[1]):
            matches = np.where(col_idx == i)[0]
            if len(matches) == 1:
                # 1:1 maps
                idx.append(row_idx[matches[0]])
            elif len(matches) > 1:
                # For many:1 maps, only the first occuring var row in the input is used
                # (note, this is a collapse of multiple rows to one row.)
                # ToDo could check if this collapse destroys variability in these rows.
                idx.append(row_idx[matches[0]])
            # 0:1 are imputed as nan.
        idx = np.array(idx)
        var_new = var.iloc[idx, :]
        var_new = var_new.reindex(var_index_target)
    else:
        x_new = x
        var_new = var

    return x_new, var_new


def convert_matrix_format(x, matrix_format: str):
    sp_cls = {"coo": scipy.sparse.coo_matrix,
              "csc": scipy.sparse.csc_matrix,
              "csr": scipy.sparse.csr_matrix}
    if isinstance(x, np.ndarray):
        if matrix_format in sp_cls.keys():
            x = sp_cls[matrix_format](x)
        elif matrix_format == "numpy":
            pass
        else:
            raise ValueError(f"could not convert {type(x)} to {matrix_format}")
    elif isinstance(x, scipy.sparse.spmatrix):
        if matrix_format in sp_cls.keys() or matrix_format == "numpy":
            x = sp_cls[matrix_format](x)
        else:
            raise ValueError(f"could not convert {type(x)} to {matrix_format}")
    return x
