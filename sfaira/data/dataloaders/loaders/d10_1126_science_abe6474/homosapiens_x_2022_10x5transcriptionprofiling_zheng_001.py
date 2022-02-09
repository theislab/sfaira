import anndata
import gzip
import os
import pandas as pd
import tarfile
import scanpy as sc
import scipy.sparse


def ir_read_from_chain_gz(fn, sep, barcodes, id_col, locus_col, v_call_col, d_call_col, j_call_col, c_call_col,
                          productive_col, junction_col, junction_aa_col, consensus_count_col, duplicate_count_col):
    """
    Util function to read the VDJ files with custom record format

    Replaces ir.io.read_airr.
    ToDo: Can this code be moved to scirpy? This file format does not seem very rare.

    E.g. of the form d10_1126_science_abe6474::GSE156728_10X_VDJ.merge.txt.gz.
    Issues with ir.io.read_airr on this file are:

        - gzip compression
        - non-standard column names
        - selected extraction of cells into AIRR format to save time
    """
    from scirpy.io._io import from_airr_cells, AirrCell, DEFAULT_AIRR_CELL_ATTRIBUTES, DEFAULT_AIRR_FIELDS

    tab = pd.read_csv(fn, sep=sep, compression="gzip")
    airr_cells = {}
    for bc in barcodes:
        airr_cells[bc] = AirrCell(
            cell_id=bc,
            cell_attribute_fields=DEFAULT_AIRR_CELL_ATTRIBUTES,
        )
        # Get all chains for this barcode:
        tab_i = tab.loc[tab[id_col].values == bc, :]
        for i in range(tab_i.shape[0]):
            chain = {
                "productive": tab_i[productive_col].values[i],
                "locus": tab_i[locus_col].values[i],
                "v_call": tab_i[v_call_col].values[i],
                "d_call": tab_i[d_call_col].values[i],
                "j_call": tab_i[j_call_col].values[i],
                "c_call": tab_i[c_call_col].values[i],
                "junction": tab_i[junction_col].values[i],
                "junction_aa": tab_i[junction_aa_col].values[i],
                "consensus_count": tab_i[consensus_count_col].values[i],
                "duplicate_count": tab_i[duplicate_count_col].values[i],
                "sequence_id": f"{bc}_{i}",
                "sequence": None,
                "rev_comp": None,
                "sequence_alignment": None,
                "germline_alignment": None,
                "v_cigar": None,
                "d_cigar": None,
                "j_cigar": None,
            }
            airr_cells[bc].add_chain(chain)

    return from_airr_cells(airr_cells.values(), include_fields=DEFAULT_AIRR_FIELDS)


def load(data_dir, sample_fn, **kwargs):
    import scirpy as ir

    fn = os.path.join(data_dir, sample_fn + ".counts.txt.gz")
    # Some of the count matrices are in a tar archive, the rest is given as individual downloads
    if sample_fn.startswith("GSE156728_RAW/"):
        fn_tar = os.path.join(data_dir, "GSE156728_RAW.tar")
        with tarfile.open(fn_tar) as tar:
            with gzip.open(tar.extractfile(sample_fn.split("GSE156728_RAW/")[-1] + ".counts.txt.gz"), "rb") as f:
                tab = pd.read_csv(f, sep="\t", index_col=0).T
            adata = anndata.AnnData(tab)
    else:
        adata = sc.read(fn).transpose()
        adata.X = scipy.sparse.csr_matrix(adata.X)
    fn_meta = os.path.join(data_dir, "GSE156728_metadata.txt.gz")
    tab_meta = pd.read_csv(fn_meta, sep="\t", compression="gzip")
    tab_meta.index = tab_meta["cellID"].values
    del tab_meta["cellID"]
    adata.obs = tab_meta.loc[adata.obs_names, :]
    fn_vdj = os.path.join(data_dir, "GSE156728_10X_VDJ.merge.txt.gz")
    adata_tcr = ir_read_from_chain_gz(
        fn=fn_vdj,
        barcodes=adata.obs_names,
        sep="\t",
        id_col="barcode",
        locus_col="chain",
        v_call_col="v_gene",
        d_call_col="d_gene",
        j_call_col="j_gene",
        c_call_col="c_gene",
        productive_col="productive",
        junction_col="cdr3_nt",
        junction_aa_col="cdr3",
        consensus_count_col="umis",
        duplicate_count_col="reads",
    )
    ir.pp.merge_with_ir(adata, adata_tcr)
    return adata
