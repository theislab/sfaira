import gzip
import os
import pandas as pd
import scanpy as sc
import tarfile


def load(data_dir, sample_fn, **kwargs):
    import scirpy as ir

    fn = os.path.join(data_dir, "GSE171037_RAW.tar")
    with tarfile.open(fn) as tar:
        # TODO debug reading of h5 from tar
        adata_rna = sc.read_10x_h5(tar.extractfile(sample_fn + "_filtered_feature_bc_matrix.h5"))
        with gzip.open(tar.extractfile(sample_fn + "_all_contig_annotations.json.gz"), "rb") as f_json:
            adata_tcr = ir.io.read_10x_vdj(f_json, filtered=True)
        ir.pp.merge_with_ir(adata_rna, adata_tcr)
        tab = pd.read_csv(tar.extractfile(sample_fn + "_filtered_contig_annotations.csv.gz"), compression="gzip")
        sample_id = "_".join(sample_fn.split("_")[1:])
        tab["barcode_adjusted"] = ["-".join([x, sample_id]) for x in tab["barcode"].values]
        adata_rna.obs["clonotype_cellranger"] = [
            tab.loc[tab["barcode_adjusted"].values == xx, 'raw_clonotype_id'].values[0] + "_" + sample_id
            if xx in tab["barcode_adjusted"].values else "none"
            for xx in adata_rna.obs_names
        ]
    fn_meta = os.path.join(data_dir, "GSE171037_TcellReversePT_metadata_munich.txt.gz")
    # TODO match meta data with RNA object
    return adata_rna
