import anndata
import os
import scipy.sparse
import tarfile
import pandas as pd
import scipy.sparse
import h5py
import yaml
import re

def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, "GSE118127_RAW.tar")
    fn_annot = os.path.join(data_dir, "cell_id_to_type.tsv")
    cell_annot = pd.read_csv(fn_annot, sep = "\t")
    cell_annot.set_index(keys = "cellID", inplace = True)
    print(cell_annot.head())
    with tarfile.open(fn) as tar:
        f = h5py.File(tar.extractfile(sample_fn), 'r')['GRCh38']
        x = scipy.sparse.csc_matrix((f['data'], f['indices'], f['indptr']), shape=f['shape']).T
        ensembl = [id.decode('UTF-8') for id in f['genes']]
        symbol = [id.decode('UTF-8') for id in f['gene_names']]
        var = pd.DataFrame({'feature_id': ensembl, 'feature_symbol': symbol})

        # Get sample IDs from the sample file name
        regex = "GSM[0-9]+_sample[_]?(.+)_filtered_gene_bc_matrices_h5.h5"
        result = re.search(regex, sample_fn)
        id = result.groups()[0]
        id_cor = id.replace("-", ".")
        id_cor = re.sub("_B.+", "", id_cor)
        if len(id_cor) == 1:
            id_cor = f"0{id_cor}"

        barcodes = []
        for brcd in f['barcodes']:
            barcode = brcd.decode('UTF-8').replace("-1", "")
            barcodes.append(id_cor + "_" + barcode)
        obs = pd.DataFrame({'barcode': barcodes})
        obs_cell_types = obs.join(cell_annot, on = 'barcode', how = 'left')
        obs_cell_types['organ'] = 'ovary'
        adata = anndata.AnnData(X=x, obs=obs_cell_types, var=var)
        
        # NaN introduced if a cell not in annotation. Remove these cells
        adata = adata[adata.obs["clusterID"].notna(), :]
    return adata
