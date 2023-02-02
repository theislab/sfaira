import os
import numpy as np
import pandas as pd
import gzip
import scipy
import tarfile
import anndata as ad
from scipy.io import mmread


# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
def load(data_dir, **kwargs):

    fn = os.path.join(data_dir, "GSE149931_RAW.tar")
    barcodes_dict = {}
    features_dict = {}
    mtx_dict = {}
    adatas_list = []

    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            sample_key = (member.name.split("_")[1] + '_' +
                          member.name.split("_")[2] + '_' +
                          member.name.split("_")[3])

            if 'barcodes' in member.name:
                barcodes_df = pd.read_csv(tar.extractfile(member), compression="gzip", header=None).set_index(0)
                barcodes_df.index.name = None
                barcodes_dict[sample_key] = barcodes_df

            elif 'features' in member.name:
                features_df = pd.read_csv(tar.extractfile(member), compression="gzip", sep='\t', header=None)
                features_df = features_df[[0, 1]].set_index(1)
                features_df.index.name = None
                features_df = features_df.rename(columns={0: "EnsemblID"})
                features_dict[sample_key] = features_df

            elif 'matrix' in member.name:
                mtx = mmread(gzip.open(tar.extractfile(member), 'rb'))
                mtx_dict[sample_key] = mtx

            else:
                print("Error: unexpected file detected")

    for key in barcodes_dict:
        temp = ad.AnnData(X=scipy.sparse.csr_matrix(mtx_dict[key]).T,
                          var=features_dict[key],
                          obs=barcodes_dict[key])

        if key == 'hStrS_d80-83_aggr':
            temp.obs['cell_line'] = 'aggregated (2242-1, 1205-4, 8858-3)'
            temp.obs['organoid_age_days'] = '80-83'
            temp.obs['age'] = '80-83'
        else:
            temp.obs['cell_line'] = key.split('_')[0]
            temp.obs['organoid_age_days'] = key.split('_')[1][1:]
            temp.obs['age'] = key.split('_')[1][1:]

        temp.obs['sample'] = key
        temp.var_names_make_unique()
        temp.obs_names_make_unique()

        adatas_list.append(temp)
        del temp

    del barcodes_dict, features_dict, mtx_dict
    adata = ad.concat(adatas_list, join='outer')

    # add back EnsemblIDs because they are lost on concatenation
    EnsemblIDs = list(features_df['EnsemblID'])
    adata.var['EnsemblID'] = EnsemblIDs

    return adata
