import anndata as ad
import os
import numpy as np
import pandas as pd
import scipy
import tarfile


# the data_dir argument will be automatically set by sfaira to the folder where your datafiles lie
def load(data_dir, **kwargs):

    fn = os.path.join(data_dir, "GSE168323_RAW.tar")

    adatas_list = []
    protocol_dict = {'standardorg': 'Fiorenzano, 2021 standard VM (doi: 10.1038/s41467-021-27464-5)',
                     'silk': 'Fiorenzano, 2021 VM with spider-silk microfibers (doi: 10.1038/s41467-021-27464-5)',
                     'silkLam': 'Fiorenzano, 2021 VM with spider-silk microfibers & full-length human laminin '
                     '(doi:10.1038/s41467-021-27464-5)'}

    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            d = pd.read_csv(tar.extractfile(member), compression="gzip", delimiter=',', header=0, index_col=0).T

            if "Silk" in member.name:
                age = member.name.split("Day")[1].split("-")[0]
                tech_sample = member.name.split('-')[1].split('.csv')[0]
                protocol = protocol_dict[member.name.split('-')[1][2:].split(".csv")[0]]
            elif "standardorg" in member.name:
                age = member.name.split("day")[1].split(".csv")[0]
                tech_sample = member.name.split("_")[1].split(".csv")[0]
                protocol = protocol_dict['standardorg']
            else:
                print("Error: organoid protocol not found in protocol dict")

            temp = ad.AnnData(X=scipy.sparse.csr_matrix(d, dtype=np.float32),
                              var=pd.DataFrame(index=d.columns.values),
                              obs=pd.DataFrame(index=d.index))
            temp.obs['organoid_age_days'] = age
            temp.obs['tech_sample'] = tech_sample
            temp.obs['protocol'] = protocol
            adatas_list.append(temp)

    adata = ad.concat(adatas_list, join='outer')

    return adata
