import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.sparse
import os
import re


def load(data_dir, sample_fn, **kwargs):
    if sample_fn == "microwell":
        fn = os.path.join(data_dir, "GSE93811_RAW.tar")
        dfs = []
        obs = []
        with tarfile.open(fn) as tar:
            for member in tar.getmembers():
                d = pd.read_csv(tar.extractfile(member), compression="gzip", header=0, index_col=None)
                for i in range(d.shape[0]):
                    obs.append(member.name.split("_")[1])
                dfs.append(d)
        adata = ad.AnnData(X=scipy.sparse.csr_matrix(np.vstack((dfs[0], dfs[1])), dtype=np.float32),
                           obs=pd.DataFrame({"location": obs}),
                           var=pd.DataFrame(index=dfs[0].columns.values))
        adata.obs["state_exact"] = "Non-fused organoids"
        adata.obs['organoid_age_days'] = "105"
        adata.obs['cell_line'] = "2242-1"
        adata.obs["protocol"] = adata.obs["location"].replace({
            "hCS": "Pasca, 2015 (doi: 10.1038/nmeth.3415)",
            "hSS": "Birey, 2017 (doi: 10.1038/nature22330)",
        })
    elif sample_fn == "Smart-seq":
        line_2_idxs = ['hSS507', 'hSS515', 'hSS536', 'hSS544', 'hSS627', 'hSS671', 'hSS683',
                       'hSS692', 'hSS693', 'hCS717', 'hCS730', 'hCS731', 'hCS732', 'hCS746',
                       'hCS752', 'hCS816', 'hCS823', 'hCS872', 'hCS889']
        fn_f2 = os.path.join(data_dir, "GSE93321_exps.txt.gz")  # fusion data, 2 weeks
        fn_f4 = os.path.join(data_dir, "GSE95837_exps.txt.gz")  # fusion data, 4 weeks
        adata_f2 = ad.AnnData(pd.read_csv(fn_f2, delimiter='\t', compression='gzip', index_col=0).T)
        adata_f4 = ad.AnnData(pd.read_csv(fn_f4, delimiter='\t', compression='gzip', index_col=0).T)
        adata_f2.obs['state_exact'] = 'Neurons from hSS and hCS oganoids fused for 2 weeks'
        adata_f4.obs['state_exact'] = 'Neurons from hSS and hCS oganoids fused for 4 weeks'
        adata_f2.obs['cell_line'] = '2242-1'
        adata_f4.obs['cell_line'] = ['2242-1' if i in line_2_idxs else '8858-1' for i in adata_f4.obs_names]
        adata_f2.obs['organoid_age_days'] = "104"
        adata_f4.obs['organoid_age_days'] = "118"
        adata = adata_f2.concatenate(adata_f4, index_unique=None)
        adata.X = scipy.sparse.csr_matrix(adata.X, dtype=np.float32).copy()
        del adata.obs["batch"]
        adata.obs['location'] = ["".join(re.findall("[a-zA-Z]+", i)) for i in adata.obs_names]
        adata.obs["protocol"] = adata.obs["location"].replace({
            "hCS": "Pasca, 2015 (doi: 10.1038/nmeth.3415)",
            "hSS": "Birey, 2017 (doi: 10.1038/nature22330)",
        })
    else:
        raise ValueError()
    return adata
