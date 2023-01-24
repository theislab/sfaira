import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.sparse
import os
import re
import scanpy as sc


def load(data_dir, sample_fn, **kwargs):
    if sample_fn == "GSE93811_RAW.tar":
        fn = os.path.join(data_dir, "GSE93811_RAW.tar")
        dfs = []
        obs = []
        with tarfile.open(fn) as tar:
            for member in tar.getmembers():
                d = pd.read_csv(tar.extractfile(member), compression="gzip", header=0, index_col=None)
                for i in range(d.shape[0]):
                    obs.append(member.name.split("_")[1])
                dfs.append(d)

        adata = ad.AnnData(X=scipy.sparse.csr_matrix(np.vstack((dfs[0], dfs[1]))),
                           obs=pd.DataFrame({"sample": obs}),
                           var=pd.DataFrame(index=dfs[0].columns.values))
        adata.obs['organoid_age_days'] = "105"
        adata.obs['cell_line'] = "2242-1"
    else:
        fn_f2 = os.path.join(data_dir, "GSE93321_exps.txt.gz")# fusion data, 2 weeks 37
        fn_f4 = os.path.join(data_dir, "GSE95837_exps.txt.gz")# fusion data, 4 weeks
        f2 = pd.read_csv(fn_f2,delimiter='\t', compression='gzip', index_col=0)
        f4 = pd.read_csv(fn_f4,delimiter='\t', compression='gzip', index_col=0)

        adata_f2 = ad.AnnData(f2.T)
        adata_f2.obs['migration_time'] = '2 weeks'
        adata_f2.obs['cell_line'] = '2242-1'
        adata_f4 = ad.AnnData(f4.T)
        adata_f4.obs['migration_time'] = '4 weeks'
        adata_f4.obs['cell_line'] = '2242-1/8858-1'
        adata = adata_f2.concatenate(adata_f4, index_unique=None)
        
        annot = sc.queries.biomart_annotations(
            "hsapiens",
            ["ensembl_gene_id", "external_gene_name"],
            host="grch37.ensembl.org",
        ).set_index("ensembl_gene_id")
        adata.var[annot.columns] = annot
        adata.obs['sample'] = 'hSS'
        adata.obs['organoid_age_days'] = "90"

        location = []
        for i in range(len(adata.obs)):
            temp = re.findall("[a-zA-Z]+", adata.obs_names[i])
            temp = "".join(temp)
            location.append(temp)
        adata.obs['location'] = location
    return adata