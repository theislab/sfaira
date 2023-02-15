import anndata
import os
import tarfile
import pandas as pd
import numpy as np



def load(data_dir, **kwargs):
    dataset_tar_path = os.path.join(data_dir, "GSE137941_RAW.tar")
    with tarfile.open(dataset_tar_path, "r") as tar:
        tar.extractall(path=data_dir)
    expression_path = os.path.join(data_dir, "GSM4094681_C1_150d_fseq12BC77_S2.deg.txt.gz")
    expression_df = pd.read_csv(expression_path,
                                compression='gzip', sep="\t", index_col=0, low_memory=False, header=0)
    adata = anndata.AnnData(X=expression_df.T,
                            obs=pd.DataFrame(expression_df.columns).set_index(0),
                            var=pd.DataFrame(expression_df.index).set_index("GENE"), dtype=np.int64)
    adata.obs.index.name = "cell"

    return adata
