from sfaira.data import DatasetBaseGroupLoadingManyFiles


class Dataset(DatasetBaseGroupLoadingManyFiles):

    def _load(self):
        import os
        import pandas as pd
        import anndata as ad
        import scipy.sparse
        import numpy as np

        fn = [
            os.path.join(self.data_dir, f"GSE114374_Human_{self.sample_fn}_expression_matrix.txt.gz"),
            os.path.join(self.data_dir, f"{self.sample_fn.lower()}_meta_data_stromal_with_donor.txt"),
        ]
        matrix = pd.read_csv(fn[0], sep="\t")
        obs = pd.read_csv(fn[1], sep="\t", index_col=3)
        adata = ad.AnnData(matrix.T)
        adata.X = scipy.sparse.csc_matrix(np.expm1(adata.X))
        adata.obs = obs
        adata.obs['state_exact'] = "healthy colon" if self.sample_fn == "HC" else "ulcerative colitis"
        s_dict = {"F": "female", "M": "male"}
        adata.obs['Sex'] = [s_dict[i] for i in adata.obs['Sex']]

        return adata
