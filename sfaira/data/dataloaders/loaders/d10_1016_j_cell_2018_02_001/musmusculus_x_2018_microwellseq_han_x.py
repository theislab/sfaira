import anndata
import numpy as np
import pandas
import zipfile
import tarfile
import os


def load(data_dir, sample_fn, **kwargs):
    fn = os.path.join(data_dir, '5435866.zip')
    with zipfile.ZipFile(fn) as archive:
        celltypes = pandas.read_csv(archive.open('MCA_CellAssignments.csv'), index_col=1)
        celltypes = celltypes.drop(["Unnamed: 0"], axis=1)

        with tarfile.open(fileobj=archive.open('MCA_500more_dge.tar.gz')) as tar:
            data = pandas.read_csv(tar.extractfile(f'500more_dge/{sample_fn}'),
                                   compression="gzip",
                                   sep=" ",
                                   header=0
                                   )

    adata = anndata.AnnData(data.T)
    annotated_cells = np.array([x in celltypes.index for x in adata.obs_names])
    # Add annotation if available for this data set:
    if np.sum(annotated_cells) > 0:
        # Subset to annotated cells if any are annotated:
        adata = adata[annotated_cells].copy()
        # Clean nans in data frame to avoid issues with cell type annotation:
        celltypes["Annotation"] = [
            x if x not in [np.nan, "nan"] else "unknown"
            for x in celltypes["Annotation"].values
        ]
        adata.obs = celltypes.loc[adata.obs_names, :]

    return adata
