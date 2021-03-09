import anndata
import os
import scipy.sparse
import numpy as np


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "")

    fn = "../data/d10_1038_s41586_019_1654_9/E-MTAB-7552.processed.3.zip"
    with zipfile.ZipFile(fn) as archive:
        x = scipy.io.mmread(archive.open('human_cell_counts_GRCh38.mtx')).T.tocsr()

    fn = "../data/d10_1038_s41586_019_1654_9/E-MTAB-7552.processed.1.zip"
    with zipfile.ZipFile(fn) as archive:
        var = pd.read_csv(archive.open('genes_GRCh38.txt'), sep="\t", index_col=1, names=['ensembl', 'genetype'])

    fn = "../data/d10_1038_s41586_019_1654_9/E-MTAB-7552.processed.7.zip"
    with zipfile.ZipFile(fn) as archive:
        obs = pd.read_csv(archive.open('metadata_human_cells.tsv'), sep="\t", index_col=0)

    adata = ad.AnnData(X=x, var=var, obs=obs)

    # remove non-protein coding genes?

    return adata
