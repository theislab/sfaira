import anndata
import os
import scipy.io
import zipfile
import pandas


def load(data_dir, **kwargs):
    fn = [
        os.path.join(data_dir, "E-MTAB-7552.processed.3.zip"),
        os.path.join(data_dir, "E-MTAB-7552.processed.1.zip"),
        os.path.join(data_dir, "E-MTAB-7552.processed.7.zip")
    ]
    with zipfile.ZipFile(fn[0]) as archive:
        x = scipy.io.mmread(archive.open('human_cell_counts_GRCh38.mtx')).T.tocsr()
    with zipfile.ZipFile(fn[1]) as archive:
        var = pandas.read_csv(archive.open('genes_GRCh38.txt'), sep="\t", index_col=1, names=['ensembl', 'genetype'])
    with zipfile.ZipFile(fn[2]) as archive:
        obs = pandas.read_csv(archive.open('metadata_human_cells.tsv'), sep="\t", index_col=0)
    adata = anndata.AnnData(X=x, var=var, obs=obs)

    # TODO: remove non-protein coding genes?

    return adata
