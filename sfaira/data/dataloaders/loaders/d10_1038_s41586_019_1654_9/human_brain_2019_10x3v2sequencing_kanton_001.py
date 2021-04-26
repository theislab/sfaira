import anndata
import os
import scipy.io
import zipfile
import pandas


def load(data_dir, **kwargs):

    cell_line_dict = {
        '409b2': '409B2',
        'H9': 'WA09',
        'Wibj2': 'HPSI0214i-wibj_2',
        'Sc102a1': 'SC102A-1',
        'Kucg2': 'HPSI0214i-kucg_2',
        'Hoik1': 'HPSI0314i-hoik_1',
        'Sojd3': 'HPSI0314i-sojd_3',
    }

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
    adata.obs["Line"] = [cell_line_dict[x] for x in adata.obs["Line"]]

    return adata
