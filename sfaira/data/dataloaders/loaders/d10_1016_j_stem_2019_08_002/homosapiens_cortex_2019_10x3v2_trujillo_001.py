import numpy as np
import pandas as pd
import tarfile
import anndata as ad
import scipy.sparse
import os
import scanpy as sc
import re
from tempfile import TemporaryDirectory
import shutil
import gzip


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE130238_RAW.tar")
    barcodes = []
    genes = []
    matrix = []
    sample = []

    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            temp = member.name
            if re.search('barcodes', temp):
                bar = pd.read_csv(tar.extractfile(member), delimiter='\t', compression="gzip", header=None, index_col=0)
                bar.index.name = None
                bar['time'] = member.name.split("_")[1]
                if member.name.split("_")[1] == "1M":
                    bar['organoid_age_days'] = "30"
                elif member.name.split("_")[1] == "3M":
                    bar['organoid_age_days'] = "90"
                elif member.name.split("_")[1] == "6M":
                    bar['organoid_age_days'] = "180"
                else:
                    bar['organoid_age_days'] = "300"
                barcodes.append(bar)
            elif re.search('genes', temp):
                gen = pd.read_csv(tar.extractfile(member), delimiter='\t', compression='gzip', index_col=0, header = None, names=["ensembl_gene_id", "gene_id"])
                genes.append(gen)
            else:
                fn_matrix = tar.extractfile(member)
                uncompressed_file_type = "mtx"
                with TemporaryDirectory() as tmpdir:
                    tmppth = tmpdir + f"/decompressed.{uncompressed_file_type}"
                    with gzip.open(fn_matrix, "rb") as input_f, open(tmppth, "wb") as output_f:
                        shutil.copyfileobj(input_f, output_f)
                    mat = scipy.io.mmread(tmppth)
                matrix.append(mat.T)

    adatas = []
    for i in range(len(matrix)):
        temp = ad.AnnData(X=matrix[i], obs=barcodes[i], var=genes[i])
        adatas.append(temp)
    adata = ad.concat(adatas, merge='first')

    return adata
