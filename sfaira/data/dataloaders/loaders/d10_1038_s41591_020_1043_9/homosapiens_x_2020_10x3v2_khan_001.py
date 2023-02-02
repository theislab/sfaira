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
    fn = os.path.join(data_dir, "GSE145122_RAW.tar")
    barcodes = []
    genes = []
    matrix = []
    with tarfile.open(fn) as tar:
        for member in tar.getmembers():
            temp = member.name
            if re.search('barcodes', temp):
                bar = pd.read_csv(tar.extractfile(member), delimiter='\t', compression="gzip", header=None, index_col=0)
                bar.index.name = None
                sample = member.name.split("_")[1] + "_" + member.name.split("_")[2]
                bar['sample'] = sample
                if sample == "Control_1":
                    bar['cell_line'] = "0524-1"
                    bar['sex'] = "female"
                    bar['disease'] = "none"
                    bar['diagnosis'] = "none"
                    bar['full_scale_IQ'] = "120"
                    bar['subject_age'] = "27 ys"
                elif sample == "Control_2":
                    bar['cell_line'] = "0307-1"
                    bar['sex'] = "male"
                    bar['disease'] = "none"
                    bar['diagnosis'] = "none"
                    bar['full_scale_IQ'] = "N/A"
                    bar['subject_age'] = "27 ys"
                elif sample == "Patient_1":
                    bar['cell_line'] = "7958-3"
                    bar['sex'] = "female"
                    bar['disease'] = "22q11DS"
                    bar['diagnosis'] = "Schizophrenia, Depressive Disorder Not Otherwise Specified"
                    bar['full_scale_IQ'] = "56"
                    bar['subject_age'] = "25 ys"
                else:
                    bar['cell_line'] = "6303-5"
                    bar['sex'] = "female"
                    bar['disease'] = "22q11DS"
                    bar['diagnosis'] = "no psychiatric diagnoses "
                    bar['full_scale_IQ'] = "65"
                    bar['subject_age'] = "22 ys"
                barcodes.append(bar)
            elif re.search('features', temp):
                gen = pd.read_csv(tar.extractfile(member), delimiter='\t', compression='gzip',
                                  index_col=0, header=None, names=["ensembl_gene_id", "gene_id", "status"])
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
    adata.obs['organoid_age_days'] = "83"
    del( adata.var['status'])

    return adata
