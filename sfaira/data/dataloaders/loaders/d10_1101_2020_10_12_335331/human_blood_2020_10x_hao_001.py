import tarfile 
import os
from typing import Union
import pandas as pd
import scipy.io
import gzip

import anndata
from sfaira.data import DatasetBase


# class Dataset(DatasetBase):
#     """
#     This is a dataloader for a multimodal PBMC data set (CITE-seq) from Hao et al (https://doi.org/10.1101/2020.10.12.335331).
    
#     :param path:
#     :param meta_path:
#     :param kwargs:
#     """
    
#     def __init__(
#             self,
#             **kwargs
#     ):
#         super().__init__(**kwargs)
#         self.download_url_data = "https://atlas.fredhutch.org/nygc/multimodal-pbmc/"  # download website(s) of data files
#         self.download_url_meta = "https://atlas.fredhutch.org/nygc/multimodal-pbmc/"  # download website(s) of meta data files

#         self.assay_sc = "10X sequencing" # (!) protocol used to sample data (e.g. smart-seq2)
#         self.author = "Hao, Yuhan"  # author (list) who sampled / created the data set
#         self.disease_obs_key = 'time'  # (!) whether sample represents a healthy organism
#         self.doi = "10.1101/2020.10.12.335331"  # doi of data set accompanying manuscript
#         self.normalisation = "raw"  # normalisation applied to raw data loaded (ideally counts, "raw")
#         self.organ = "blood"  # (!) organ
#         self.organism = "human"  # (!) species / organism
#         self.sample_source = 'primary_tissue' # (!) sub-tissue name, otherwise organ
#         self.year = 2020

#         self.development_stage = "Adult"  # (!) developmental stage of organism

#         self.var_symbol_col = "names"
#         self.cellontology_original_obs_key = "celltype.l3"

#         self.set_dataset_id(idx=1)


def load(data_dir, **kwargs):
    fn = os.path.join(data_dir, "GSE164378_RAW.tar")  #defined file in streamlined directory structure

    adatas = []        
    with tarfile.open(fn) as tar:
        samples = ['GSM5008737_RNA_3P', 'GSM5008738_ADT_3P'] # first sample is rna, second is protein data
        
        for sample in samples:
            print(sample)
            with gzip.open(tar.extractfile(sample + '-matrix.mtx.gz'), 'rb') as mm:
                print('Loading matrix')
                X = scipy.io.mmread(mm).T.tocsr()
            obs = pd.read_csv(tar.extractfile(sample + '-barcodes.tsv.gz'), compression='gzip', header=None, sep='\t', index_col=0)
            obs.index.name = None
            var = pd.read_csv(tar.extractfile(sample + '-features.tsv.gz'), compression='gzip', header=None, sep='\t').iloc[:, :1]
            var.columns = ['names']
            var.index = var['names'].values
            adata = anndata.AnnData(X=X, obs=obs, var=var)

            adata.var_names_make_unique()
            adatas.append(adata)
        
        tar.close()

    adata = adatas[0]
    protein = adatas[1]

    meta = pd.read_csv(os.path.join(data_dir, 'GSE164378_sc.meta.data_3P.csv.gz'), index_col=0)
    adata.obs = adata.obs.join(meta)

    adata.obsm['protein_expression'] = pd.DataFrame(protein.X.A, columns=protein.var_names, index=protein.obs_names)
    
    return adata

