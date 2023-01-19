import anndata
import os
import scipy.sparse
import pandas as pd
import numpy as np
import h5py


def load(data_dir, **kwargs):
    f = h5py.File(os.path.join(data_dir, 'HumanFetalBrainPool.h5'), 'r')
    dset = f['shoji']

    X = scipy.sparse.csr_matrix(dset['Expression'][()], dtype=np.float32)
    obs = pd.DataFrame({'Age': np.array(dset['Age'], dtype=np.float32),
                        'CellClass': np.array(dset['CellClass'], dtype=str),
                        'CellCycleFraction': np.array(dset['CellCycleFraction'], dtype=np.float32),
                        'CellID': np.array(dset['CellID'], dtype=str),
                        'Chemistry': np.array(dset['Chemistry'], dtype=str),
                        'Clusters': np.array(dset['Clusters'], dtype=np.uint32),
                        'Donor': np.array(dset['Donor'], dtype=str),
                        'DoubletFlag': np.array(dset['DoubletFlag'], dtype=bool),
                        'DoubletScore': np.array(dset['DoubletScore'], dtype=np.float32),
                        'DropletClass': np.array(dset['DropletClass'], dtype=np.uint32),
                        'MitoFraction': np.array(dset['MitoFraction'], dtype=np.float32),
                        'NGenes': np.array(dset['NGenes'], dtype=np.uint32),
                        'PrevClusters': np.array(dset['PrevClusters'], dtype=np.uint32),
                        'Region': np.array(dset['Region'], dtype=str),
                        'SampleID': np.array(dset['SampleID'], dtype=str),
                        'Sex': np.array(dset['Sex'], dtype=str),
                        'Subdivision': np.array(dset['Subdivision'], dtype=str),
                        'Subregion': np.array(dset['Subregion'], dtype=str),
                        'Tissue': np.array(dset['Tissue'], dtype=str),
                        'TopLevelCluster': np.array(dset['TopLevelCluster'], dtype=np.uint32),
                        'TotalUMIs': np.array(dset['TotalUMIs'], dtype=np.uint32),
                        'UnsplicedFraction': np.array(dset['UnsplicedFraction'], dtype=np.float32),
                        'ValidCells': np.array(dset['ValidCells'], dtype=bool)},
                       index=np.array(dset['CellID'], dtype=str))

    var = pd.DataFrame({'Accession': np.array(dset['Accession'], dtype=str),
                        'Chromosome': np.array(dset['Chromosome'], dtype=str),
                        'End': np.array(dset['End'], dtype=str),
                        'Gene': np.array(dset['Gene'], dtype=str),
                        'GeneNonzeros': np.array(dset['GeneNonzeros'], dtype=np.uint32),
                        'GeneTotalUMIs': np.array(dset['GeneTotalUMIs'], dtype=np.uint32),
                        'SelectedFeatures': np.array(dset['SelectedFeatures'], dtype=bool),
                        'Start': np.array(dset['Start'], dtype=str),
                        'StdevExpression': np.array(dset['StdevExpression'], dtype=np.float32),
                        'ValidGenes': np.array(dset['ValidGenes'], dtype=bool)},
                       index=np.array(dset['Gene'], dtype=str))

    uns = {'AnnotationDefinition': np.array(dset['AnnotationDefinition'], dtype=str),
           'AnnotationDescription': np.array(dset['AnnotationDescription'], dtype=str),
           'AnnotationName': np.array(dset['AnnotationName'], dtype=str),
           'AnnotationPosterior': np.array(dset['AnnotationPosterior'], dtype=np.float32),
           'Linkage': np.array(dset['Linkage'], dtype=np.float32),
           'OverallTotalUMIs': np.array(dset['OverallTotalUMIs'], dtype=np.uint64),
           'Recipe': np.array(dset['Recipe'], dtype=str),
           'Species': np.array(dset['Species'], dtype=str)}

    obsm = {'X_umap': np.array(dset['Embedding'], dtype=np.float32),
            'X_factors': np.array(dset['Factors'], dtype=np.float32)}

    varm = {'Enrichment': np.array(dset['Enrichment'], dtype=np.float32).T,
            'Loadings': np.array(dset['Loadings'], dtype=np.float32),
            'Nonzeros': np.array(dset['Nonzeros'], dtype=np.uint64).T,
            'Trinaries': np.array(dset['Trinaries'], dtype=np.float32).T}

    f.close()

    adata = anndata.AnnData(X=X, obs=obs, var=var, obsm=obsm, varm=varm, uns=uns)

    return adata
