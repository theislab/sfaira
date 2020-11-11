import anndata
try:
    from kipoi.model import BaseModel
except ImportError:
    BaseModel = None
import numpy as np
import pandas as pd
import os
from typing import List, Union

from .external import EstimatorKerasEmbedding, EstimatorKerasCelltype
from .model_zoo import ModelZooEmbedding, ModelZooCelltype


class UserInterface:
    """
    This class performs data set handling and coordinates estimators for the different model types.
    Example code to obtain a UMAP embedding plot of the embedding created from your data with cell-type labels:
    ```
    import sfaira.api as sfaira
    import anndata
    import scanpy

    # initialise your sfaira instance with a model lookuptable.
    # instead of setting `custom_repo` when initialising the UI you can also use `sfaira_repo=True` to use public weights
    ui = sfaira.ui.UserInterface(custom_repo="/path/to/local/repo/folder/or/zenodo/repo/URL", sfaira_repo=False)
    ui.zoo_embedding.set_latest(species, organ, model_type, organisation, model_topology)
    ui.zoo_celltype.set_latest(species, organ, model_type, organisation, model_topology)
    ui.load_data(anndata.read("/path/to/file.h5ad"))  # load your dataset into sfaira
    ui.load_model_embedding()
    ui.load_model_celltype()
    ui.compute_all()
    adata = ui.data
    scanpy.pp.neighbors(adata, use_rep="X_sfaira")
    scanpy.tl.umap(adata)
    scanpy.pl.umap(adata, color="celltype_sfaira", show=True, save="UMAP_sfaira.png")
    ```
    """

    estimator_embedding: Union[EstimatorKerasEmbedding, None]
    estimator_celltype: Union[EstimatorKerasCelltype, None]
    model_kipoi_embedding: Union[None]
    model_kipoi_celltype: Union[BaseModel, None]
    zoo_embedding: Union[ModelZooEmbedding, None]
    zoo_celltype: Union[ModelZooCelltype, None]
    data: Union[anndata.AnnData]
    model_lookuptable: Union[pd.DataFrame, None]

    def __init__(
            self,
            custom_repo: Union[list, str, None] = None,
            sfaira_repo: bool = False,
            cache_path: str = os.path.join('cache', '')
    ):
        self.model_kipoi_embedding = None
        self.model_kipoi_celltype = None
        self.estimator_embedding = None
        self.estimator_celltype = None
        self.use_sfaira_repo = sfaira_repo
        self.cache_path = os.path.join(cache_path, '')

        if sfaira_repo:  # check if public sfaira repository should be accessed
            self.model_lookuptable = self._load_lookuptable("https://sandbox.zenodo.org/record/647061/files/")   #TODO: this still points to zenodo sandbox

        if custom_repo:
            if isinstance(custom_repo, str):
                custom_repo = [custom_repo]

            for repo in custom_repo:
                if os.path.exists(repo) and not os.path.exists(os.path.join(repo, 'model_lookuptable.csv')):
                    self.write_lookuptable(repo)

                if hasattr(self, 'model_lookuptable'):
                    self.model_lookuptable = self.model_lookuptable.append(
                        self._load_lookuptable(repo)
                    )
                else:
                    self.model_lookuptable = self._load_lookuptable(repo)

        else:
            if not sfaira_repo:
                raise ValueError("please either provide a custom folder/repository with model weights or specify "
                                 "`sfaira_repo=True` to access the public weight repository")

        self.zoo_embedding = ModelZooEmbedding(self.model_lookuptable)
        self.zoo_celltype = ModelZooCelltype(self.model_lookuptable)

    def _load_lookuptable(
            self,
            repo_path: str
    ):
        """
        checks if there is a csv file that lists the model_id and path of models in the directory
        returns model_lookuptable that connects model_id with the link to the model

        :param repo_path:
        :return: model_lookuptable
        """
        model_lookuptable = pd.read_csv(os.path.join(repo_path, 'model_lookuptable.csv'), header=0, index_col=0)

        # check for any duplicated model_ids
        if hasattr(self, 'model_lookuptable'):
            old_ids = self.model_lookuptable['model_id'].values
            new_ids = model_lookuptable['model_id'].values
            if any(ix in old_ids for ix in new_ids):
                ixs = new_ids[[ix in old_ids for ix in new_ids]]
                raise ValueError('Model ids ' + ','.join(ixs) + ' already exist in a loaded lookuptable.'
                                                                'Please remove any duplicate model_ids')
        return model_lookuptable

    def write_lookuptable(
            self,
            repo_path: str
    ):
        """
        :param repo_path:
        :return:
        """
        import hashlib

        file_names = []
        file_paths = []
        md5 = []
        for subdir, dirs, files in os.walk(repo_path):
            for file in files:
                if os.path.isfile(os.path.join(subdir, file)) and (
                        file.endswith('_weights.h5') or file.endswith('_weights.data-00000-of-00001')) and (
                        file.startswith('embedding') or file.startswith('celltype')):
                    file_paths.append(subdir)
                    file_names.append(file)
                    with open(os.path.join(subdir, file), 'rb') as f:
                        md5.append(hashlib.md5(f.read()).hexdigest())
        s = [i.split('_')[0:7] for i in file_names]
        ids = ['_'.join(i) for i in s]

        if ids:
            pd.DataFrame(
                list(zip(ids, file_paths, md5)),
                columns=['model_id', 'model_path', 'md5']
            ).sort_values('model_id').to_csv(os.path.join(repo_path, 'model_lookuptable.csv'))
        else:
            raise ValueError(f'No model weights found in {repo_path} '
                             'Weights need to have .h5 or .data-00000-of-00001 extension'
                             'to be recognised'
                             )

    def load_data(
            self,
            data: anndata.AnnData
    ):
        """

        :return:
        """
        self.data = data

    def filter_cells(self):
        """
        Filters cells with a basic pre-defined filter.

        :return:
        """
        # call cell_filter()
        raise NotImplementedError()

    def load_model_embedding(self):
        """
        Initialise embedding model and load parameters from public parameter repository.

        Loads model defined in self.model_id.

        :return: Model ID loaded.
        """
        assert self.zoo_embedding.model_id is not None, "choose embedding model first"
        model_dir = self.model_lookuptable.model_path[self.model_lookuptable.model_id == self.zoo_embedding.model_id].iloc[0]
        model_dir = self.path.join(model_dir, '')
        md5 = self.model_lookuptable.md5[self.model_lookuptable.model_id == self.zoo_embedding.model_id].iloc[0]
        self.estimator_embedding = EstimatorKerasEmbedding(
            data=self.data,
            model_dir=model_dir,
            model_id=self.zoo_embedding.model_id,
            species=self.zoo_embedding.species,
            organ=self.zoo_embedding.organ,
            model_type=self.zoo_embedding.model_type,
            model_topology=self.zoo_embedding.model_topology,
            weights_md5=md5,
            cache_path=self.cache_path
        )
        self.estimator_embedding.init_model()
        self.estimator_embedding.load_pretrained_weights()

    def load_model_celltype(self):
        """
        Initialise cell type model and load parameters from public parameter repository.

        Loads model defined in self.model_id.

        :return: Model ID loaded.
        """
        assert self.zoo_celltype.model_id is not None, "choose cell type model first"
        model_dir = self.model_lookuptable.model_path[self.model_lookuptable.model_id == self.zoo_celltype.model_id].iloc[0]
        model_dir = self.path.join(model_dir, '')
        md5 = self.model_lookuptable.md5[self.model_lookuptable.model_id == self.zoo_celltype.model_id].iloc[0]
        self.estimator_celltype = EstimatorKerasCelltype(
            data=self.data,
            model_dir=model_dir,
            model_id=self.zoo_celltype.model_id,
            species=self.zoo_celltype.species,
            organ=self.zoo_celltype.organ,
            model_type=self.zoo_celltype.model_type,
            model_topology=self.zoo_celltype.model_topology,
            weights_md5=md5,
            cache_path=self.cache_path
        )
        self.estimator_celltype.init_model()
        self.estimator_celltype.load_pretrained_weights()

    def _adata_write_celltype(
            self,
            labels: List[str],
            key: str
    ):
        """
        Writes a list of cell type labels into the column of adata.obs indicated
        :return:
        """
        self.data.obs[key] = [self.zoo_celltype.celltypes[i][0] for i in np.argmax(labels, axis=1)]

    def _adata_write_embedding(
            self,
            embedding: List[str],
            key: str
    ):
        """
        Writes the embedding matrix into adata.obsm with the key indicated.
        :return:
        """
        self.data.obsm[key] = embedding

    def _adata_write_denoised_data(
            self,
            denoised_data: np.ndarray,
            key: str
    ):
        """
        Writes the denoised expression matrix into adata.obsm with the key indicated.
        :return:
        """
        self.data.layers[key] = denoised_data

    def compute_celltype(self):
        """
        Run local cell type prediction model and add predictions to adata.obs.

        :return:
        """
        if self.zoo_celltype is not None:
            self._adata_write_celltype(
                labels=self.estimator_celltype.predict(),
                key="celltype_sfaira"
            )
        else:
            raise ValueError("celltype zoo has to be set before local model can be run.")

    def compute_embedding(self):
        """
        Run local embedding prediction model and add embedding to adata.obsm.

        :return:
        """
        if self.zoo_embedding is not None:
            self._adata_write_embedding(
                embedding=self.estimator_embedding.predict_embedding(),
                key="X_sfaira"
            )
        else:
            raise ValueError("embedding zoo has to be set before local model can be run.")

    def compute_all(self):
        """
        Run local cell type prediction and embedding models and add results of both to adata.

        :return:
        """
        self.compute_embedding()
        self.compute_celltype()

    def compute_denoised_expression(self):
        """
        Run local embedding prediction model and add denoised expression to new adata layer.

        :return:
        """
        if self.zoo_embedding is not None:
            self._adata_write_denoised_data(
                denoised_data=self.estimator_embedding.predict(),
                key="denoised_sfaira"
            )
        else:
            raise ValueError("embedding zoo has to be set before local model can be run.")

    def compute_celltype_kipoi(self):
        """
        Run executable cell type prediction model from kipoi_experimental and add prediction to adata.obs.

        :return:
        """
        if self.zoo_celltype is not None:
            self.model_kipoi_celltype = self.zoo_celltype.get_kipoi_model()
            self._adata_write_celltype(
                labels=self.model_kipoi_celltype.pipeline.predict(dict(adata=self.data)),
                key="celltype_sfaira"
            )
        else:
            raise ValueError("celltype zoo has to be set before kipoi_experimental model can be run.")

    def compute_embedding_kipoi(self):
        """
        Run executable embedding prediction model from kipoi_experimental and add embedding to adata.obsm.

        :return:
        """
        if self.zoo_embedding is not None:
            self.model_kipoi_embedding = self.zoo_embedding.get_kipoi_model()
            self._adata_write_embedding(
                embedding=self.model_kipoi_embedding.pipeline.predict_embedding(dict(adata=self.data)),
                key="X_sfaira"
            )
        else:
            raise ValueError("embedding zoo has to be set before kipoi_experimental model can be run.")

    def compute_all_kipoi(self):
        """
        Run executable cell type prediction and embedding models from kipoi_experimental and add results to adata.

        :return:
        """
        self.compute_embedding_kipoi()
        self.compute_celltype_kipoi()

    def compute_denoised_expression_kipoi(self):
        """
        Run executable embedding prediction model from kipoi_experimental and add denoised expression to adata layer.

        :return:
        """
        if self.zoo_embedding is not None:
            self.model_kipoi_embedding = self.zoo_embedding.get_kipoi_model()
            self._adata_write_denoised_data(
                denoised_data=self.model_kipoi_embedding.pipeline.predict(dict(adata=self.data)),
                key="denoised_sfaira"
            )
        else:
            raise ValueError("embedding zoo has to be set before local model can be run.")

    def celltype_summary(self):
        """
        Return type with frequencies of predicted cell types.

        :return:
        """
        return self.data.obs['celltype_sfaira'].value_counts()

    def get_references(self):
        """
        Return papers to cite when using the embedding model.

        Collects references from the estimators of each model type.

        :return:
        """
        return self.estimator_embedding.get_citations()
