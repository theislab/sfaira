import anndata
import numpy as np
import pandas as pd
import os
from typing import List, Union
import warnings

from sfaira.data import DatasetInteractive
from sfaira.estimators import EstimatorKerasEmbedding, EstimatorKerasCelltype
from sfaira.interface.model_zoo import ModelZoo
from sfaira.consts import AdataIdsSfaira, AdataIds


class UserInterface:
    """
    This class performs data set handling and coordinates estimators for the different model types.
    Example code to obtain a UMAP embedding plot of the embedding created from your data with cell-type labels:
    ```
    import sfaira
    import anndata
    import scanpy

    # initialise your sfaira instance with a model lookuptable.
    # instead of setting `custom_repo` when initialising the UI you can also use `sfaira_repo=True` to use public weights
    ui = sfaira.ui.UserInterface(custom_repo="/path/to/local/repo/folder/or/zenodo/repo/URL", sfaira_repo=False)
    ui.zoo_embedding.set_latest(organism, organ, model_type, organisation, model_topology)
    ui.zoo_celltype.set_latest(organism, organ, model_type, organisation, model_topology)
    ui.load_data(anndata.read("/path/to/file.h5ad"))  # load your dataset into sfaira
    ui.load_model_embedding()
    ui.load_model_celltype()
    ui.predict_all()
    adata = ui.data.adata
    scanpy.pp.neighbors(adata, use_rep="X_sfaira")
    scanpy.tl.umap(adata)
    scanpy.pl.umap(adata, color="celltypes_sfaira", show=True, save="UMAP_sfaira.png")
    ```
    """

    estimator_embedding: Union[EstimatorKerasEmbedding, None]
    estimator_celltype: Union[EstimatorKerasCelltype, None]
    zoo_embedding: Union[ModelZoo, None]
    zoo_celltype: Union[ModelZoo, None]
    data: Union[DatasetInteractive, None]
    model_lookuptable: Union[pd.DataFrame, None]
    adata_ids: AdataIds

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
        self.data = None
        self.use_sfaira_repo = sfaira_repo
        self.cache_path = os.path.join(cache_path, '')
        self.adata_ids = AdataIdsSfaira()

        if sfaira_repo:  # check if public sfaira repository should be accessed
            self.model_lookuptable = self._load_lookuptable("https://zenodo.org/record/4304660/files/")

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

        self.zoo_embedding = ModelZoo(model_lookuptable=self.model_lookuptable, model_class="embedding")
        self.zoo_celltype = ModelZoo(model_lookuptable=self.model_lookuptable, model_class="celltype")

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
        model_paths = []
        file_paths = []
        md5 = []
        for subdir, dirs, files in os.walk(repo_path):
            for file in files:
                if os.path.isfile(os.path.join(subdir, file)) and (
                        file.endswith('.h5') or file.endswith('.data-00000-of-00001')) and (
                        file.startswith('embedding_') or file.startswith('celltype_')):
                    model_paths.append(os.path.join(subdir, ""))
                    file_paths.append(os.path.join(subdir, file))
                    file_names.append(file)
                    with open(os.path.join(subdir, file), 'rb') as f:
                        md5.append(hashlib.md5(f.read()).hexdigest())
        ids = ['_'.join(i.split('_')[0:3]) for i in file_names]
        ids_cleaned = [i.replace('.h5', '').replace('.data-00000-of-00001', '') for i in ids]  # remove file extensions

        if ids:
            pd.DataFrame(
                list(zip(ids_cleaned, model_paths, file_paths, md5)),
                columns=['model_id', 'model_path', 'model_file_path', 'md5']
            )\
                .sort_values('model_id')\
                .reset_index(drop=True)\
                .to_csv(os.path.join(repo_path, 'model_lookuptable.csv'))
        else:
            raise ValueError(f'No model weights found in {repo_path} '
                             'Weights need to have .h5 or .data-00000-of-00001 extension'
                             'to be recognised'
                             )

    def deposit_zenodo(
            self,
            zenodo_access_token: str,
            title: str,
            authors: list,
            description: str,
            metadata: dict = {},
            publish: bool = False,
            sandbox: bool = False
    ):
        """
        Deposit all models in model lookup table on Zenodo. If publish is set to false, files will be uploaded to a
        deposition draft, which can be further edited (additional metadata, files etc.). Returns the DOI link if
        publish=True or a link to the deposition draft if publish=False.

        :param zenodo_access_token: Your personal Zenodo API access token. Create one here: https://zenodo.org/account/settings/applications/tokens/new/
        :param title: Title of the Zenodo deposition
        :param authors: List of dicts, where each dict defines one author (dict keys:
         name: Name of creator in the format "Family name, Given names",
         affiliation: Affiliation of creator (optional), orcid: ORCID identifier of creator (optional),
         gnd: GND identifier of creator (optional)
        :param description: Description of the Zenodo deposition.
        :param metadata: Dictionary with further metadata attributes of the deposit.
         See the Zenodo API refenrece for accepted keys: https://developers.zenodo.org/#representation
        :param publish: Set this to True to directly publish the weights on Zenodo.
         When set to False a draft will be created, which can be edited in the browser before publishing.
        :param sandbox: If True, use the Zenodo testing platform at https://sandbox.zenodo.org for your deposition.
         We recommend testing your upload with sandbox first as depositions cannot be deleted from the main Zenodo platfowm once created.
        """

        import requests
        import json
        headers = {"Content-Type": "application/json"}
        params = {'access_token': zenodo_access_token}
        sandbox = 'sandbox.' if sandbox else ''

        # Verify access token
        r = requests.get(f'https://{sandbox}zenodo.org/api/deposit/depositions', params=params)
        if r.status_code != 200:
            raise ValueError(
                "Your Zenodo access token was not accepted by the API. Please provide a valid access token.")

        # Create empty deposition
        r = requests.post(f'https://{sandbox}zenodo.org/api/deposit/depositions',
                          params=params,
                          json={},
                          headers=headers)

        # Obtain bucket URL and deposition ID
        bucket_url = r.json()["links"]["bucket"]
        deposition_id = r.json()['id']

        # Loop over files in model lookup table and upload them one by one
        for i, weight_path in enumerate(self.model_lookuptable['model_file_path']):
            filename = os.path.basename(weight_path)
            with open(weight_path, "rb") as fp:
                r = requests.put(
                    f"{bucket_url}/{filename}",
                    data=fp,
                    params=params,
                )
            # Verify checksum after upload
            if r.json()['checksum'][4:] != self.model_lookuptable['md5'][i]:
                warnings.warn(f"The md5 checksum in your model_lookuptable for {self.model_lookuptable['model_id'][i]} "
                              f"does not match the md5 checksum of the uploaded file.")

        # Add model lookup table to zenodo
        df = self.model_lookuptable.copy()
        df['model_path'] = f"https://{sandbox}zenodo.org/record/{deposition_id}/files/"
        df['model_file_path'] = [f"https://{sandbox}zenodo.org/record/{deposition_id}/files/{os.path.basename(f)}" for f
                                 in self.model_lookuptable['model_file_path']]
        r = requests.put(
            f"{bucket_url}/model_lookuptable.csv",
            data=df.to_csv(),
            params=params,
        )

        # Add metadata
        meta_core = {
            'title': title,
            'creators': authors,
            'description': description,
            'license': 'cc-by-4.0',
            'upload_type': 'dataset',
            'access_right': 'open'
        }
        meta = {**meta_core, **metadata}
        r = requests.put(f'https://{sandbox}zenodo.org/api/deposit/depositions/{deposition_id}',
                         params=params,
                         data=json.dumps({
                             'metadata': meta
                         }),
                         headers=headers)

        if not publish:
            print(f"Zenodo deposition draft has been created: {r.json()['links']['latest_draft_html']}")
            return r.json()['links']['latest_draft_html']
        else:
            # Publish the deposition
            r = requests.post(f'https://{sandbox}zenodo.org/api/deposit/depositions/{deposition_id}/actions/publish',
                              params=params)
            if r.status_code == 202:
                if sandbox:
                    print(f"Weights referenced in model_lookuptable have been sucessfully published on Zenodo: "
                          f"{r.json()['links']['latest_html']}")
                    return r.json()['links']['latest_html']
                else:
                    print(f"Weights referenced in model_lookuptable have been sucessfully published on Zenodo: "
                          f"{r.json()['links']['conceptdoi']}")
                    return r.json()['links']['conceptdoi']
            else:
                try:
                    m = r.json()['message']
                except KeyError:
                    m = f"Submission failed with html status code {r.status_code}"
                raise ValueError(m)

    def load_data(
            self,
            data: anndata.AnnData,
            gene_symbol_col: Union[str, None] = None,
            gene_ens_col: Union[str, None] = None,
            obs_key_celltypes: Union[str, None] = None,
            class_maps: dict = {},
    ):
        """
        Loads the provided AnnData object into sfaira.

        If genes in the provided AnnData object are annotated as gene symbols,
         please provide the name of the corresponding var column (or 'index') through the gene_symbol_col argument.
        If genes in the provided AnnData object are annotated as ensembl ids,
         please provide the name of the corresponding var column (or 'index') through the gene_ens_col argument.
        You need to provide at least one of the two.
        :param data: AnnData object to load
        :param gene_symbol_col: Var column name (or 'index') which contains gene symbols
        :param gene_ens_col: ar column name (or 'index') which contains ensembl ids
        :param obs_key_celltypes: .obs column name which contains cell type labels.
        :param class_maps: Cell type class maps.
        """
        if self.zoo_embedding.model_organism is not None:
            organism = self.zoo_embedding.model_organism
            organ = self.zoo_embedding.model_organ
        elif self.zoo_celltype.model_organism is not None:
            organism = self.zoo_celltype.model_organism
            organ = self.zoo_celltype.model_organ
        else:
            raise ValueError("Please first set which model_id to use via the model zoo before loading the data")

        if gene_ens_col is None and gene_symbol_col is None:
            raise ValueError("Please provide either the gene_ens_col or the gene_symbol_col argument.")

        self.data = DatasetInteractive(
            data=data,
            organism=organism,
            organ=organ,
            gene_symbol_col=gene_symbol_col,
            gene_ens_col=gene_ens_col,
            obs_key_celltypes=obs_key_celltypes,
            class_maps=class_maps,
        )
        # Align to correct featurespace
        self.data.streamline_features(
            match_to_reference=self.zoo_embedding.topology_container.gc.assembly,
            subset_genes_to_type=list(set(self.zoo_embedding.topology_container.gc.biotype))
        )

    def load_model_embedding(self):
        """
        Initialise embedding model and load parameters from public parameter repository.

        Loads model defined in self.model_id.

        :return: Model ID loaded.
        """
        assert self.zoo_embedding.model_id is not None, "choose embedding model first"
        model_dir = self.model_lookuptable["model_file_path"].loc[self.model_lookuptable["model_id"] == self.zoo_embedding.model_id].iloc[0]
        md5 = self.model_lookuptable["md5"].loc[self.model_lookuptable["model_id"] == self.zoo_embedding.model_id].iloc[0]
        self.estimator_embedding = EstimatorKerasEmbedding(
            data=self.data.adata,
            model_dir=model_dir,
            model_id=self.zoo_embedding.model_id,
            model_topology=self.zoo_embedding.topology_container,
            weights_md5=md5,
            cache_path=self.cache_path,
            adata_ids=self.adata_ids
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
        model_dir = self.model_lookuptable["model_file_path"].loc[self.model_lookuptable["model_id"] == self.zoo_celltype.model_id].iloc[0]
        md5 = self.model_lookuptable["md5"].loc[self.model_lookuptable["model_id"] == self.zoo_celltype.model_id].iloc[0]
        self.estimator_celltype = EstimatorKerasCelltype(
            data=self.data.adata,
            model_dir=model_dir,
            model_id=self.zoo_celltype.model_id,
            model_topology=self.zoo_celltype.topology_container,
            weights_md5=md5,
            cache_path=self.cache_path,
            remove_unlabeled_cells=False,
            adata_ids=self.adata_ids
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
        self.data.adata.obs[key] = [self.zoo_celltype.celltypes[i] for i in np.argmax(labels, axis=1)]
        self.data.adata.obs[key] = self.data.adata.obs[key].astype('category')

    def _adata_write_embedding(
            self,
            embedding: List[str],
            key: str
    ):
        """
        Writes the embedding matrix into adata.obsm with the key indicated.
        :return:
        """
        self.data.adata.obsm[key] = embedding

    def _adata_write_denoised_data(
            self,
            denoised_data: np.ndarray,
            key: str
    ):
        """
        Writes the denoised expression matrix into adata.obsm with the key indicated.
        :return:
        """
        self.data.adata.layers[key] = denoised_data

    def predict_celltypes(self):
        """
        Run local cell type prediction model and add predictions to adata.obs.

        :return:
        """
        if self.zoo_celltype is not None:
            self._adata_write_celltype(
                labels=self.estimator_celltype.predict(),
                key="celltypes_sfaira"
            )
        else:
            raise ValueError("celltype zoo has to be set before local model can be run.")

    def predict_embedding(self):
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

    def predict_all(self):
        """
        Run local cell type prediction and embedding models and add results of both to adata.

        :return:
        """
        self.predict_embedding()
        self.predict_celltypes()

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

    def celltype_summary(self):
        """
        Return type with frequencies of predicted cell types.

        :return:
        """
        assert "celltypes_sfaira" in self.data.adata.obs.keys(), \
            "Column celltypes_sfaira not found in the data. Please run UserInterface.predict_celltypes() first."
        return self.data.adata.obs['celltypes_sfaira'].value_counts()
