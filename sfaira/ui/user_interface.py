import anndata
import numpy as np
import os
import re
import pandas as pd
import pickle
from typing import List, Union
import warnings
import time

from sfaira import settings
from sfaira.consts import AdataIdsSfaira, AdataIds, OCS
from sfaira.data import DatasetInteractive
from sfaira.estimators.keras.base import EstimatorKerasEmbedding, EstimatorKerasCelltype
from sfaira.ui.model_zoo import ModelZoo
from sfaira.versions.topologies import TopologyContainer


class UserInterface:
    """
    This class performs data set handling and coordinates estimators for the different model types.
    Example code to obtain a UMAP embedding plot of the embedding created from your data with cell-type labels::

        import sfaira
        import anndata
        import scanpy

        # initialise your sfaira instance with a model lookuptable.
        ui = sfaira.ui.UserInterface(custom_repo="/path/to/local/repo/folder/or/zenodo/repo/URL", sfaira_repo=False)
        ui.zoo_embedding.model_id = 'embedding_human-blood-ae-0.2-0.1_theislab'  # pick desired model here
        ui.zoo_celltype.model_id = 'celltype_human-blood-mlp-0.1.3-0.1_theislab'  # pick desired model here
        ui.load_data(anndata.read("/path/to/file.h5ad"), gene_symbol_col='index', gene_ens_col='gene_ids')
        ui.load_model_embedding()
        ui.load_model_celltype()
        ui.predict_all()
        adata = ui.data.adata
        scanpy.pp.neighbors(adata, use_rep="X_sfaira")
        scanpy.tl.umap(adata)
        scanpy.pl.umap(adata, color="celltypes_sfaira", show=True, save="UMAP_sfaira.png")

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
        self.estimator_embedding = None
        self.estimator_celltype = None
        self.data = None
        self.use_sfaira_repo = sfaira_repo
        self.cache_path = os.path.join(cache_path, '')
        self.adata_ids = AdataIdsSfaira()

        if sfaira_repo:  # check if public sfaira repository should be accessed
            self.model_lookuptable = self._load_lookuptable(settings.sfaira_repo_url)

        if custom_repo:
            if isinstance(custom_repo, str):
                custom_repo = [custom_repo]

            for repo in custom_repo:
                if not os.path.exists(repo):
                    raise OSError(f"provided repo directory does not exist, please create it first: {repo}")

                if not os.path.exists(os.path.join(repo, 'model_lookuptable.csv')):
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
            update_existing_deposition: Union[None, str] = None,
            publish: bool = False,
            sandbox: bool = False,
            deposit_topologies: bool = True
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
        :param update_existing_deposition: If None, a new deposition will be created.
        If an existing deposition ID is provided as a sting, than this deposition will be updated with a new version.
        :param publish: Set this to True to directly publish the weights on Zenodo.
         When set to False a draft will be created, which can be edited in the browser before publishing.
        :param sandbox: If True, use the Zenodo testing platform at https://sandbox.zenodo.org for your deposition.
         We recommend testing your upload with sandbox first as depositions cannot be deleted from the main Zenodo platform once created.
         :param deposit_topologies: If true, an associated topology file for every weights file will be uploaded to zenodo.
         The naming format for the topology files is <model_id>_topology.pickle
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

        if update_existing_deposition is None:
            # Create empty deposition
            r = requests.post(f'https://{sandbox}zenodo.org/api/deposit/depositions',
                              params=params,
                              json={},
                              headers=headers)
            # Obtain bucket URL and deposition ID
            bucket_url = r.json()["links"]["bucket"]
            deposition_id = r.json()['id']
        else:
            update_existing_deposition = str(update_existing_deposition) if isinstance(update_existing_deposition, int)\
                else update_existing_deposition
            # Create a new version of the existing deposition
            r = requests.post(
                f'https://{sandbox}zenodo.org/api/deposit/depositions/{update_existing_deposition}/actions/newversion',
                params=params)
            try:
                deposition_id = r.json()["links"]["latest_draft"].split("/")[-1]
            except json.decoder.JSONDecodeError:
                time.sleep(10)
                r = requests.post(
                    f'https://{sandbox}zenodo.org/api/deposit/depositions/{update_existing_deposition}/actions/newversion',
                    params=params)
                deposition_id = r.json()["links"]["latest_draft"].split("/")[-1]
            if r.status_code != 201:
                raise ValueError(
                    f"A new version of deposition {update_existing_deposition} could not be created, "
                    f"please make sure your API key is associated with the account that owns this deposition.")
            r = requests.get(f'https://{sandbox}zenodo.org/api/deposit/depositions/{deposition_id}', params=params)
            bucket_url = r.json()["links"]["bucket"]

            # Delete all existing files from new version
            r_files = requests.get(f'https://{sandbox}zenodo.org/api/deposit/depositions/{deposition_id}/files',
                                   params=params)
            while len(r_files.json()) > 0:
                for file_dict in r_files.json():
                    requests.delete(
                        f'https://{sandbox}zenodo.org/api/deposit/depositions/{deposition_id}/files/{file_dict["id"]}',
                        params=params)
                r_files = requests.get(f'https://{sandbox}zenodo.org/api/deposit/depositions/{deposition_id}/files',
                                       params=params)
                while isinstance(r_files.json(), dict):
                    print("Pausing due to Zenodo API rate limitng")
                    time.sleep(10)
                    r_files = requests.get(f'https://{sandbox}zenodo.org/api/deposit/depositions/{deposition_id}/files',
                                           params=params)

        # Loop over files in model lookup table and upload them one by one
        for i, weight_path in enumerate(self.model_lookuptable['model_file_path']):
            basepath, filename_weights = os.path.split(weight_path)
            with open(weight_path, "rb") as fp:
                r = requests.put(
                    f"{bucket_url}/{filename_weights}",
                    data=fp,
                    params=params,
                )
            while r.status_code != 200:
                print(f"Upload of {weight_path} was not successful (status code {r.status_code}), retrying")
                time.sleep(10)
                with open(weight_path, "rb") as fp:
                    r = requests.put(
                        f"{bucket_url}/{filename_weights}",
                        data=fp,
                        params=params,
                    )
            # Verify checksum after upload
            if r.json()['checksum'][4:] != self.model_lookuptable['md5'][i]:
                warnings.warn(f"The md5 checksum in your model_lookuptable for {self.model_lookuptable['model_id'][i]} "
                              f"does not match the md5 checksum of the uploaded file.")
            if deposit_topologies:  # Deposit associated topology file
                filename_topology = ".".join(filename_weights.split(".")[:-1])
                filename_topology = re.sub(r"_weights$", "", filename_topology)
                filename_topology += "_topology.pickle"
                topology_path = os.path.join(basepath, filename_topology)
                assert os.path.isfile(topology_path), f"topology file {topology_path} not found. " \
                                                      f"consider deactivating the deposition of topology files."
                with open(topology_path, "rb") as fp:
                    r = requests.put(
                        f"{bucket_url}/{filename_topology}",
                        data=fp,
                        params=params,
                    )

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
        """
        if self.zoo_embedding.model_organism is not None and self.zoo_celltype.model_organism is not None:
            assert self.zoo_embedding.model_organism == self.zoo_celltype.model_organism, \
                "Model ids set for embedding and celltype model need to correspond to the same organism"
            assert self.zoo_embedding.model_organ == self.zoo_celltype.model_organ, \
                "Model ids set for embedding and celltype model need to correspond to the same organ"
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

        # handle organ names with stripped spaces
        if organ not in OCS.organ.node_names:
            organ_dict = {i.replace(" ", ""): i for i in OCS.organ.node_names}
            assert organ in organ_dict, f"Organ {organ} is not a valid nodename in the UBERON organ ontology"
            organ = {i.replace(" ", ""): i for i in OCS.organ.node_names}[organ]

        self.data = DatasetInteractive(
            data=data,
            feature_symbol_col=gene_symbol_col,
            feature_id_col=gene_ens_col,
        )
        self.data.organism = organism
        self.data.organ = organ
        self.data.cell_type_obs_key = obs_key_celltypes
        # Align to correct featurespace
        self.data.streamline_features(
            match_to_release=self.zoo_embedding.topology_container.gc.release,
            subset_genes_to_type=list(set(self.zoo_embedding.topology_container.gc.biotype))
        )
        # Transfer required metadata from the Dataset instance to the adata object
        self.data.streamline_metadata(
            clean_obs=False,
            clean_var=True,
            clean_uns=False,
            clean_obs_names=False,
        )

    def _load_topology_dict(self, model_weights_file) -> dict:
        topology_filepath = ".".join(model_weights_file.split(".")[:-1])
        topology_filepath = re.sub(r"_weights$", "", topology_filepath)
        topology_filepath += "_topology.pickle"
        if topology_filepath.startswith('http'):
            # Download into cache if file is on a remote server.
            if not os.path.exists(self.cache_path):
                os.makedirs(self.cache_path)
            import urllib.request
            from urllib.error import HTTPError
            try:
                urllib.request.urlretrieve(
                    topology_filepath,
                    os.path.join(self.cache_path, os.path.basename(topology_filepath))
                )
                topology_filepath = os.path.join(self.cache_path, os.path.basename(topology_filepath))
            except HTTPError:
                raise FileNotFoundError(f"cannot find remote topology file {topology_filepath}")
        with open(topology_filepath, "rb") as f:
            topology = pickle.load(f)
        return topology

    def load_model_embedding(self):
        """
        Initialise embedding model and load parameters from public parameter repository.

        Loads model defined in self.model_id.

        :return: Model ID loaded.
        """
        assert self.zoo_embedding.model_id is not None, "choose embedding model first"
        if self.zoo_celltype.topology_container.gc.release is not None:
            assert self.zoo_embedding.topology_container.gc.release == \
                self.zoo_celltype.topology_container.gc.release, \
                "genome assemblies defined in the topology " \
                "containers if the embedding and the celltype " \
                "prediction model are not equivalent " \
                f"({self.zoo_embedding.topology_container.gc.release} " \
                f"and {self.zoo_celltype.topology_container.gc.release} " \
                f"respectively, aborting.)"
        model_weights_file = self.model_lookuptable["model_file_path"].loc[self.model_lookuptable["model_id"] ==
                                                                           self.zoo_embedding.model_id].iloc[0]
        md5 = self.model_lookuptable["md5"].loc[self.model_lookuptable["model_id"] ==
                                                self.zoo_embedding.model_id].iloc[0]
        tc = TopologyContainer(
            topology=self._load_topology_dict(model_weights_file=model_weights_file),
            topology_id=self.zoo_embedding.topology_container.topology_id
        )
        self.estimator_embedding = EstimatorKerasEmbedding(
            data=self.data.adata,
            model_dir=model_weights_file,
            model_id=self.zoo_embedding.model_id,
            model_topology=tc,
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
        if self.zoo_embedding.topology_container.gc.release is not None:
            assert self.zoo_embedding.topology_container.gc.release == \
                self.zoo_celltype.topology_container.gc.release, \
                "genome assemblies defined in the topology " \
                "containers if the embedding and the celltype " \
                "prediction model are not equivalent " \
                f"({self.zoo_embedding.topology_container.gc.release} " \
                f"and {self.zoo_celltype.topology_container.gc.release} " \
                f"respectively, aborting.)"
        model_weights_file = self.model_lookuptable["model_file_path"].loc[self.model_lookuptable["model_id"] ==
                                                                           self.zoo_celltype.model_id].iloc[0]
        md5 = self.model_lookuptable["md5"].loc[self.model_lookuptable["model_id"] ==
                                                self.zoo_celltype.model_id].iloc[0]
        tc = TopologyContainer(
            topology=self._load_topology_dict(model_weights_file=model_weights_file),
            topology_id=self.zoo_celltype.topology_container.topology_id
        )
        self.estimator_celltype = EstimatorKerasCelltype(
            data=self.data.adata,
            model_dir=model_weights_file,
            model_id=self.zoo_celltype.model_id,
            model_topology=tc,
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
        key_id = key + "_id"
        self.data.adata.obs[key] = [self.estimator_celltype.ontology_names[i] for i in np.argmax(labels, axis=1)]
        self.data.adata.obs[key] = self.data.adata.obs[key].astype('category')
        self.data.adata.obs[key_id] = [self.estimator_celltype.ontology_ids[i] for i in np.argmax(labels, axis=1)]
        self.data.adata.obs[key_id] = self.data.adata.obs[key_id].astype('category')

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
