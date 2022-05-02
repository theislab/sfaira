import datetime
from IPython.display import display, display_javascript, display_html
import json
import os
import pickle
import pydoc
from typing import Union
import uuid

from sfaira import settings
from sfaira.data.dataloaders.base import DatasetGroup, DatasetSuperGroup
from sfaira.data.dataloaders.databases.cellxgene.cellxgene_loader import Dataset
from sfaira.data.dataloaders.databases.cellxgene.rest_helpers import get_collection, get_collections


class DatasetGroupCellxgene(DatasetGroup):

    collection_id: str

    def __init__(
            self,
            collection_id: str = "default",
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            cache_metadata: bool = False,
            verbose: int = 0,
    ):
        self._collection = None
        self._cache_metadata = cache_metadata
        self._data_dir = data_path
        # _collection_id needs to be set here so that dataset_ids can be accessed from cache:
        self.collection_id = collection_id
        dataset_ids = [x["id"] for x in self.collection['datasets']]
        if len(dataset_ids) == 0:
            print(f"WARNING: Zero data sets retrieved for cellxgene collection {collection_id}.")
        loader_pydoc_path_sfaira = "sfaira.data.dataloaders.databases.cellxgene.cellxgene_loader"
        load_func = pydoc.locate(loader_pydoc_path_sfaira + ".load")
        datasets = [
            Dataset(
                collection_id=collection_id,
                data_path=data_path,
                meta_path=meta_path,
                cache_path=cache_path,
                load_func=load_func,
                sample_fn=x,
                sample_fns=dataset_ids,
                cache_metadata=cache_metadata,
                verbose=verbose,
            )
            for x in dataset_ids
        ]
        keys = [x.id for x in datasets]
        super().__init__(datasets=dict(zip(keys, datasets)), collection_id=collection_id)

    @property
    def _collection_cache_dir(self):
        """
        Collection meta data caching directory.

        Caches are writtgen into specific directory under HOME/.cache with other sfaira meta data unless
        path_data was supplied to constructor. This co-localisation of collection meta data with data downloads eases
        transfer of cache files between machines.
        """
        if self._data_dir is None:
            return settings.cachedir_databases_cellxgene
        else:
            return os.path.join(self._data_dir, "meta")

    @property
    def _collection_cache_fn(self):
        return os.path.join(self._collection_cache_dir, self.collection_id + ".pickle")

    @property
    def collection(self):
        if self._collection is None:
            self._collection = get_collection(collection_id=self.collection_id, cache=self._collection_cache_fn)
        return self._collection

    def show_summary(self):
        uuid_session = str(uuid.uuid4())
        display_html('<div id="{}" style="height: 600px; width:100%;"></div>'.format(uuid_session), raw=True)
        display_javascript("""
        require(["https://rawgit.com/caldwell/renderjson/master/renderjson.js"], function() {
          document.getElementById('%s').appendChild(renderjson(%s))
        });
        """ % (uuid_session, json.dumps(self.collection)), raw=True)


class DatasetSuperGroupCellxgene(DatasetSuperGroup):

    def __init__(
            self,
            data_path: Union[str, None] = None,
            meta_path: Union[str, None] = None,
            cache_path: Union[str, None] = None,
            cache_metadata: bool = False,
            verbose: int = 0,
    ):
        self._collections = None
        self._cache_metadata = cache_metadata
        self._data_dir = data_path
        # Get all collection IDs and instantiate one data set group per collection.
        # Throw warning in collection are not obtained:
        collections = self.collections
        if len(collections) == 0:
            print("WARNING: Zero cellxgene collections retrieved.")
        # Note that the collection itself is not passed to DatasetGroupCellxgene but only the ID string.
        dataset_groups = [
            DatasetGroupCellxgene(
                collection_id=x["id"],
                data_path=data_path,
                meta_path=meta_path,
                cache_path=cache_path,
                cache_metadata=cache_metadata,
                verbose=verbose,
            )
            for x in collections
        ]
        super().__init__(dataset_groups=dataset_groups)

    @property
    def _collections_cache_dir(self):
        """
        Collections summary caching directory.

        Caches are written into specific directory under HOME/.cache with other sfaira meta data unless
        path_data was supplied to constructor. This co-localisation of collection meta data with data downloads eases
        transfer of cache files between machines.
        """
        if self._data_dir is None:
            return settings.cachedir_databases_cellxgene
        else:
            return os.path.join(self._data_dir, "meta")

    @property
    def _collections_cache_fn(self):
        return os.path.join(self._collections_cache_dir, "collections.pickle")

    @property
    def collections(self):
        """
        Yield all collections available from REST API.
        """
        if self._collections is None:
            self._collections = get_collections(cache=self._collections_cache_fn)
        return self._collections

    def show_summary(self):
        """
        Prints overview of all collections available.
        """
        display("There are " + str(len(self.collections)) + " public collections sorting by newest first:")
        for collection in sorted(self.collections, key=lambda key: key['created_at'], reverse=True):
            display("id: " + collection['id'] + ' created on: ' +
                    datetime.date.fromtimestamp(collection['created_at']).strftime("%m/%d/%y"))
