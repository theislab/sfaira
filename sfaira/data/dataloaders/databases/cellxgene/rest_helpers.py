"""
Helper functionalities to interact with cellxgene REST API.
"""

import os
import pathlib
import pickle
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

CELLXGENE_PRODUCTION_ENDPOINT = 'https://api.cellxgene.cziscience.com'
COLLECTIONS = "/dp/v1/collections/"
DOWNLOAD_DATASET = "/dp/v1/datasets/"
MAX_RETRIES = 5
TIMEOUT_COLLECTION = 5
TIMEOUT_DATA = 10
HTTP_ERROR_LIST = [429, 502, 504]


class CustomHTTPAdapter(HTTPAdapter):
    def __init__(self, timeout, **kwargs):
        self.timeout = timeout
        super().__init__(**kwargs)

    def send(self, request, **kwargs):
        kwargs["timeout"] = self.timeout
        return super().send(request, **kwargs)


def rest_api_collection_request(url):
    retry_strategy = Retry(
        backoff_factor=0,
        method_whitelist=["GET"],
        status_forcelist=HTTP_ERROR_LIST,
        total=MAX_RETRIES,
    )
    adapter = CustomHTTPAdapter(timeout=TIMEOUT_COLLECTION, max_retries=retry_strategy)
    https = requests.Session()
    https.mount("https://", adapter)
    r = https.get(url)
    r.raise_for_status()
    return r.json()


def rest_api_data_request(presigned_url):
    retry_strategy = Retry(
        backoff_factor=0,
        method_whitelist=["GET"],
        status_forcelist=HTTP_ERROR_LIST,
        total=MAX_RETRIES,
    )
    adapter = CustomHTTPAdapter(timeout=TIMEOUT_DATA, max_retries=retry_strategy)
    https = requests.Session()
    https.mount("https://", adapter)
    r = https.get(presigned_url)
    return r.content


def get_collection(collection_id, cache=None):
    """
    Meta data of a specific collection in cellxgene.

    :param collection_id: Collection to query meta data for.
    :param cache: File name of pickle to cache this summary to.
    :return: Dictionary over meta data fields.
    """
    # Check if cached:
    if cache is not None and os.path.exists(cache):
        with open(cache, "rb") as f:
            collection = pickle.load(f)
    else:
        # Download and cache:
        collection = rest_api_collection_request(url=CELLXGENE_PRODUCTION_ENDPOINT + COLLECTIONS + collection_id)
        if cache:
            parent_dir = pathlib.Path(cache).parents[0]
            pathlib.Path(parent_dir).mkdir(parents=True, exist_ok=True)
            with open(cache, "wb") as f:
                pickle.dump(obj=collection, file=f)
    return collection


def get_collections(cache=None):
    """
    Summary of all collections in cellxgene.

    :param cache: File name of pickle to cache this summary to.
    :return: Summary as list over collections, each as a dictionary over meta data fields.
    """
    # Check if cached:
    if cache is not None and os.path.exists(cache):
        with open(cache, "rb") as f:
            collections = pickle.load(f)
    else:
        # Download and cache:
        collections = rest_api_collection_request(url=CELLXGENE_PRODUCTION_ENDPOINT + COLLECTIONS)['collections']
        if cache:
            parent_dir = pathlib.Path(cache).parents[0]
            pathlib.Path(parent_dir).mkdir(parents=True, exist_ok=True)
            with open(cache, "wb") as f:
                pickle.dump(obj=collections, file=f)
    return collections


def get_data(presigned_url):
    return rest_api_data_request(presigned_url=presigned_url)
