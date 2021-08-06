"""
Helper functionalities to interact with cellxgene REST API.
"""

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


def get_collection(collection_id):
    return rest_api_collection_request(url=CELLXGENE_PRODUCTION_ENDPOINT + COLLECTIONS + collection_id)


def get_collections():
    return rest_api_collection_request(url=CELLXGENE_PRODUCTION_ENDPOINT + COLLECTIONS)['collections']


def get_data(presigned_url):
    return rest_api_data_request(presigned_url=presigned_url)
