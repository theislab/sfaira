import cgi
import os
import ssl
from tqdm import tqdm
import urllib.request
import urllib.parse
import urllib.error
import warnings


class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download(urls, data_dir, directory_formatted_doi, dataset_id, **kwargs):
    for url in urls:
        if url is None:
            continue
        # Special case for data that is not publicly available
        if url.split(",")[0] == 'private':
            if "," in url:
                fn = ','.join(url.split(',')[1:])
                if os.path.isfile(os.path.join(data_dir, fn)):
                    print(f"File {fn} already found on disk, skipping download.")
                else:
                    warnings.warn(f"Dataset {dataset_id} is not available for automatic download, please manually "
                                  f"copy the file {fn} to the following location: {data_dir}")
            else:
                warnings.warn(f"A file for dataset {dataset_id} is not available for automatic download, please"
                              f"manually copy the associated file to the following location: {data_dir}")
        # Special case for data from the synapse portal
        elif url.split(",")[0].startswith('syn'):
            fn = ",".join(url.split(",")[1:])
            if os.path.isfile(os.path.join(data_dir, fn)):
                print(f"File {fn} already found on disk, skipping download.")
            else:
                _download_synapse(synapse_entity=url.split(",")[0], data_dir=data_dir, fn=fn, dataset_id=dataset_id,
                                  **kwargs)
        # Special case for data from the braod single cell portal
        elif url.split(",")[0].startswith('SCP'):
            fn = url
            if os.path.isdir(os.path.join(data_dir, url)):
                print(f"SCP directory {fn} already found on disk, assuming complete content, skipping download.")
            else:
                warnings.warn(
                    f"Dataset {dataset_id} is not available for automatic download as it is served on the Broad Institute Single Cell Portal. "
                    f"To retreive the data, follow the steps below.\n"
                    f"Step 1: Log in to the Broad Single Cell Portal at https://singlecell.broadinstitute.org\n"
                    f"Step 2: Create a bulk download curl command using the 'Bulk download' button on the following site: "
                    f"https://singlecell.broadinstitute.org/single_cell/study/{fn}#study-download\n"
                    f"Step 3: Change directory in your terminal to {data_dir} and execute the curl command obtained in Step 2 there."
                )
        # Special case for public data that is labelled as not automatically downloadable
        elif url.split(",")[0] == 'manual':
            u = ",".join(url.split(",")[2:])
            fn = url.split(",")[1]
            if os.path.isfile(os.path.join(data_dir, fn)):
                print(f"File {fn} already found on disk, skipping download.")
            else:
                print(f"Data file {fn} for dataset {dataset_id} cannot be retrieved automatically. "
                      f"Please download it from {u} and copy to {os.path.join(data_dir, fn)}")
        # All other cases
        else:
            if url.split(",")[0] == 'rename':
                rename = url.split(",")[1]
                url = ",".join(url.split(",")[2:])
            else:
                rename = None
            url = urllib.parse.unquote(url)
            try:
                urllib.request.urlopen(url)
            except urllib.error.HTTPError as err:
                # modify headers if urllib useragent is blocked (eg.10x datasets)
                if err.code == 403:
                    opener = urllib.request.build_opener()
                    opener.addheaders = [('User-Agent', 'Mozilla/5.0 (Windows NT 6.1; WOW64)')]
                    urllib.request.install_opener(opener)
            except urllib.error.URLError:
                # Catch SSLCertVerificationError: [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable
                # to get local issuer certificate (_ssl.c:1124)
                ssl._create_default_https_context = ssl._create_unverified_context

            if rename is not None:
                fn = rename
            elif 'Content-Disposition' in urllib.request.urlopen(url).info().keys():
                fn = cgi.parse_header(urllib.request.urlopen(url).info()['Content-Disposition'])[1]["filename"]
            else:
                fn = url.split("/")[-1]
            # Only download if file not already downloaded:
            if os.path.isfile(os.path.join(data_dir, fn)):
                print(f"File {fn} already found on disk, skipping download.")
            else:
                print_dir = f"{directory_formatted_doi}:{fn}"
                with DownloadProgressBar(unit='B', unit_scale=True, miniters=1, desc=print_dir) as pbar:
                    urllib.request.urlretrieve(url, os.path.join(data_dir, fn), reporthook=pbar.update_to)


def _download_synapse(synapse_entity, data_dir, fn, dataset_id, **kwargs):
    try:
        import synapseclient
    except ImportError:
        warnings.warn("synapseclient python package not found. This package is required to download some of the "
                      "selected datasets. Run `pip install synapseclient` to install it. Skipping download of the "
                      f"following dataset: {dataset_id}")
        return
    import logging
    logging.captureWarnings(False)  # required to properly display warning messages below with sypaseclient loaded

    if "synapse_user" not in kwargs.keys():
        warnings.warn(f"No synapse username provided, skipping download of synapse dataset {fn}."
                      f"Provide your synapse username as the `synapse_user` argument to the download method.")
        return
    if "synapse_pw" not in kwargs.keys():
        warnings.warn(f"No synapse password provided, skipping download of synapse dataset {fn}."
                      f"Provide your synapse password as the `synapse_pw` argument to the download method.")
        return

    print(f"Downloading from synapse: {fn}")
    syn = synapseclient.Synapse()
    syn.login(kwargs['synapse_user'], kwargs['synapse_pw'])
    syn.get(entity=synapse_entity, downloadLocation=os.path.join(data_dir, fn))
    syn.logout()
