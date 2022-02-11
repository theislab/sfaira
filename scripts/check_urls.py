import itertools
import urllib
from typing import List
from urllib.error import HTTPError, URLError

import sfaira

u = sfaira.data.Universe()
urls: List[str] = []
for x in u.datasets.values():
    a = x.download_url_data
    if a is not None:
        urls.append(a) if isinstance(a, str) else urls.extend(a)
    b = x.download_url_meta
    if b is not None:
        urls.append(b) if isinstance(b, str) else urls.extend(b)

flat_urls = list(filter(lambda url: url is not None, list(itertools.chain(*urls))))
failed_urls: List[List[str]] = []
for url in flat_urls:
    try:
        if not url.startswith(("syn", "manual", "private")):
            assert urllib.request.urlopen(url).getcode() == 200
    except (AssertionError, ValueError, HTTPError, URLError) as e:
        failed_urls.append([url, e])

if len(failed_urls) != 0:
    for fail in failed_urls:
        print(fail)
    raise AssertionError(f"{len(failed_urls)} URLs failed.")
