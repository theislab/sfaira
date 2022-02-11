import itertools
import urllib.request

import sfaira

u = sfaira.data.Universe()
urls = []
for x in u.datasets.values():
    a = x.download_url_data
    if a is not None:
        urls.append(a) if isinstance(a, str) else urls.extend(a)
    b = x.download_url_meta
    if b is not None:
        urls.append(b) if isinstance(b, str) else urls.extend(b)

flat_urls = list(filter(lambda url: url is not None, list(itertools.chain(*urls))))
for url in flat_urls:
    if not url.startswith(("syn", "manual", "private")):
        assert urllib.request.urlopen(url).getcode() == 200
