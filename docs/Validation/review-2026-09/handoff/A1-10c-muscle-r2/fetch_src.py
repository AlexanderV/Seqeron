#!/usr/bin/env python3
"""Fetches the source tree of Debian's muscle3 3.8.1551-2 (the final MUSCLE 3.8 release; flat directory of .cpp/.h files) from sources.debian.org into the
current directory (JSON API + raw files). Neutral User-Agent. Prints the sha256 of every file to fetch_src.sha256."""
import hashlib, json, urllib.parse, urllib.request

BASE = 'https://sources.debian.org'
VER = '3.8.1551-2'


def get(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'Seqeron-review/1.0'}), timeout=60).read()


listing = json.loads(get(f'{BASE}/api/src/muscle3/{VER}/'))['content']
sums = []
for e in listing:
    if e['type'] != 'file':
        continue
    data = get(f"{BASE}/data/main/m/muscle3/{VER}/{urllib.parse.quote(e['name'])}")
    with open(e['name'], 'wb') as f:
        f.write(data)
    sums.append(f"{hashlib.sha256(data).hexdigest()}  {e['name']}")
open('fetch_src.sha256', 'w').write('\n'.join(sums) + '\n')
print('fetched', len(sums), 'files')
