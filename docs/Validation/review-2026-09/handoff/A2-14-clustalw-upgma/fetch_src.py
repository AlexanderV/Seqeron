#!/usr/bin/env python3
"""Fetches the C++ sources of Debian's clustalw 2.1+lgpl-7 (src/**.cpp, src/**.h) from sources.debian.org into ./src (neutral User-Agent) and writes the two generated headers the
autotools build would produce (config.h empty, clustalw_version.h = "2.1"). Prints the sha256 of every file to fetch_src.sha256."""
import hashlib, json, os, urllib.parse, urllib.request

BASE = 'https://sources.debian.org'
VER = '2.1+lgpl-7'


def get(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'Seqeron-review/1.0'}), timeout=60).read()


def walk(path):
    for e in json.loads(get(f'{BASE}/api/src/clustalw/{VER}/{path}'))['content']:
        if e['type'] == 'directory':
            yield from walk(path + e['name'] + '/')
        else:
            yield path + e['name']


sums, n = [], 0
for f in walk('src/'):
    if f.endswith(('.cpp', '.h')):
        os.makedirs(os.path.dirname(f), exist_ok=True)
        data = get(f"{BASE}/data/main/c/clustalw/{VER}/{urllib.parse.quote(f)}")
        open(f, 'wb').write(data)
        sums.append(f'{hashlib.sha256(data).hexdigest()}  {f}')
        n += 1
open('src/config.h', 'w').close()
open('src/clustalw_version.h', 'w').write('#define CLUSTALW_VERSION "2.1"\n#define VERSION "2.1"\n')
open('fetch_src.sha256', 'w').write('\n'.join(sums) + '\n')
print('fetched', n, 'files')
