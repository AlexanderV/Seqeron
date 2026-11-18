#!/usr/bin/env python3
"""Fetches src/ of Debian's clustalw 2.1+lgpl-7 (sources.debian.org JSON API + raw files) into the current directory. Neutral User-Agent."""
import json, os, urllib.parse, urllib.request

BASE = 'https://sources.debian.org'
VER = '2.1%2Blgpl-7'


def get(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'Seqeron-review/1.0'}), timeout=60).read()


def walk(path):
    for e in json.loads(get(f'{BASE}/api/src/clustalw/{VER}/{path}/'))['content']:
        p = f"{path}/{e['name']}"
        if e['type'] == 'directory':
            walk(p)
        else:
            os.makedirs(os.path.dirname(p), exist_ok=True)
            with open(p, 'wb') as f:
                f.write(get(f'{BASE}/data/main/c/clustalw/{VER}/{urllib.parse.quote(p)}'))


walk('src')
print('fetched')
