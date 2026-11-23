#!/usr/bin/env python3
"""Fetches DaliLite v3 from the Internet Archive snapshot of its download page and tries to build it. Everything is logged to stdout; never raises.
Layout written: <out>/src/<archive>, <out>/DaliLite.v3 (unpacked), <out>/DaliLite.wrap (a bash wrapper around the built `DaliLite` that logs every call to /tmp/dali_log), <out>/BUILD_OK or <out>/BUILD_FAILED."""
import os, re, subprocess, sys, urllib.parse, urllib.request

UA = {'User-Agent': 'Seqeron-review/1.0'}
out = sys.argv[1]
os.makedirs(out + '/src', exist_ok=True)
PAGES = ['http://web.archive.org/web/20140722051239/http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/v3',
         'http://web.archive.org/web/2016/http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/v3',
         'http://web.archive.org/web/2012/http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/v3/',
         'http://web.archive.org/web/2011/http://ekhidna.biocenter.helsinki.fi/dali_lite/']


def get(url, binary=False):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read()
    return data if binary else data.decode('latin-1')


got = []
for page in PAGES:
    try:
        html = get(page)
        print('PAGE', page, len(html))
    except Exception as e:  # noqa: BLE001
        print('PAGE FAILED', page, e)
        continue
    links = sorted(set(re.findall(r'href="([^"]+)"', html)))
    print('LINKS', links)
    for l in links:
        if re.search(r'\.(tar\.gz|tgz|zip|tar)$', l, re.I) or re.search(r'DaliLite', l):
            url = urllib.parse.urljoin(page, l) if hasattr(urllib, 'parse') else l
            if not url.startswith('http'): url = 'http://web.archive.org' + l
            name = os.path.basename(url.split('?')[0]) or 'download'
            try:
                data = get(url, True)
                open(out + '/src/' + name, 'wb').write(data)
                print('DOWNLOADED', url, len(data))
                got.append(out + '/src/' + name)
            except Exception as e:  # noqa: BLE001
                print('DOWNLOAD FAILED', url, e)
    if got: break
for f in got:
    sh = lambda c: subprocess.run(c, shell=True, capture_output=True, text=True, cwd=out)  # noqa: E731
    print(sh(f'sha256sum {f}; file {f}').stdout)
    if re.search(r'\.zip$', f): print(sh(f'unzip -o -q {f} -d {out}/DaliLite.v3').stderr)
    else: print(sh(f'mkdir -p {out}/DaliLite.v3 && tar xf {f} -C {out}/DaliLite.v3').stderr)
    print(sh(f'find {out}/DaliLite.v3 -maxdepth 3 | head -60').stdout)
    r = sh(f'cd $(dirname $(find {out}/DaliLite.v3 -name Makefile -o -name makefile | head -1)) && make 2>&1 | tail -30')
    print(r.stdout, r.stderr)
    exe = sh(f'find {out}/DaliLite.v3 -type f -name DaliLite -perm -u+x | head -1').stdout.strip()
    if exe:
        open(out + '/DaliLite.wrap', 'w').write('#!/bin/bash\nmkdir -p /tmp/dali_log\nL=/tmp/dali_log/call_$$_$(date +%s%N).txt\n{ echo "ARGS: $*"; echo "CWD: $PWD"; } > $L\n'
                                               f'{exe} "$@" > /tmp/dali_out_$$.txt 2>/tmp/dali_err_$$.txt; rc=$?\n'
                                               'echo "RC: $rc" >> $L; echo "--- STDOUT" >> $L; cat /tmp/dali_out_$$.txt >> $L; echo "--- STDERR" >> $L; cat /tmp/dali_err_$$.txt >> $L\n'
                                               'cat /tmp/dali_out_$$.txt; cat /tmp/dali_err_$$.txt >&2; rm -f /tmp/dali_out_$$.txt /tmp/dali_err_$$.txt; exit $rc\n')
        os.chmod(out + '/DaliLite.wrap', 0o755)
        open(out + '/BUILD_OK', 'w').write(exe + '\n')
        print('BUILD_OK', exe)
        break
else:
    open(out + '/BUILD_FAILED', 'w').write('no DaliLite v3 executable\n')
    print('BUILD_FAILED')
