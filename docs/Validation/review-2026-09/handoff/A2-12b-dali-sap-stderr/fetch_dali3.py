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


def curl_resume(url, dest, attempts=40):
    """Run 3 of the package: the archive cuts the transfer (IncompleteRead after 13 032 of 1 354 361 bytes of DaliLite.tar.gz). Resume with `curl -C -` until the file is complete (curl exits 0)
    or the attempts are used up. Only the program archive is wanted (the DAT / GTG / dump archives are 0.7 - 6 GB databases of structures, not needed for `DaliLite -pairwise`)."""
    for k in range(attempts):
        r = subprocess.run(['curl', '-sS', '-L', '-A', 'Seqeron-review/1.0', '-C', '-', '--connect-timeout', '30', '--max-time', '300', '-o', dest, url], capture_output=True, text=True)
        size = os.path.getsize(dest) if os.path.exists(dest) else 0
        print('CURL attempt', k + 1, 'rc', r.returncode, 'bytes', size, (r.stderr or '').strip()[:160])
        if r.returncode == 0 and size > 0:
            return True
    return False


# Step 00: the program archive alone, by the URLs that run 3 showed to exist (2012id_ and the 2014 snapshot of the page), resumed with curl until complete.
for ts in ('2012', '20140722051239', '2016', '2011'):
    url = 'http://web.archive.org/web/%sid_/http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/v3/DaliLite.tar.gz' % ts
    dest = out + '/src/DaliLite.tar.gz'
    if os.path.exists(dest): os.remove(dest)
    if curl_resume(url, dest):
        with open(dest, 'rb') as fh: head = fh.read(8)
        if head[:2] == b'\x1f\x8b':
            print('DOWNLOADED (curl, resumed)', url, os.path.getsize(dest))
            got.append(dest)
            break
        print('NOT A GZIP', url, head)
        os.remove(dest)
# Step 0 (run 3 showed: CDX returned 0 lines; run 2 of the package showed that the links of the snapshot page resolve to captured-less URLs, "http:/" joined and the /v3/ part lost): ask the CDX index of the archive which
# tarball URLs of the DaliLite download area were captured with status 200 and fetch the newest capture of each in raw mode (timestamp + "id_"). Logged; never raises.
try:
    if got: raise RuntimeError('program archive already downloaded')
    cdx = get('http://web.archive.org/cdx/search/cdx?url=ekhidna.biocenter.helsinki.fi/dali_lite/downloads/*&matchType=prefix&output=txt&fl=timestamp,original,statuscode&filter=statuscode:200&limit=2000')
    print('CDX', len(cdx.splitlines()), 'lines')
    newest = {}
    for ln in cdx.splitlines():
        f = ln.split()
        if len(f) >= 2 and re.search(r'(tar\.gz|tgz|\.zip)$', f[1], re.I):
            if f[1] not in newest or f[0] > newest[f[1]]: newest[f[1]] = f[0]
    print('CDX tarballs', newest)
    for orig, ts in sorted(newest.items()):
        url = 'http://web.archive.org/web/%sid_/%s' % (ts, orig)
        name = '_'.join(orig.split('/')[-2:]) if '/v3/' in orig else os.path.basename(orig)
        try:
            data = get(url, True)
            if data[:5].lower() in (b'<!doc', b'<html'): raise ValueError('HTML instead of an archive')
            open(out + '/src/' + name, 'wb').write(data)
            print('DOWNLOADED', url, len(data))
            got.append(out + '/src/' + name)
        except Exception as e:  # noqa: BLE001
            print('DOWNLOAD FAILED', url, e)
except Exception as e:  # noqa: BLE001
    print('CDX FAILED', e)
for page in ([] if got else PAGES):
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
            url = re.sub(r'(web\.archive\.org/web/\d+)/http:/(?!/)', r'\1id_/http://', url)
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
