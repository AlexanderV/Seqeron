#!/usr/bin/env python3
"""Looks for OLDER DaliLite releases (v1 / v2.x, before v3.0 of Feb 2008) and inspects what it finds. Everything is logged to stdout; never raises; downloads nothing but candidate archives and never executes them.
Layout written: <out>/found/<n>_<name> (the archives), <out>/inspect/<n>/ (text files read out of an archive: CHANGES, README, INSTALL, Makefile, the `DaliLite` script) and <out>/report.txt (one line per archive: sha256, file list size, whether the `DaliLite` script
contains the strings Query / Sbjct, and the lines of the script that print to standard output in the -pairwise branch)."""
import hashlib, io, os, re, subprocess, sys, tarfile, urllib.parse, urllib.request, zipfile

UA = {'User-Agent': 'Seqeron-review/1.0'}
out = sys.argv[1]
os.makedirs(out + '/found', exist_ok=True)
os.makedirs(out + '/inspect', exist_ok=True)
report = open(out + '/report.txt', 'w')


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)


def get(url, binary=False, timeout=120):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = r.read()
    return data if binary else data.decode('latin-1')


def curl_resume(url, dest, attempts=40):
    for k in range(attempts):
        r = subprocess.run(['curl', '-sS', '-L', '-A', 'Seqeron-review/1.0', '-C', '-', '--connect-timeout', '30', '--max-time', '300', '-o', dest, url], capture_output=True, text=True)
        size = os.path.getsize(dest) if os.path.exists(dest) else 0
        log('CURL attempt', k + 1, 'rc', r.returncode, 'bytes', size, (r.stderr or '').strip()[:160])
        if r.returncode == 0 and size > 0:
            return True
    return False


# 1. the live servers: names that DaliLite releases had (HEAD only; a listing of what answers 200)
LIVE = ['http://ekhidna2.biocenter.helsinki.fi/dali/', 'http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/', 'http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/v2/', 'http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/v1/',
        'http://www.ebi.ac.uk/DaliLite/', 'http://www.ebi.ac.uk/dali/DaliLite/']
NAMES = ['DaliLite_%s.tar.gz' % v for v in ('1.0', '1', '2.0', '2.1', '2.2', '2.3', '2.4', '2.4.1', '2.4.2', '2.4.3', '2.4.4', '2.4.5', '3.0', '3.1', '3.2', '3.3')] + ['DaliLite.tar.gz', 'DaliLite.zip']
live_hits = []
for base in LIVE:
    for n in NAMES:
        url = base + n
        try:
            req = urllib.request.Request(url, headers=UA, method='HEAD')
            with urllib.request.urlopen(req, timeout=40) as r:
                log('LIVE 200', url, r.headers.get('Content-Length'))
                live_hits.append(url)
        except Exception as e:  # noqa: BLE001
            pass
log('LIVE hits', live_hits)

# 2. the Internet Archive: every captured URL (status 200) of the places DaliLite was distributed from whose name looks like an archive
cands = {}
PREFIXES = ['ekhidna.biocenter.helsinki.fi/dali_lite/', 'ekhidna.biocenter.helsinki.fi/dali/', 'ekhidna.biocenter.helsinki.fi/dali_server/', 'ekhidna.biocenter.helsinki.fi/downloads/', 'ekhidna.biocenter.helsinki.fi/DaliLite/',
            'ekhidna2.biocenter.helsinki.fi/dali/', 'www.ebi.ac.uk/DaliLite/', 'www.ebi.ac.uk/dali/', 'www.ebi.ac.uk/Tools/dalilite/', 'ftp.ebi.ac.uk/pub/databases/dali/', 'ftp.ebi.ac.uk/pub/software/', 'www.bioinformatics.org/project/',
            'www.bioinformatics.org/dali/', 'www.embl-ebi.ac.uk/DaliLite/', 'ftp.embl-ebi.ac.uk/pub/software/', 'www.ebi.ac.uk/~holm/', 'www.fi.embnet.org/', 'ekhidna.biocenter.helsinki.fi/pub/']
for pre in PREFIXES:
    url = 'http://web.archive.org/cdx/search/cdx?url=' + urllib.parse.quote(pre + '*', safe='*/.') + '&matchType=prefix&output=txt&fl=timestamp,original,statuscode,length&filter=statuscode:200&collapse=urlkey&limit=20000'
    try:
        txt = get(url)
    except Exception as e:  # noqa: BLE001
        log('CDX FAILED', pre, e)
        continue
    lines = [l for l in txt.splitlines() if l.strip()]
    log('CDX', pre, len(lines), 'captured URLs')
    for l in lines:
        f = l.split()
        if len(f) >= 2 and re.search(r'dali', f[1], re.I) and re.search(r'\.(tar\.gz|tgz|tar|zip|tar\.Z|Z|gz)$', f[1], re.I):
            if f[1] not in cands or f[0] > cands[f[1]]:
                cands[f[1]] = f[0]
log('CANDIDATES', len(cands))
for orig, ts in sorted(cands.items()):
    log('  ', ts, orig)

# 3. download every candidate (newest capture, raw mode) and inspect it as TEXT
n = 0
for orig, ts in sorted(cands.items()) + [(u, '') for u in live_hits]:
    n += 1
    name = '%02d_%s' % (n, re.sub(r'[^A-Za-z0-9._-]', '_', os.path.basename(orig)))
    dest = out + '/found/' + name
    url = ('http://web.archive.org/web/%sid_/%s' % (ts, orig)) if ts else orig
    if not curl_resume(url, dest, attempts=10):
        log('DOWNLOAD FAILED', url)
        continue
    data = open(dest, 'rb').read()
    sha = hashlib.sha256(data).hexdigest()
    line = '%s\t%s\t%s\tbytes=%d' % (name, url, sha, len(data))
    members = []
    texts = {}
    try:
        if data[:2] == b'\x1f\x8b' or data[:5] == b'ustar' or data[257:262] == b'ustar':
            with tarfile.open(fileobj=io.BytesIO(data), mode='r:*') as t:
                for m in t:
                    members.append(m.name)
                    if m.isfile() and m.size < 3_000_000 and re.search(r'(^|/)(CHANGES|README[^/]*|INSTALL|Makefile|DaliLite|LICENCE|[^/]*\.pl)$', m.name):
                        texts[m.name] = t.extractfile(m).read().decode('latin-1')
        elif data[:2] == b'PK':
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                for m in z.namelist():
                    members.append(m)
                    if re.search(r'(^|/)(CHANGES|README[^/]*|INSTALL|Makefile|DaliLite|LICENCE|[^/]*\.pl)$', m) and z.getinfo(m).file_size < 3_000_000:
                        texts[m] = z.read(m).decode('latin-1')
        else:
            line += '\tNOT AN ARCHIVE (first bytes %r)' % data[:12]
    except Exception as e:  # noqa: BLE001
        line += '\tUNREADABLE %s' % e
    line += '\tmembers=%d' % len(members)
    d = '%s/inspect/%02d' % (out, n)
    os.makedirs(d, exist_ok=True)
    open(d + '/members.txt', 'w').write('\n'.join(members))
    for k, v in texts.items():
        open(d + '/' + re.sub(r'[^A-Za-z0-9._-]', '_', k), 'w').write(v)
    script = [v for k, v in texts.items() if re.search(r'(^|/)DaliLite$', k)]
    ver = re.search(r'Changes in DaliLite v\.?\s*([0-9.]+)', '\n'.join(v for k, v in texts.items() if k.endswith('CHANGES')))
    line += '\tnewest_CHANGES=%s' % (ver.group(1) if ver else '?')
    line += '\tscript_has_Query=%s\tscript_has_Sbjct=%s' % tuple(str(any(s in sc for sc in script)) for s in ('Query', 'Sbjct'))
    report.write(line + '\n')
    report.flush()
    log('INSPECTED', line)
log('DONE', n, 'candidates')
report.close()
