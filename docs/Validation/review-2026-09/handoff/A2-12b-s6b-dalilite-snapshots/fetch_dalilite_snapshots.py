#!/usr/bin/env python3
"""Probes Internet Archive SNAPSHOTS directly (no CDX index: it was unreliable in run 1 of this question) for older DaliLite releases, keeps what is a real gzip / zip / compress / bzip2 archive (magic bytes), reads CHANGES / README / INSTALL / Makefile / the
`DaliLite` script and the Perl scripts as TEXT ONLY (python tarfile / zipfile in memory: nothing is extracted to disk, compiled or run) and reports whether any release prints `Query` / `Sbjct` to standard output in its `-pairwise` branch.
Every request is logged (one line: id, url, rc / HTTP status, final URL, size, magic); a failing request never stops the run. Neutral User-Agent `Seqeron-review/1.0`; no account data.
Usage: fetch_dalilite_snapshots.py OUTDIR"""
import concurrent.futures as cf, hashlib, io, json, os, re, subprocess, sys, tarfile, threading, time, urllib.parse, zipfile

UA = 'Seqeron-review/1.0'
out = sys.argv[1]
for d in ('found', 'inspect', 'listings', 'avail'):
    os.makedirs(f'{out}/{d}', exist_ok=True)
lock = threading.Lock()
logf = open(f'{out}/requests.log', 'w')
counter = [0]


def log(msg):
    with lock:
        logf.write(msg + '\n')
        logf.flush()
        print(msg, flush=True)


def nextid():
    with lock:
        counter[0] += 1
        return counter[0]


def curl(url, dest=None, head=False, max_time=240, rng=None):
    """One request through curl (follows redirects). Returns (rc, http_status, final_url, content_type, size). Never raises."""
    rid = nextid()
    cmd = ['curl', '-sS', '-L', '-A', UA, '--connect-timeout', '25', '--max-time', str(max_time), '--max-redirs', '8', '-w', '\n@@%{http_code}|%{url_effective}|%{content_type}|%{size_download}', '-o', dest or '/dev/null']
    if head:
        cmd.insert(1, '-I')
    if rng:
        cmd[1:1] = ['-r', rng]
    cmd.append(url)
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=max_time + 30)
        tail = (r.stdout or '').split('@@')[-1] if '@@' in (r.stdout or '') else '||'
        parts = (tail.split('|') + ['', '', '', ''])[:4]
        res = (r.returncode, parts[0], parts[1], parts[2], int(parts[3]) if parts[3].strip().isdigit() else 0)
        log(f'REQ {rid}\t{"HEAD" if head else "GET"}\t{url}\trc={res[0]}\thttp={res[1]}\tfinal={res[2]}\ttype={res[3]}\tbytes={res[4]}\t{(r.stderr or "").strip()[:120]}')
        return res
    except Exception as e:  # noqa: BLE001
        log(f'REQ {rid}\t{"HEAD" if head else "GET"}\t{url}\tEXCEPTION {type(e).__name__} {e}')
        return (-1, '000', '', '', 0)


def magic(path):
    try:
        b = open(path, 'rb').read(512)
    except OSError:
        return 'none'
    if b[:2] == b'\x1f\x8b': return 'gzip'
    if b[:4] == b'PK\x03\x04': return 'zip'
    if b[:2] == b'\x1f\x9d': return 'compress'
    if b[:3] == b'BZh': return 'bzip2'
    if b[257:262] == b'ustar': return 'tar'
    if re.match(rb'\s*<', b) or b'<html' in b.lower(): return 'html'
    return 'other'


# 1. candidate release names x base URLs
VERSIONS = ['1.0', '1.1', '1.2', '2.0', '2.1', '2.2', '2.3', '2.4', '2.4.1', '2.4.2', '2.4.3', '2.4.4', '2.4.5', '3.0', '3.1', '3.2', '3.3']
NAMES = []
for v in VERSIONS:
    NAMES += [f'DaliLite_{v}.tar.gz', f'DaliLite.v{v}.tar.gz', f'DaliLite-{v}.tar.gz', f'DaliLite_{v}.tgz']
NAMES += ['DaliLite.tar.gz', 'DaliLite.zip', 'DaliLite.tgz', 'DaliLite.tar.Z', 'dalilite.tar.gz', 'DaliLite_old.tar.gz', 'DaliLite.v2.tar.gz', 'DaliLite.v3.tar.gz', 'DaliLite.v4.tar.gz']
BASES = ['http://www.ebi.ac.uk/dali/DaliLite/', 'http://www.ebi.ac.uk/DaliLite/', 'http://www.ebi.ac.uk/dali/', 'http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/', 'http://ekhidna.biocenter.helsinki.fi/dali_lite/',
         'http://ekhidna.biocenter.helsinki.fi/dali/', 'http://ekhidna.biocenter.helsinki.fi/dali/DaliLite/', 'http://ekhidna.biocenter.helsinki.fi/DaliLite/', 'http://ekhidna2.biocenter.helsinki.fi/dali/',
         'http://ekhidna2.biocenter.helsinki.fi/dali/downloads/', 'http://ekhidna2.biocenter.helsinki.fi/dali_lite/downloads/', 'http://www.bioinformatics.org/dali/', 'http://www.bioinformatics.org/project/',
         'ftp://ftp.ebi.ac.uk/pub/databases/dali/', 'http://ftp.ebi.ac.uk/pub/databases/dali/', 'ftp://ftp.ebi.ac.uk/pub/software/dali/', 'http://ftp.ebi.ac.uk/pub/software/dali/', 'http://www.embl-ebi.ac.uk/DaliLite/', 'http://www.ebi.ac.uk/~holm/']
YEARS = list(range(2003, 2014))
if os.environ.get('DALI_QUICK'):                       # self-test only: 2 names x 2 bases
    NAMES, BASES = NAMES[:2], BASES[:2]
log(f'PLAN names={len(NAMES)} bases={len(BASES)} years={len(YEARS)} -> {len(NAMES) * len(BASES)} URLs, 3 stage-A probes each, 11 more for each URL that answers')

kept = {}       # sha256 -> (name, url, magic, size)
tried = set()


def probe(url, year):
    """web.archive.org/web/<year>id_/<url> redirects to the capture nearest to that year. A HEAD request first (cheap); the body is fetched only when the answer is 2xx and not HTML; a real archive has the magic bytes of an archive.
    Returns the HTTP status of the HEAD (so the caller can tell 'something is archived under this URL' from 'nothing')."""
    key = (url, year)
    if key in tried:
        return None
    tried.add(key)
    snap = f'https://web.archive.org/web/{year}id_/{url}'
    rc, http, final, ctype, size = curl(snap, head=True, max_time=90)
    if not http.startswith('2') or 'html' in ctype.lower():
        return http
    tmp = f'{out}/found/.tmp_{threading.get_ident()}'
    rc, http, final, ctype, size = curl(snap, dest=tmp, max_time=900)
    m = magic(tmp)
    if http.startswith('2') and size > 0 and m in ('gzip', 'zip', 'compress', 'bzip2', 'tar'):
        sha = hashlib.sha256(open(tmp, 'rb').read()).hexdigest()
        with lock:
            fresh = sha not in kept
            if fresh:
                n = len(kept) + 1
                name = f'{n:02d}_{re.sub(r"[^A-Za-z0-9._-]", "_", os.path.basename(urllib.parse.urlparse(url).path))}'
                os.replace(tmp, f'{out}/found/{name}')
                kept[sha] = (name, snap, final, m, size)
        log(f'KEPT\t{"new" if fresh else "duplicate"}\t{snap}\tfinal={final}\t{m}\t{size}\tsha256={sha}')
    elif os.path.exists(tmp):
        os.remove(tmp)
    return http


# 2. the availability API for every (base, name) (what archive.org says it holds) - both the documented endpoint and the one named in the task
def avail(url):
    for api in ('https://archive.org/wayback/available?url=', 'https://archive.org/wait/available?url='):
        dest = f'{out}/avail/{hashlib.sha1((api + url).encode()).hexdigest()[:12]}.json'
        rc, http, final, ctype, size = curl(api + urllib.parse.quote(url, safe=''), dest=dest, max_time=60)
        try:
            d = json.load(open(dest))
            snap = d.get('archived_snapshots', {}).get('closest', {})
            if snap.get('available'):
                log(f'AVAILABLE\t{url}\t{snap.get("timestamp")}\t{snap.get("url")}')
        except Exception:  # noqa: BLE001
            pass


urls = [b + n for b in BASES for n in NAMES]
with cf.ThreadPoolExecutor(8) as ex:
    list(ex.map(avail, urls))

# 3. directory listings of every base for every year: the real names of the releases (read as text, nothing else)
found_names = set()


def listing(args):
    base, year = args
    dest = f'{out}/listings/{hashlib.sha1(f"{base}{year}".encode()).hexdigest()[:10]}_{year}.html'
    rc, http, final, ctype, size = curl(f'https://web.archive.org/web/{year}id_/{base}', dest=dest, max_time=120)
    try:
        t = open(dest, 'rb').read(2_000_000).decode('latin-1')
    except OSError:
        return
    for h in re.findall(r'href="([^"]+)"', t, re.I):
        if re.search(r'dali', h, re.I) and re.search(r'\.(tar\.gz|tgz|tar|zip|tar\.Z|Z|gz)$', h, re.I):
            found_names.add((base, os.path.basename(urllib.parse.unquote(h.split('?')[0]))))


with cf.ThreadPoolExecutor(8) as ex:
    list(ex.map(listing, [(b, y) for b in BASES for y in YEARS]))
log(f'LISTINGS found {len(found_names)} archive names: {sorted(found_names)}')

# 4. the snapshot probes. Stage A: three years (2003, 2008, 2013) for every URL - the archive answers with the capture nearest to the year, so one hit anywhere in time shows up. Stage B: ALL years 2003-2013 for every URL for which
#    stage A got anything but a plain 404 / 403 / connection failure (an archived redirect, an archived page or file), and for every name the listings showed.
extra = [b + n for b, n in sorted(found_names)]
allurls = list(dict.fromkeys(urls + extra))
stageA = {}


def a(u):
    stageA[u] = [probe(u, y) for y in (2003, 2008, 2013)]


with cf.ThreadPoolExecutor(8) as ex:
    list(ex.map(a, allurls))
interesting = [u for u in allurls if u in extra or any(h not in (None, '404', '403', '000', '400') for h in stageA.get(u, []))]
log(f'STAGE A {len(allurls)} URLs, {len(interesting)} interesting (an answer other than 404 / 403 / 000)')
with cf.ThreadPoolExecutor(8) as ex:
    list(ex.map(lambda t: probe(*t), [(u, y) for u in interesting for y in YEARS]))
log(f'KEPT {len(kept)} distinct archives')


# 5. inspect as TEXT (in memory)
def members_of(path, m):
    data = open(path, 'rb').read()
    if m == 'compress':
        data = subprocess.run(['gzip', '-dc'], input=data, capture_output=True).stdout
        m = 'tar'
    if m == 'bzip2':
        import bz2
        data = bz2.decompress(data)
        m = 'tar'
    names, texts = [], {}
    if m in ('gzip', 'tar'):
        with tarfile.open(fileobj=io.BytesIO(data), mode='r:*') as t:
            for mem in t:
                names.append(mem.name)
                if mem.isfile() and mem.size < 3_000_000 and re.search(r'(^|/)(CHANGES|README[^/]*|INSTALL|Makefile|DaliLite|LICEN[CS]E|[^/]*\.pl)$', mem.name):
                    texts[mem.name] = t.extractfile(mem).read().decode('latin-1')
    elif m == 'zip':
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for mem in z.namelist():
                names.append(mem)
                if re.search(r'(^|/)(CHANGES|README[^/]*|INSTALL|Makefile|DaliLite|LICEN[CS]E|[^/]*\.pl)$', mem) and z.getinfo(mem).file_size < 3_000_000:
                    texts[mem] = z.read(mem).decode('latin-1')
    return names, texts


report = open(f'{out}/report.txt', 'w')
for sha, (name, snap, final, m, size) in sorted(kept.items(), key=lambda kv: kv[1][0]):
    line = f'{name}\t{snap}\tfinal={final}\t{m}\tbytes={size}\tsha256={sha}'
    try:
        names, texts = members_of(f'{out}/found/{name}', m)
        d = f'{out}/inspect/{name}'
        os.makedirs(d, exist_ok=True)
        open(d + '/members.txt', 'w').write('\n'.join(names))
        for k, v in texts.items():
            open(d + '/' + re.sub(r'[^A-Za-z0-9._-]', '_', k), 'w').write(v)
        script = [v for k, v in texts.items() if re.search(r'(^|/)DaliLite$', k)]
        chg = '\n'.join(v for k, v in texts.items() if k.endswith('CHANGES'))
        ver = re.findall(r'(?:Changes in )?DaliLite v\.?\s*([0-9.]+)', chg)
        line += f'\tmembers={len(names)}\tnewest_CHANGES={ver[0] if ver else "?"}'
        line += '\tscript_has_Query=%s\tscript_has_Sbjct=%s' % tuple(str(any(s in sc for sc in script)) for s in ('Query', 'Sbjct'))
        # the lines of the script's -pairwise branch that print (print / echo / system "cat" / cat) - the reader must see them as text
        pr = [l.strip() for sc in script for l in sc.splitlines() if re.search(r'Query|Sbjct|print\s+STDOUT|print\s+"|echo\s+"|\bcat\b', l)][:40]
        open(d + '/pairwise_print_lines.txt', 'w').write('\n'.join(pr))
    except Exception as e:  # noqa: BLE001
        line += f'\tUNREADABLE {type(e).__name__} {e}'
    report.write(line + '\n')
    report.flush()
    log('INSPECTED ' + line)
report.close()
log(f'DONE requests={counter[0]} kept={len(kept)}')
