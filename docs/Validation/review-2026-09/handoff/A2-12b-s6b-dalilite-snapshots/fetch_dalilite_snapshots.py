#!/usr/bin/env python3
"""RUN 2 (B09 WP81, F371). Probes Internet Archive SNAPSHOTS for older DaliLite releases, keeps what is a real gzip / zip / compress / bzip2 archive (magic bytes), reads CHANGES / README / INSTALL / Makefile / the `DaliLite` script and the Perl
scripts as TEXT ONLY (python tarfile / zipfile in memory: nothing is extracted to disk, compiled or run) and reports whether any release prints `Query` / `Sbjct` to standard output in its `-pairwise` branch.

Why run 2: run 1 was DEFECTIVE. It sent 8 parallel requests with no delay and no retry: archive.org answered the availability API with HTTP 429 (1 363 of 1 463 requests), then refused connections altogether (curl rc=7 for 4 302 of 4 389 snapshot probes
and 189 of 209 directory listings, every one after request 2 947), and the script counted "no connection" as "not archived" ("0 interesting") and wrote an empty report.txt. Only 20 listings, 100 availability answers and 87 probes (all 404) are real answers of run 1.
Run 2 asks politely: ONE request at a time, at least DALI_DELAY seconds (default 2.5) apart, retry with exponential back-off on 429 / 5xx / connection failures (Retry-After honoured), a pause when 10 requests in a row fail, and a three-way result for
every request (OK / NOT-FOUND (HTTP 404) / FAILED): a FAILED request is never read as "not archived". report.txt always starts with a SUMMARY line.
Order of work: (A) leads: the pages that run 1 showed to be captured (`www.ebi.ac.uk/~holm/DaliLite` = the 2002 redirect target of `www.ebi.ac.uk/dali/DaliLite/`, `www.ebi.ac.uk/DaliLite/`) for the years 2002 - 2013, their links followed one level;
(B) the CDX index of every base URL (one request each; known to be unreliable, cheap); (C) one HEAD probe per (base URL, release name) at web/2008id_/ (the archive answers with the capture nearest in time of the exact URL); (D) all years for every URL
that answered with anything but 404. Every request is logged in requests.log (id, URL, curl rc, HTTP status, final URL, content type, bytes, attempt). Neutral User-Agent `Seqeron-review/1.0`; no account data.
Usage: fetch_dalilite_snapshots.py OUTDIR"""
import hashlib, io, json, os, re, subprocess, sys, tarfile, time, urllib.parse, zipfile

UA = 'Seqeron-review/1.0'
out = sys.argv[1]
for d in ('found', 'inspect', 'listings', 'cdx'):
    os.makedirs(f'{out}/{d}', exist_ok=True)
logf = open(f'{out}/requests.log', 'a')
DELAY = float(os.environ.get('DALI_DELAY', '2.5'))
counter = [0]
stats = {'ok': 0, 'notfound': 0, 'failed': 0}
last_request = [0.0]
fail_streak = [0]
SCALE = float(os.environ.get('DALI_BACKOFF_SCALE', '1'))     # 1 in the package; the self-test of the script uses a small number
BACKOFF = [30 * SCALE, 60 * SCALE, 120 * SCALE, 300 * SCALE, 600 * SCALE, 900 * SCALE]
ARCHIVE = os.environ.get('DALI_ARCHIVE_HOST', 'https://web.archive.org')       # self-test only: a local mock of the archive


def log(msg):
    logf.write(msg + '\n')
    logf.flush()
    print(msg, flush=True)


def once(url, dest=None, head=False, max_time=240):
    """One curl attempt (follows redirects). Returns (rc, http, final_url, content_type, size, retry_after)."""
    wait = DELAY - (time.time() - last_request[0])
    if wait > 0:
        time.sleep(wait)
    last_request[0] = time.time()
    cmd = ['curl', '-sS', '-L', '-A', UA, '--connect-timeout', '25', '--max-time', str(max_time), '--max-redirs', '8', '-D', '-', '-w', '\n@@%{http_code}|%{url_effective}|%{content_type}|%{size_download}', '-o', dest or '/dev/null']
    if head:
        cmd.insert(1, '-I')
        cmd.remove('-o'); cmd.remove(dest or '/dev/null')
    cmd.append(url)
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=max_time + 30)
        text = r.stdout or ''
        tail = text.split('@@')[-1] if '@@' in text else '||'
        parts = (tail.split('|') + ['', '', '', ''])[:4]
        ra = re.findall(r'(?im)^retry-after:\s*(\d+)', text)
        return (r.returncode, parts[0], parts[1], parts[2], int(parts[3]) if parts[3].strip().isdigit() else 0, int(ra[-1]) if ra else 0, (r.stderr or '').strip()[:120])
    except Exception as e:  # noqa: BLE001
        return (-1, '000', '', '', 0, 0, f'EXCEPTION {type(e).__name__} {e}')


def curl(url, dest=None, head=False, max_time=240):
    """Polite request: serial, delayed, retried with back-off. Returns (status, http, final_url, content_type, size) where status is OK (2xx / 3xx resolved), NOT-FOUND (404 / 403 / 410 after redirects) or FAILED."""
    res = None
    for attempt in range(len(BACKOFF) + 1):
        counter[0] += 1
        rid = counter[0]
        rc, http, final, ctype, size, retry_after, err = once(url, dest, head, max_time)
        log(f'REQ {rid}\t{"HEAD" if head else "GET"}\t{url}\trc={rc}\thttp={http}\tfinal={final}\ttype={ctype}\tbytes={size}\tattempt={attempt}\t{err}')
        transient = rc != 0 or http in ('429', '500', '502', '503', '504', '000')
        if not transient:
            fail_streak[0] = 0
            status = 'OK' if http.startswith('2') else ('NOT-FOUND' if http in ('404', '403', '410', '400') else 'OK')
            stats['ok' if status == 'OK' else 'notfound'] += 1
            return (status, http, final, ctype, size)
        fail_streak[0] += 1
        if attempt == len(BACKOFF):
            break
        pause = max(BACKOFF[attempt], retry_after)
        if fail_streak[0] >= 10:
            pause = max(pause, 900 * SCALE)
        log(f'BACKOFF {pause}s after rc={rc} http={http}')
        time.sleep(pause)
    stats['failed'] += 1
    return ('FAILED', http, final, ctype, size)


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


VERSIONS = ['1.0', '1.1', '1.2', '2.0', '2.1', '2.2', '2.3', '2.4', '2.4.1', '2.4.2', '2.4.3', '2.4.4', '2.4.5', '3.0', '3.1', '3.2', '3.3']
NAMES = []
for v in VERSIONS:
    NAMES += [f'DaliLite_{v}.tar.gz', f'DaliLite.v{v}.tar.gz', f'DaliLite-{v}.tar.gz', f'DaliLite_{v}.tgz']
NAMES += ['DaliLite.tar.gz', 'DaliLite.zip', 'DaliLite.tgz', 'DaliLite.tar.Z', 'dalilite.tar.gz', 'DaliLite_old.tar.gz', 'DaliLite.v2.tar.gz', 'DaliLite.v3.tar.gz', 'DaliLite.v4.tar.gz']
BASES = ['http://www.ebi.ac.uk/~holm/DaliLite/', 'http://www.ebi.ac.uk/~holm/', 'http://www.ebi.ac.uk/dali/DaliLite/', 'http://www.ebi.ac.uk/DaliLite/', 'http://www.ebi.ac.uk/dali/', 'http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/',
         'http://ekhidna.biocenter.helsinki.fi/dali_lite/', 'http://ekhidna.biocenter.helsinki.fi/dali/', 'http://ekhidna.biocenter.helsinki.fi/dali/DaliLite/', 'http://ekhidna.biocenter.helsinki.fi/DaliLite/',
         'http://ekhidna2.biocenter.helsinki.fi/dali/', 'http://ekhidna2.biocenter.helsinki.fi/dali/downloads/', 'http://ekhidna2.biocenter.helsinki.fi/dali_lite/downloads/', 'http://www.bioinformatics.org/dali/',
         'http://www.bioinformatics.org/project/', 'ftp://ftp.ebi.ac.uk/pub/databases/dali/', 'http://ftp.ebi.ac.uk/pub/databases/dali/', 'ftp://ftp.ebi.ac.uk/pub/software/dali/', 'http://ftp.ebi.ac.uk/pub/software/dali/', 'http://www.embl-ebi.ac.uk/DaliLite/']
YEARS = list(range(2002, 2014))
if os.environ.get('DALI_QUICK'):                       # self-test only
    NAMES, BASES, YEARS = NAMES[:2], BASES[:2], YEARS[:2]
log(f'PLAN run 2: delay {DELAY}s, names={len(NAMES)} bases={len(BASES)} years={len(YEARS)} -> {len(NAMES) * len(BASES)} name URLs, 1 probe each, {len(YEARS) - 1} more for each URL that answers')

ARCHIVE_EXT = re.compile(r'\.(tar\.gz|tgz|tar|zip|tar\.Z|Z|gz|tar\.bz2|bz2)$', re.I)
kept = {}
tried = set()
candidates = {}        # url -> why


def snap_url(year, url):
    return f'{ARCHIVE}/web/{year}id_/{url}'


def fetch_text(url, tag):
    dest = f'{out}/listings/{hashlib.sha1(url.encode()).hexdigest()[:10]}_{tag}.html'
    status, http, final, ctype, size = curl(url, dest=dest, max_time=120)
    text = ''
    if status == 'OK':
        try:
            text = open(dest, 'rb').read(2_000_000).decode('latin-1')
        except OSError:
            pass
    return status, final, text


def links_of(text, base_final):
    res = []
    m = re.search(r'/web/\d+id_/(https?://[^ ]+)$', base_final or '')
    origin = m.group(1) if m else ''
    for h in re.findall(r'(?:href|src)\s*=\s*["\']([^"\']+)["\']', text, re.I):
        h = h.split('#')[0]
        if not h or h.startswith(('mailto:', 'javascript:')):
            continue
        h = re.sub(r'^(?:https?://web\.archive\.org)?/web/\d+(?:id_|if_|im_)?/', '', h)
        res.append(urllib.parse.urljoin(origin, h) if origin and not re.match(r'^[a-z]+://', h) else h)
    return res


# (A) the leads
LEADS = ['http://www.ebi.ac.uk/~holm/DaliLite', 'http://www.ebi.ac.uk/~holm/DaliLite/', 'http://www.ebi.ac.uk/dali/DaliLite/', 'http://www.ebi.ac.uk/DaliLite/', 'http://www.ebi.ac.uk/~holm/']
pages = {}
for u in LEADS:
    for y in YEARS:
        status, final, text = fetch_text(snap_url(y, u), f'{y}')
        key = (u, final)
        if status == 'OK' and key not in pages:
            pages[key] = text
            log(f'LEAD-PAGE {u} year {y} -> {final} ({len(text)} bytes)')
            for h in links_of(text, final):
                if ARCHIVE_EXT.search(h.split('?')[0]) or re.search(r'dali', h, re.I) and re.search(r'download|\.tar|\.zip|\.gz', h, re.I):
                    candidates.setdefault(h, f'link on {u} ({y})')
log(f'LEADS: {len(pages)} distinct captured pages, {len(candidates)} candidate links: {sorted(candidates)[:30]}')

# (B) CDX
for b in BASES:
    q = urllib.parse.quote(b, safe=':/~')
    dest = f'{out}/cdx/{hashlib.sha1(b.encode()).hexdigest()[:10]}.txt'
    status, http, final, ctype, size = curl(f'{ARCHIVE}/cdx/search/cdx?url={q}&matchType=prefix&output=txt&fl=timestamp,original,statuscode,mimetype&collapse=urlkey&limit=5000', dest=dest, max_time=180)
    if status == 'OK':
        for line in open(dest, errors='replace'):
            f = line.split()
            if len(f) >= 2 and (ARCHIVE_EXT.search(f[1].split('?')[0]) or re.search(r'dali', f[1], re.I) and ARCHIVE_EXT.search(f[1])):
                candidates.setdefault(f[1], f'CDX {b}')
        log(f'CDX {b}: {sum(1 for _ in open(dest, errors="replace"))} lines')


def probe(url, year):
    """HEAD at web/<year>id_/<url>; the body only after a 2xx that is not HTML; kept only if the magic bytes are an archive's."""
    key = (url, year)
    if key in tried:
        return None
    tried.add(key)
    snap = snap_url(year, url)
    status, http, final, ctype, size = curl(snap, head=True, max_time=90)
    if status != 'OK' or not http.startswith('2') or 'html' in ctype.lower():
        return status if status != 'OK' else http
    tmp = f'{out}/found/.tmp'
    status, http, final, ctype, size = curl(snap, dest=tmp, max_time=900)
    m = magic(tmp)
    if status == 'OK' and http.startswith('2') and size > 0 and m in ('gzip', 'zip', 'compress', 'bzip2', 'tar'):
        sha = hashlib.sha256(open(tmp, 'rb').read()).hexdigest()
        fresh = sha not in kept
        if fresh:
            name = f'{len(kept) + 1:02d}_{re.sub(r"[^A-Za-z0-9._-]", "_", os.path.basename(urllib.parse.urlparse(url).path))}'
            os.replace(tmp, f'{out}/found/{name}')
            kept[sha] = (name, snap, final, m, size)
        log(f'KEPT\t{"new" if fresh else "duplicate"}\t{snap}\tfinal={final}\t{m}\t{size}\tsha256={sha}')
    elif os.path.exists(tmp):
        os.remove(tmp)
    return http


# (C) one probe per (base, name) + the candidates of (A) / (B)
urls = list(dict.fromkeys(list(candidates) + [b + n for b in BASES for n in NAMES]))
answers = {}
for u in urls:
    answers[u] = probe(u, 2008)
interesting = [u for u in urls if u in candidates or answers.get(u) not in (None, '404', '403', '410', '400', 'NOT-FOUND', 'FAILED')]
failed = [u for u in urls if answers.get(u) == 'FAILED']
log(f'STAGE C {len(urls)} URLs: {len(interesting)} interesting, {len(failed)} FAILED (no answer even after retries: NOT evidence of absence)')
# (D) all years for the interesting ones; the failed ones are asked once more
for u in interesting:
    for y in YEARS:
        probe(u, y)
for u in failed:
    tried.discard((u, 2008))
    probe(u, 2008)
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
report.write(f'SUMMARY requests={counter[0]} ok={stats["ok"]} notfound={stats["notfound"]} failed={stats["failed"]} urls_probed={len(urls)} interesting={len(interesting)} failed_urls={len(failed)} archives_kept={len(kept)} captured_lead_pages={len(pages)}\n')
report.flush()
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
json.dump({'requests': counter[0], **stats, 'urls': len(urls), 'interesting': interesting, 'failed': failed, 'kept': len(kept), 'candidates': candidates}, open(f'{out}/summary.json', 'w'), indent=1)
log(f'DONE requests={counter[0]} kept={len(kept)} stats={stats}')
