#!/usr/bin/env python3
"""RUN 3 (B09 WP87, F378). Probes Internet Archive SNAPSHOTS for older DaliLite releases, keeps what is a real gzip / zip / compress / bzip2 archive (magic bytes), reads CHANGES / README / INSTALL / LICENCE / COPYING / the `DaliLite` script and the
Perl scripts as TEXT ONLY (python tarfile / zipfile in memory: nothing is extracted to disk, compiled or run) and reports whether any release prints `Query` / `Sbjct` to standard output.

History. Run 1 (F371) was defective (parallel requests -> HTTP 429 / refused connections counted as "not archived"). Run 2 asked politely and answered stage C completely (1539 HEAD probes = 404 for 77 release names x 20 bases; kept DaliLite_3.3 and
DaliLite.v4.1) but was stopped by the user after ~17 h: the candidate filter admitted every archive under the bases (7 multi-GB DALI database files, 0.35 MB/s, 900 s timeouts retried six times), the years of one URL resolved to the same capture,
the CDX of ekhidna.../dali/ was cut at limit=5000, nothing reached the console and no summary was written on interruption (returned/A2-12b-s6b-run2-STOPPED-notes.md, defects 1 - 9). Run 3 repairs all nine:
 1 candidate filter: the BASENAME must look like a DaliLite release (`dali.?lite` + archive extension; a looser `^dali...` + archive extension only if the size is known to be small); every rejected URL is logged once (REJECT);
 2 size guard: `--max-filesize 20M` (curl rc 63) AND a poll of the growing file (kills the GET above 20 MB even when the server sends no length), HEAD Content-Length above the cap skips the GET; final status TOO-LARGE, never retried;
 3 a timeout with data flowing (rc 28, HTTP 2xx, bytes > 0) is slow/huge, not unreachable: final status TIMEOUT, no retry; HEAD rc 28: at most one retry; other partial transfers (rc 92 / 56 / 18): at most one retry;
 4 stage D only for URLs whose probe KEPT an archive;
 5 captures are de-duplicated: stage D asks the CDX for the distinct content digests of the exact URL and fetches each distinct capture once (by its timestamp); without CDX rows the years are swept but a capture timestamp already fetched is skipped;
 6 CDX with showResumeKey / resumeKey pagination (limit 5000 per page), two queries per base (everything except 404 / query strings, and archive extensions only), `TRUNCATED` in the log if the page cap is reached;
 7 progress: every log line goes to the console (in_container.sh tees it to out/fetch.log) and a PROGRESS line every 25 requests;
 8 summary: report.txt / summary.json (SUMMARY line) rewritten after every stage, every 25 requests and after every kept archive; each kept archive is inspected the moment it is kept; SIGTERM / SIGINT write the summary and exit (run.sh uses `docker run --init`);
 9 run.sh packs out/ without out/found (third-party archives: academic licence, public repository, GitHub 100 MB limit) and without .tmp files.
New in run 3 (what the run-2 data shows is still missing): (a) PAGE LEADS: the 200 text pages that the CDX lists under the bases and whose path names DaliLite / download / licence / README (e.g. `~holm/DaliLite/dalilicence` 2001, `dali_lite/downloads` 2005,
`downloads/download.html` 2017) are fetched at their exact capture and their links followed one level (they may name release archives the name list does not know); (b) the Helsinki directory layout `dali_lite/downloads/v1 ... v5/`, `old/`, `/downloads/`,
`/software/` as extra bases and the release names 4.0 / 4.1 / 5.0, `DaliLite.v1` / `.v5`; (c) the 1539 (base, name) pairs that run 2 answered 404 are not asked again (run2_known_404.txt; DALI_REPROBE=1 asks them again).
Requests are serial, >= DALI_DELAY seconds apart (default 2.5), retried with back-off 30 s ... 15 min on 429 / 5xx / connection faults; every request is OK / NOT-FOUND / FAILED / TOO-LARGE / TIMEOUT and FAILED is never read as "not archived".
Neutral User-Agent `Seqeron-review/1.0`; no account data. Usage: fetch_dalilite_snapshots.py OUTDIR"""
import hashlib, io, json, os, re, signal, subprocess, sys, tarfile, time, urllib.parse, zipfile

UA = 'Seqeron-review/1.0'
out = sys.argv[1]
for d in ('found', 'inspect', 'listings', 'cdx', 'pages'):
    os.makedirs(f'{out}/{d}', exist_ok=True)
logf = open(f'{out}/requests.log', 'a')
DELAY = float(os.environ.get('DALI_DELAY', '2.5'))
SCALE = float(os.environ.get('DALI_BACKOFF_SCALE', '1'))     # 1 in the package; the self-test uses a small number
BACKOFF = [30 * SCALE, 60 * SCALE, 120 * SCALE, 300 * SCALE, 600 * SCALE, 900 * SCALE]
ARCHIVE = os.environ.get('DALI_ARCHIVE_HOST', 'https://web.archive.org')       # self-test only: a local mock of the archive
CAP = int(os.environ.get('DALI_MAX_BYTES', str(20 * 1024 * 1024)))               # largest file that is downloaded (the real releases are 0.8 - 1.4 MB)
CDX_LIMIT = int(os.environ.get('DALI_CDX_LIMIT', '5000'))
CDX_MAX_PAGES = int(os.environ.get('DALI_CDX_PAGES', '12'))
GET_TIME = int(os.environ.get('DALI_GET_TIME', '600'))
KNOWN_404 = os.environ.get('DALI_KNOWN_404', os.path.join(os.path.dirname(os.path.abspath(out)), 'run2_known_404.txt'))
T0 = time.time()
counter = [0]
stats = {'ok': 0, 'notfound': 0, 'failed': 0, 'toolarge': 0, 'timeout': 0}
last_request = [0.0]
fail_streak = [0]
stage = ['init']
report_lines = {}     # archive name -> report line
state = {'cdx_truncated': [], 'rejected': 0, 'pages': 0, 'skipped_known': 0, 'same_capture': 0}
kept = {}             # sha256 -> (name, snap, final, magic, size)
tried = set()
captures_done = {}    # (url, ts) -> result
candidates = {}       # url -> why
cdx_ok = [False]      # did the last CDX query get an answer
rejected = set()
pages_seen = {}       # (url, ts) -> text
version_mentions = {}


class Stopped(BaseException):
    pass


def on_signal(signum, _frame):
    raise Stopped(signal.Signals(signum).name)


signal.signal(signal.SIGTERM, on_signal)
signal.signal(signal.SIGINT, on_signal)


def log(msg):
    logf.write(msg + '\n')
    logf.flush()
    print(msg, flush=True)


def write_summary(final=False):
    """report.txt (SUMMARY line + one line per kept archive), summary.json, found_index.tsv: rewritten as things happen."""
    head = (f'SUMMARY{" FINAL" if final else " (in progress)"} stage={stage[0]} elapsed_s={int(time.time() - T0)} requests={counter[0]} ok={stats["ok"]} notfound={stats["notfound"]} failed={stats["failed"]} '
            f'toolarge={stats["toolarge"]} timeout={stats["timeout"]} archives_kept={len(kept)} rejected_candidates={state["rejected"]} text_pages={state["pages"]} skipped_known_404={state["skipped_known"]} '
            f'same_capture_skips={state["same_capture"]} cdx_truncated={len(state["cdx_truncated"])}')
    with open(f'{out}/report.txt', 'w') as f:
        f.write(head + '\n')
        for n in sorted(report_lines):
            f.write(report_lines[n] + '\n')
    with open(f'{out}/summary.json', 'w') as f:
        json.dump({'final': final, 'stage': stage[0], 'elapsed_s': int(time.time() - T0), 'requests': counter[0], **stats, **{k: v for k, v in state.items()}, 'kept': len(kept), 'candidates': candidates,
                   'version_mentions': {k: sorted(v) for k, v in version_mentions.items()}}, f, indent=1)
    with open(f'{out}/found_index.tsv', 'w') as f:
        f.write('name\tsha256\tbytes\tmagic\tsnapshot\tfinal\n')
        for sha, (name, snap, final_url, m, size) in sorted(kept.items(), key=lambda kv: kv[1][0]):
            f.write(f'{name}\t{sha}\t{size}\t{m}\t{snap}\t{final_url}\n')


def pace():
    wait = DELAY - (time.time() - last_request[0])
    if wait > 0:
        time.sleep(wait)
    last_request[0] = time.time()


def once(url, dest=None, head=False, max_time=240, cap=None):
    """One curl attempt (follows redirects). Returns dict(rc, http, final, ctype, size, retry_after, clen, err, killed). A GET is polled: it is killed when the file grows above `cap`."""
    pace()
    cmd = ['curl', '-sS', '-L', '-A', UA, '--connect-timeout', '25', '--max-time', str(max_time), '--max-redirs', '8', '-D', '-', '-w', '\n@@%{http_code}|%{url_effective}|%{content_type}|%{size_download}']
    if head:
        cmd.insert(1, '-I')
    else:
        cmd += ['--max-filesize', str(cap or CAP), '-o', dest or '/dev/null']
    cmd.append(url)
    killed = None
    so = se = ''
    p = None
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, errors='replace')
        t0 = time.time()
        while True:
            try:
                so, se = p.communicate(timeout=0.5)
                break
            except subprocess.TimeoutExpired:
                if not head and dest and os.path.exists(dest) and os.path.getsize(dest) > (cap or CAP):
                    killed = 'TOO-LARGE'
                elif time.time() - t0 > max_time + 20:
                    killed = 'TIMEOUT'
                if killed:
                    p.kill()
                    so, se = p.communicate()
                    break
    except OSError as e:
        return dict(rc=-1, http='000', final='', ctype='', size=0, retry_after=0, clen=0, err=f'EXCEPTION {type(e).__name__} {e}', killed=None)
    finally:
        if p is not None and p.poll() is None:
            p.kill()
            p.wait()
    so = so or ''
    tail = so.split('@@')[-1] if '@@' in so else ''
    parts = (tail.split('|') + ['', '', '', ''])[:4]
    http = parts[0] if parts[0].strip() else (re.findall(r'(?m)^HTTP/\S+ (\d{3})', so) or ['000'])[-1]
    size = int(parts[3]) if parts[3].strip().isdigit() else (os.path.getsize(dest) if dest and os.path.exists(dest) else 0)
    ra = re.findall(r'(?im)^retry-after:\s*(\d+)', so)
    cl = re.findall(r'(?im)^content-length:\s*(\d+)', so)
    rc = p.returncode if p is not None else -1
    if killed == 'TOO-LARGE':
        rc = 63
    elif killed == 'TIMEOUT':
        rc = 28
    return dict(rc=rc, http=http, final=parts[1], ctype=parts[2], size=size, retry_after=int(ra[-1]) if ra else 0, clen=int(cl[-1]) if cl else 0, err=(se or '').strip()[:120], killed=killed)


def fetch(url, dest=None, head=False, max_time=240, cap=None):
    """Polite request: serial, delayed, retried with back-off. Returns dict(status, http, final, ctype, size, clen, ts); status is OK (2xx / 3xx resolved), NOT-FOUND (404 / 403 / 410 / 400),
    TOO-LARGE (above the cap, final), TIMEOUT (rc 28 with data flowing, final) or FAILED (no answer even after the retries)."""
    r = None
    partial_retries = 0
    head_timeouts = 0
    for attempt in range(len(BACKOFF) + 1):
        counter[0] += 1
        rid = counter[0]
        r = once(url, dest, head, max_time, cap)
        log(f'REQ {rid}\t{"HEAD" if head else "GET"}\t{url}\trc={r["rc"]}\thttp={r["http"]}\tfinal={r["final"]}\ttype={r["ctype"]}\tbytes={r["size"]}\tattempt={attempt}\t{r["err"]}')
        if counter[0] % 25 == 0:
            log(f'PROGRESS stage={stage[0]} requests={counter[0]} elapsed={int(time.time() - T0)}s ok={stats["ok"]} notfound={stats["notfound"]} failed={stats["failed"]} kept={len(kept)}')
            write_summary()
        ts = (re.search(r'/web/(\d{14})id_/', r['final'] or '') or [None, ''])[1]
        res = dict(status='', http=r['http'], final=r['final'], ctype=r['ctype'], size=r['size'], clen=r['clen'], ts=ts)
        data_flowing = (not head) and r['http'].startswith('2') and r['size'] > 0
        if r['rc'] == 63:
            res['status'] = 'TOO-LARGE'
        elif r['rc'] == 28 and data_flowing:
            res['status'] = 'TIMEOUT'
        elif r['rc'] != 0 and data_flowing:
            partial_retries += 1                       # rc 92 / 56 / 18: the connection broke with data flowing: one more try, not seven
            if partial_retries > 1:
                res['status'] = 'FAILED'
        elif r['rc'] == 28 and head:
            head_timeouts += 1                         # a HEAD that times out on a redirect: at most one retry
            if head_timeouts > 1:
                res['status'] = 'FAILED'
        if not res['status']:
            transient = r['rc'] != 0 or r['http'] in ('429', '500', '502', '503', '504', '000')
            if not transient:
                res['status'] = 'OK' if r['http'].startswith('2') else ('NOT-FOUND' if r['http'] in ('404', '403', '410', '400') else 'OK')
        if res['status']:
            if res['status'] in ('OK', 'NOT-FOUND'):
                fail_streak[0] = 0
            key = {'OK': 'ok', 'NOT-FOUND': 'notfound', 'FAILED': 'failed', 'TOO-LARGE': 'toolarge', 'TIMEOUT': 'timeout'}[res['status']]
            stats[key] += 1
            if res['status'] in ('FAILED', 'TOO-LARGE', 'TIMEOUT'):
                log(f'RESULT {res["status"]} {url} rc={r["rc"]} http={r["http"]} bytes={r["size"]}')
            return res
        fail_streak[0] += 1
        if attempt == len(BACKOFF):
            break
        pause = max(BACKOFF[attempt], r['retry_after'])
        if fail_streak[0] >= 10:
            pause = max(pause, 900 * SCALE)
        log(f'BACKOFF {pause}s after rc={r["rc"]} http={r["http"]}')
        time.sleep(pause)
    stats['failed'] += 1
    res = dict(status='FAILED', http=r['http'], final=r['final'], ctype=r['ctype'], size=r['size'], clen=r['clen'], ts='')
    log(f'RESULT FAILED {url} rc={r["rc"]} http={r["http"]}')
    return res


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
NEW_NAMES = []         # run 3: the releases after 3.3 and the major-version names
for v in ('4.0', '4.1', '5.0'):
    NEW_NAMES += [f'DaliLite_{v}.tar.gz', f'DaliLite.v{v}.tar.gz', f'DaliLite-{v}.tar.gz', f'DaliLite_{v}.tgz']
NEW_NAMES += ['DaliLite.v1.tar.gz', 'DaliLite.v5.tar.gz']
OLD_BASES = ['http://www.ebi.ac.uk/~holm/DaliLite/', 'http://www.ebi.ac.uk/~holm/', 'http://www.ebi.ac.uk/dali/DaliLite/', 'http://www.ebi.ac.uk/DaliLite/', 'http://www.ebi.ac.uk/dali/', 'http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/',
             'http://ekhidna.biocenter.helsinki.fi/dali_lite/', 'http://ekhidna.biocenter.helsinki.fi/dali/', 'http://ekhidna.biocenter.helsinki.fi/dali/DaliLite/', 'http://ekhidna.biocenter.helsinki.fi/DaliLite/',
             'http://ekhidna2.biocenter.helsinki.fi/dali/', 'http://ekhidna2.biocenter.helsinki.fi/dali/downloads/', 'http://ekhidna2.biocenter.helsinki.fi/dali_lite/downloads/', 'http://www.bioinformatics.org/dali/',
             'http://www.bioinformatics.org/project/', 'ftp://ftp.ebi.ac.uk/pub/databases/dali/', 'http://ftp.ebi.ac.uk/pub/databases/dali/', 'ftp://ftp.ebi.ac.uk/pub/software/dali/', 'http://ftp.ebi.ac.uk/pub/software/dali/', 'http://www.embl-ebi.ac.uk/DaliLite/']
# run 2 found the release at ekhidna.../dali_lite/downloads/v3/DaliLite.tar.gz: the version directories were not probed
NEW_BASES = [f'http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/{d}/' for d in ('v1', 'v2', 'v3', 'v4', 'v5', 'old')] + ['http://ekhidna.biocenter.helsinki.fi/downloads/', 'http://ekhidna.biocenter.helsinki.fi/software/']
BASES = OLD_BASES + NEW_BASES
YEARS = list(range(2002, 2014))
LEADS = ['http://www.ebi.ac.uk/~holm/DaliLite', 'http://www.ebi.ac.uk/~holm/DaliLite/', 'http://www.ebi.ac.uk/dali/DaliLite/', 'http://www.ebi.ac.uk/DaliLite/', 'http://www.ebi.ac.uk/~holm/',
         'http://ekhidna.biocenter.helsinki.fi/dali_lite/', 'http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/', 'http://ekhidna.biocenter.helsinki.fi/downloads', 'http://ekhidna.biocenter.helsinki.fi/software']
if os.environ.get('DALI_QUICK'):                       # self-test only
    NAMES, NEW_NAMES, BASES, OLD_BASES, NEW_BASES, YEARS, LEADS = NAMES[:2], NEW_NAMES[:1], BASES[:3] + BASES[5:6], BASES[:3], BASES[5:6], YEARS[:2], LEADS[:2]

ARCHIVE_EXT = re.compile(r'\.(tar\.gz|tgz|tar|zip|tar\.Z|Z|gz|tar\.bz2|bz2)$', re.I)
STRICT_NAME = re.compile(r'dali.?lite', re.I)
LOOSE_NAME = re.compile(r'^dali', re.I)
KNOWN_DB = re.compile(r'^(DAT|dump|GTG|AF-Digest|benchmark|pdb_and_scope|pdb90|DaliQuiz|dalicon)\b', re.I)       # database / data files of the DALI server: never a release
TEXT_EXT = re.compile(r'\.(gif|png|jpe?g|ico|css|js|pdf|ps|pdb|brk|swf|xml|tar|gz|tgz|zip|Z|bz2)$', re.I)
PAGE_PATH = re.compile(r'dali.?lite|download|licen[cs]e|readme|install|software|dali_lite', re.I)
print_plan = (f'PLAN run 3: delay {DELAY}s, cap {CAP} bytes, bases={len(BASES)} ({len(NEW_BASES)} new), names={len(NAMES)} (+{len(NEW_NAMES)} new), years={len(YEARS)}, '
              f'CDX limit {CDX_LIMIT} x {CDX_MAX_PAGES} pages, known 404 file {KNOWN_404}')


def candidate_kind(url):
    """'strict' (DaliLite release name), 'loose' (starts with dali) or None; logs every rejected archive-like URL once."""
    path = urllib.parse.unquote(url.split('?')[0].split('#')[0])
    base = os.path.basename(path.rstrip('/'))
    if not ARCHIVE_EXT.search(base):
        return None
    if STRICT_NAME.search(base):
        return 'strict'
    if LOOSE_NAME.search(base) and not KNOWN_DB.search(base):
        return 'loose'
    if url not in rejected:
        rejected.add(url)
        state['rejected'] += 1
        log(f'REJECT {url} (basename {base!r} is not a DaliLite release name)')
    return None


def snap_url(when, url):
    return f'{ARCHIVE}/web/{when}id_/{url}'


def fetch_text(url, tag, cap=5_000_000):
    dest = f'{out}/listings/{hashlib.sha1(url.encode()).hexdigest()[:10]}_{tag}.html'
    r = fetch(url, dest=dest, max_time=120, cap=cap)
    text = ''
    if r['status'] == 'OK':
        try:
            text = open(dest, 'rb').read(2_000_000).decode('latin-1')
        except OSError:
            pass
    return r, text


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


def note_versions(src, text):
    vs = set(re.findall(r'DaliLite[ _.-]*v?\.?\s*(\d+(?:\.\d+)+)', text))
    if vs:
        version_mentions[src] = vs


def follow_links(text, final, why, page_queue):
    for h in links_of(text, final):
        if candidate_kind(h):
            if candidates.setdefault(h, why) is why:
                log(f'CANDIDATE {h} <- {why}')
        elif page_queue is not None and not TEXT_EXT.search(h.split('?')[0]) and '?' not in h and PAGE_PATH.search(h) and re.match(r'https?://', h):
            page_queue.setdefault(h, why)


# ---- CDX (resumeKey pagination) -----------------------------------------------------------------------------------------------------------------------------------------------------------------
def flt(field, regex):
    return 'filter=' + urllib.parse.quote(f'{field}:{regex}', safe='')


def cdx(params, tag, label):
    """All rows (timestamp, statuscode, mimetype, digest, original) of a CDX query, paginated by resumeKey; TRUNCATED when the page cap stops it. Raw text saved as cdx/<tag>.txt."""
    rows = []
    key = ''
    raw = open(f'{out}/cdx/{tag}.txt', 'w')
    ok = False
    cdx_ok[0] = False
    for page in range(CDX_MAX_PAGES):
        q = f'{ARCHIVE}/cdx/search/cdx?{params}&output=txt&fl=timestamp,statuscode,mimetype,digest,original&limit={CDX_LIMIT}&showResumeKey=true'
        if key:
            q += '&resumeKey=' + urllib.parse.quote(key, safe='')
        dest = f'{out}/cdx/.page'
        r = fetch(q, dest=dest, max_time=180, cap=60_000_000)
        if r['status'] != 'OK':
            log(f'CDX {label}: page {page + 1} {r["status"]}')
            break
        ok = True
        cdx_ok[0] = True
        text = open(dest, errors='replace').read()
        lines = text.split('\n')
        nxt = ''
        if '' in lines:                                  # rows, an empty line, the resume key
            i = lines.index('')
            rest = [x for x in lines[i + 1:] if x.strip()]
            nxt = rest[0].strip() if rest else ''
            lines = lines[:i]
        n = 0
        for line in lines:
            f = line.split(None, 4)
            if len(f) == 5:
                rows.append((f[0], f[1], f[2], f[3], f[4]))
                raw.write(line + '\n')
                n += 1
        if not nxt or n == 0:
            key = ''
            break
        key = nxt
    raw.close()
    if os.path.exists(f'{out}/cdx/.page'):
        os.remove(f'{out}/cdx/.page')
    if key:
        state['cdx_truncated'].append(label)
        log(f'CDX {label}: TRUNCATED after {CDX_MAX_PAGES} pages ({len(rows)} rows; resumeKey still returned): entries beyond the cap were NOT seen')
    else:
        log(f'CDX {label}: {len(rows)} rows{"" if ok else " (NO ANSWER: not evidence of absence)"}')
    return rows


def get_capture(url, when):
    """GET one capture (`when` = 14-digit timestamp or a year); the body is kept only if its magic bytes are an archive's and the size is within the cap."""
    snap = snap_url(when, url)
    tmp = f'{out}/found/.tmp'
    r = fetch(snap, dest=tmp, max_time=GET_TIME)
    m = magic(tmp)
    result = r['status'] if r['status'] != 'OK' else r['http']
    if r['status'] == 'OK' and r['http'].startswith('2') and r['size'] > 0 and m in ('gzip', 'zip', 'compress', 'bzip2', 'tar'):
        sha = hashlib.sha256(open(tmp, 'rb').read()).hexdigest()
        fresh = sha not in kept
        if fresh:
            name = f'{len(kept) + 1:02d}_{re.sub(r"[^A-Za-z0-9._-]", "_", os.path.basename(urllib.parse.urlparse(url).path))}'
            os.replace(tmp, f'{out}/found/{name}')
            kept[sha] = (name, snap, r['final'], m, r['size'])
        log(f'KEPT\t{"new" if fresh else "duplicate"}\t{snap}\tfinal={r["final"]}\t{m}\t{r["size"]}\tsha256={sha}')
        if fresh:
            inspect_kept(sha)
            write_summary()
        result = 'KEPT'
    elif os.path.exists(tmp):
        os.remove(tmp)
    if r['ts']:
        captures_done[(url, r['ts'])] = result
    return result, r


def probe(url, year):
    """HEAD at web/<year>id_/<url> (the archive answers with the nearest capture of the exact URL); the body only after a 2xx that is not HTML, within the cap, and for a capture not yet fetched."""
    key = (url, year)
    if key in tried:
        return None
    tried.add(key)
    r = fetch(snap_url(year, url), head=True, max_time=90)
    if r['status'] != 'OK':
        return r['status']
    if not r['http'].startswith('2'):
        return r['http']
    if 'html' in r['ctype'].lower():
        return 'HTML'
    if r['clen'] > CAP:
        stats['toolarge'] += 1
        log(f'RESULT TOO-LARGE {url} (HEAD content-length {r["clen"]} > {CAP})')
        return 'TOO-LARGE'
    if r['ts'] and (url, r['ts']) in captures_done:
        state['same_capture'] += 1
        log(f'SAME-CAPTURE {url} year {year} -> {r["ts"]} already fetched: skipped')
        return captures_done[(url, r['ts'])]
    result, _ = get_capture(url, r['ts'] or year)
    return result


# ---- inspection (as TEXT, in memory) ---------------------------------------------------------------------------------------------------------------------------------------------------------------
TEXT_MEMBER = re.compile(r'(^|/)(CHANGES|README[^/]*|INSTALL|Makefile|DaliLite|LICEN[CS]E|COPYING|COPYRIGHT|[^/]*\.pl)$', re.I)


def members_of(path, m):
    data = open(path, 'rb').read()
    if m == 'compress':
        p = subprocess.Popen(['gzip', '-dc'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        data = p.communicate(data)[0][:400_000_000]
        m = 'tar'
    if m == 'bzip2':
        import bz2
        data = bz2.BZ2Decompressor().decompress(data, 400_000_000)
        m = 'tar'
    names, texts, hits = [], {}, []
    if m in ('gzip', 'tar'):
        with tarfile.open(fileobj=io.BytesIO(data), mode='r:*') as t:
            for mem in t:
                names.append(f'{mem.size}\t{mem.name}')
                if mem.isfile() and mem.size < 3_000_000:
                    raw = t.extractfile(mem).read()
                    if TEXT_MEMBER.search(mem.name):
                        texts[mem.name] = raw.decode('latin-1')
                    if b'\x00' not in raw[:4096]:
                        hits += [f'{mem.name}:{i}: {l.strip()[:140]}' for i, l in enumerate(raw.decode('latin-1').splitlines(), 1) if re.search(r'\b(Query|Sbjct)\b', l)]
    elif m == 'zip':
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for mem in z.namelist():
                info = z.getinfo(mem)
                names.append(f'{info.file_size}\t{mem}')
                if info.file_size < 3_000_000:
                    raw = z.read(mem)
                    if TEXT_MEMBER.search(mem):
                        texts[mem] = raw.decode('latin-1')
                    if b'\x00' not in raw[:4096]:
                        hits += [f'{mem}:{i}: {l.strip()[:140]}' for i, l in enumerate(raw.decode('latin-1').splitlines(), 1) if re.search(r'\b(Query|Sbjct)\b', l)]
    return names, texts, hits


def inspect_kept(sha):
    name, snap, final, m, size = kept[sha]
    line = f'{name}\t{snap}\tfinal={final}\t{m}\tbytes={size}\tsha256={sha}'
    try:
        names, texts, hits = members_of(f'{out}/found/{name}', m)
        d = f'{out}/inspect/{name}'
        os.makedirs(d, exist_ok=True)
        open(d + '/members.txt', 'w').write('\n'.join(names))
        for k, v in texts.items():
            open(d + '/' + re.sub(r'[^A-Za-z0-9._-]', '_', k), 'w').write(v)
        script = [v for k, v in texts.items() if re.search(r'(^|/)DaliLite$', k)]
        chg = '\n'.join(v for k, v in texts.items() if k.endswith('CHANGES'))
        vers = re.findall(r'(?:Changes in )?DaliLite v\.?\s*([0-9.]+)', chg)
        top = sorted({n.split('\t')[1].split('/')[0] for n in names})[:3]
        lic = [k for k in texts if re.search(r'(LICEN[CS]E|COPYING|COPYRIGHT)', k, re.I)]
        lic_head = ' | '.join(re.sub(r'\s+', ' ', texts[lic[0]]).strip()[:160].split('|')) if lic else ''
        line += f'\ttop_dirs={",".join(top)}\tmembers={len(names)}\tCHANGES_versions={",".join(vers) if vers else "?"}\tlicence_file={lic[0] if lic else "-"}\tlicence_head={lic_head!r}'
        line += '\tscript_has_Query=%s\tscript_has_Sbjct=%s' % tuple(str(any(s in sc for sc in script)) for s in ('Query', 'Sbjct'))
        line += f'\tQuery_Sbjct_lines_in_text_members={len(hits)}'
        pr = [l.strip() for sc in script for l in sc.splitlines() if re.search(r'Query|Sbjct|print\s+STDOUT|print\s+"|echo\s+"|\bcat\b|-pairwise', l)][:60]
        open(d + '/pairwise_print_lines.txt', 'w').write('\n'.join(pr))
        open(d + '/Query_Sbjct_hits.txt', 'w').write('\n'.join(hits[:200]))
    except Exception as e:  # noqa: BLE001
        line += f'\tUNREADABLE {type(e).__name__} {e}'
    report_lines[name] = line
    log('INSPECTED ' + line)


# ---- main -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def run():
    log(print_plan)
    known = set()
    if os.path.exists(KNOWN_404) and not os.environ.get('DALI_REPROBE'):
        known = {l.strip() for l in open(KNOWN_404) if l.strip()}
        log(f'KNOWN-404 {len(known)} (base, name) URLs that run 2 answered 404 are not asked again (DALI_REPROBE=1 asks them)')
    write_summary()

    # (A) the leads: pages that were captured; links followed one level
    stage[0] = 'A leads'
    page_queue = {}
    pages = {}
    for u in LEADS:
        for y in YEARS:
            r, text = fetch_text(snap_url(y, u), f'{y}')
            key = (u, r['final'])
            if r['status'] == 'OK' and key not in pages:
                pages[key] = text
                log(f'LEAD-PAGE {u} year {y} -> {r["final"]} ({len(text)} bytes)')
                note_versions(r['final'], text)
                follow_links(text, r['final'], f'link on {u} ({y})', page_queue)
    log(f'LEADS: {len(pages)} distinct captured pages, {len(candidates)} candidate links, {len(page_queue)} page leads: {sorted(candidates)[:30]}')
    write_summary()

    # (B) CDX of every base: (1) everything except 404 and query strings, (2) archive extensions only (a 3xx capture of an archive too); resumeKey pagination
    stage[0] = 'B cdx'
    cdx_pages = {}
    for b in BASES:
        q = urllib.parse.quote(b, safe=':/~')
        tag = hashlib.sha1(b.encode()).hexdigest()[:10]
        rows1 = cdx(f'url={q}&matchType=prefix&collapse=urlkey&{flt("!statuscode", "404")}&{flt("!original", r".*(\?|daliquery).*")}', f'{tag}_all', f'{b} (all, no 404, no query strings)')
        ok1 = cdx_ok[0]
        rows2 = cdx(f'url={q}&matchType=prefix&collapse=urlkey&{flt("original", r".*\.(tar\.gz|tgz|tar|zip|Z|gz|bz2)")}', f'{tag}_archives', f'{b} (archive extensions)')
        if not (ok1 and cdx_ok[0]):   # the filter syntax was refused (HTTP 400) or a query failed: the plain query of run 2, paginated
            rows1 = cdx(f'url={q}&matchType=prefix&collapse=urlkey', f'{tag}_plain', f'{b} (fallback: no filters)')
        for ts, status, mime, digest, orig in rows1 + rows2:
            if candidate_kind(orig):
                candidates.setdefault(orig, f'CDX {b} {ts} {status}')
            if status == '200' and re.match(r'text/(html|plain)', mime) and not TEXT_EXT.search(orig.split('?')[0]) and '?' not in orig and PAGE_PATH.search(urllib.parse.urlparse(orig).path):
                cdx_pages.setdefault((orig, ts), b)
    write_summary()

    # (B2) page leads from the CDX (exact capture) and the pages their links name; saved as text under pages/
    stage[0] = 'B2 pages'
    todo = [(u, ts) for (u, ts) in sorted(cdx_pages)][:150]
    seen_p = set()
    n_pages = 0
    for round_ in range(2):
        nxt = []
        for u, ts in todo:
            if (u, ts) in seen_p:
                continue
            seen_p.add((u, ts))
            r, text = fetch_text(snap_url(ts, u), f'p{ts}')
            if r['status'] != 'OK' or not text:
                continue
            n_pages += 1
            state['pages'] = n_pages
            fn = f'{out}/pages/{n_pages:03d}_{ts}_{re.sub(r"[^A-Za-z0-9._-]", "_", os.path.basename(urllib.parse.urlparse(u).path) or "index")}.txt'
            open(fn, 'w').write(f'{u}\t{ts}\t{r["final"]}\n' + text)
            log(f'PAGE {u} {ts} ({len(text)} bytes) -> {os.path.basename(fn)}')
            note_versions(u, text)
            q2 = {}
            follow_links(text, r['final'], f'link on page {u} {ts}', q2 if round_ == 0 else None)
            nxt += [(h, '2008') for h in q2 if len(nxt) < 60]
        todo = nxt
        if not todo:
            break
    # also the lead-page links (A) as pages, one level
    for h in list(page_queue)[:40]:
        if (h, '2008') in seen_p:
            continue
        seen_p.add((h, '2008'))
        r, text = fetch_text(snap_url('2008', h), 'pl')
        if r['status'] == 'OK' and text:
            n_pages += 1
            state['pages'] = n_pages
            open(f'{out}/pages/{n_pages:03d}_lead_{re.sub(r"[^A-Za-z0-9._-]", "_", os.path.basename(urllib.parse.urlparse(h).path) or "index")}.txt', 'w').write(f'{h}\t2008\t{r["final"]}\n' + text)
            log(f'PAGE {h} lead ({len(text)} bytes)')
            note_versions(h, text)
            follow_links(text, r['final'], f'link on lead page {h}', None)
    write_summary()

    # (C) one HEAD probe per (base, name) that run 2 did not answer, plus every candidate
    stage[0] = 'C probes'
    matrix = [b + n for b in OLD_BASES for n in NEW_NAMES] + [b + n for b in NEW_BASES for n in NAMES + NEW_NAMES]
    if os.environ.get('DALI_REPROBE'):
        matrix += [b + n for b in OLD_BASES for n in NAMES]
    skipped = [u for u in list(candidates) + matrix if u in known]
    state['skipped_known'] = len(skipped)
    urls = [u for u in dict.fromkeys(list(candidates) + matrix) if u not in known]
    log(f'STAGE C: {len(urls)} URLs to probe ({len(candidates)} candidates, {len(skipped)} (base, name) pairs skipped as answered 404 in run 2)')
    answers = {}
    for i, u in enumerate(urls, 1):
        answers[u] = probe(u, 2008)
    kept_urls = [u for u in urls if answers.get(u) == 'KEPT']
    failed = [u for u in urls if answers.get(u) == 'FAILED']
    toolarge = [u for u in urls if answers.get(u) in ('TOO-LARGE', 'TIMEOUT')]
    log(f'STAGE C {len(urls)} URLs: {len(kept_urls)} kept an archive, {len(toolarge)} too large / timed out, {len(failed)} FAILED (no answer even after retries: NOT evidence of absence)')
    write_summary()

    # (D) only for the URLs that kept an archive: every DISTINCT capture (by CDX content digest) once; the FAILED ones are asked once more
    stage[0] = 'D distinct captures'
    for u in kept_urls[:40]:
        rows = cdx(f'url={urllib.parse.quote(u, safe=":/~")}&matchType=exact&{flt("statuscode", "2..")}', hashlib.sha1(u.encode()).hexdigest()[:10] + '_exact', f'{u} (exact, 2xx)')
        done_digests = {dg for (ts, st, mi, dg, o) in rows if (u, ts) in captures_done}
        if not rows:
            for y in YEARS:                                   # CDX gave nothing: years, but a capture timestamp already fetched is skipped (probe)
                probe(u, y)
            continue
        for ts, st, mi, dg, o in rows:
            if dg in done_digests or (u, ts) in captures_done:
                continue
            done_digests.add(dg)
            get_capture(u, ts)
    stage[0] = 'D2 failed asked again'
    for u in failed:
        tried.discard((u, 2008))
        probe(u, 2008)
    log(f'KEPT {len(kept)} distinct archives')


try:
    run()
    stage[0] = 'done'
    write_summary(final=True)
    log(f'DONE requests={counter[0]} kept={len(kept)} stats={stats}')
except Stopped as e:
    stage[0] += ' (STOPPED)'
    log(f'STOPPED by {e} in stage {stage[0]} after {counter[0]} requests: summary written (everything above is real answers; the stage was not finished)')
    write_summary(final=True)
    for tmp in (f'{out}/found/.tmp', f'{out}/cdx/.page'):
        if os.path.exists(tmp):
            os.remove(tmp)
    sys.exit(143 if str(e) == 'SIGTERM' else 130)
