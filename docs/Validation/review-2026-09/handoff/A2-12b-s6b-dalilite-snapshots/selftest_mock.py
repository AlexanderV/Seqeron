#!/usr/bin/env python3
"""Self-test of fetch_dalilite_snapshots.py / in_container.sh / pack_out.sh against a LOCAL MOCK of the Internet Archive (no network; the real archive is not reachable from the audit sandbox).
The mock serves /cdx/search/cdx (matchType prefix / exact, filter=[!]field:regex, collapse, limit, showResumeKey / resumeKey) and /web/<when>id_/<url> (nearest capture, 302 to the full timestamp, HEAD, Content-Length or not,
429 every Nth request, stalls). Fixture releases are tiny tar.gz files built in memory. Run: python3 selftest_mock.py   (prints SELFTEST OK; about 1 - 2 minutes)
Scenarios checked: (1) full run: truncated CDX (limit 5 x 3 pages) with TRUNCATED logged, an archive that only the second CDX page lists, a DB tarball rejected by the basename filter (no GET), a huge file with Content-Length
(HEAD -> TOO-LARGE, no GET), a huge chunked file without length (GET killed at the cap -> TOO-LARGE), duplicate captures (3 captures of one URL, 2 distinct digests -> 2 GETs, the same digest never twice), an archive reachable only through a
link on a CDX-listed page, a licence text page saved under pages/, 429 retried, the run-2 known-404 list honoured, inspection of every kept archive written as it is kept, SUMMARY line; (2) DALI_REPROBE asks the run-2 pairs again;
(3) SIGTERM while a GET stalls: the script exits 143 within seconds with report.txt / summary.json written and no .tmp left; (4) in_container.sh forwards SIGTERM to the script; (5) pack_out.sh excludes out/found, .tmp and big files."""
import gzip, hashlib, io, os, re, signal, subprocess, sys, tarfile, tempfile, threading, time, urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, 'fetch_dalilite_snapshots.py')
H = 'http://ekhidna.biocenter.helsinki.fi/dali_lite/downloads/'
BIG = 6_000_000


def release(version, changes_extra='', licence='ACADEMIC LICENCE AGREEMENT FOR THE PROGRAM DaliLite'):
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode='w:gz') as t:
        def add(name, text):
            b = text.encode()
            ti = tarfile.TarInfo(f'DaliLite_{version}/{name}')
            ti.size = len(b)
            t.addfile(ti, io.BytesIO(b))
        add('CHANGES', f'Changes in DaliLite v.{version} (Feb 2010)\n{changes_extra}\n')
        add('LICENCE', licence + '\n')
        add('README', 'DaliLite fixture\n')
        add('DaliLite', '#!/usr/bin/perl\nprint "Running";\n# DaliLite -pairwise <pdbfile1> <pdbfile2> <jobid>\nprint OUT "<h3> Query=$cd1 Sbjct=$cd2 </h3>";\n')
    return buf.getvalue()


A, B, C, D = release('3.3', 'a'), release('3.3', 'b - a different capture'), release('2.4.1'), release('2.2')


def data(scenario):
    """url -> list of captures dict(ts, status, mime, body | size, length, hide_cdx, stall)."""
    d = {}

    def cap(url, ts, body=None, status=200, mime='application/x-gzip', size=0, length=True, hide=False, stall=False):
        d.setdefault(url, []).append(dict(ts=ts, status=status, mime=mime, body=body, size=size, length=length, hide=hide, stall=stall))
    cap('http://www.ebi.ac.uk/~holm/DaliLite/', '20021006014055', b'<HTML><A HREF="dalilicence">licenced</A> <A HREF="DaliLite.pdf">pdf</A> DaliLite -align code1 code2</HTML>', mime='text/html')
    cap('http://www.ebi.ac.uk/~holm/DaliLite/dalilicence', '20010430025854', b'ACADEMIC LICENCE for DaliLite 2001 text', mime='text/plain', hide=True)
    cap(H + 'download.html', '20170709182157', b'<html>DaliLite_2.4.1 <a href="DaliLite_2.4.1.tar.gz">old</a> <a href="v3/DAT.tar.gz">db</a> <a href="v3/DaliLite.tar.gz">3.3</a></html>', mime='text/html')
    cap(H + 'DaliLite_2.4.1.tar.gz', '20050301000000', C, hide=True)                      # reachable only through the link on download.html
    cap(H + 'v3/DaliLite.tar.gz', '20140101000000', A)
    cap(H + 'v3/DaliLite.tar.gz', '20170709235707', A)                                    # same digest as the first capture
    cap(H + 'v3/DaliLite.tar.gz', '20200101000000', B)
    cap(H + 'v3/DAT.tar.gz', '20170709235754', size=BIG)                                  # database file: rejected by the basename filter, must never be fetched
    cap(H + 'v2/DaliLite_2.4.5.tar.gz', '20070101000000', size=BIG, length=True)           # huge with Content-Length: HEAD -> TOO-LARGE
    cap(H + 'v4/DaliLite_4.0.tar.gz', '20180101000000', size=BIG, length=False)           # huge without length: the poll kills the GET
    cap(H + 'zzold/DaliLite_2.2.tar.gz', '20040101000000', D)                             # listed on page 2 of the CDX only
    for i in range(40):                                                                   # junk under .../dali/: more than 3 pages of 5 -> TRUNCATED
        cap(f'http://ekhidna.biocenter.helsinki.fi/dali/x{i:03d}.html', '20150101000000', b'x', mime='text/html')
    if scenario == 'stall':
        cap(H + 'v5/DaliLite_5.0.tar.gz', '20190101000000', size=1_000_000, stall=True)
    return d


def norm(u):
    return u.replace(':80/', '/').rstrip('/')


class Mock(BaseHTTPRequestHandler):
    server_version = 'MockArchive'

    def log_message(self, *a):
        pass

    def reply(self, code, body=b'', ctype='text/html', length=True, extra=()):
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        if length:
            self.send_header('Content-Length', str(len(body)))
        for k, v in extra:
            self.send_header(k, v)
        self.end_headers()
        if self.command != 'HEAD':
            try:
                self.wfile.write(body)
            except OSError:
                pass

    def do_HEAD(self):
        self.do_GET()

    def do_GET(self):
        srv = self.server
        with srv.lock:
            srv.n += 1
            srv.hits.append((self.command, self.path))
            n = srv.n
        if srv.delay:
            time.sleep(srv.delay)
        if n % srv.every429 == 0:
            return self.reply(429, b'slow down', extra=[('Retry-After', '0')])
        u = urllib.parse.urlparse(self.path)
        q = urllib.parse.parse_qs(u.query)
        if u.path == '/cdx/search/cdx':
            if srv.refuse_filters and 'filter' in q:
                return self.reply(400, b'bad filter')
            return self.cdx(q)
        m = re.match(r'^/web/(\d+)id_/(.*)$', u.path + ('?' + u.query if u.query else ''))
        if not m:
            return self.reply(404, b'no')
        when, url = m.group(1), urllib.parse.unquote(m.group(2))
        caps = srv.data.get(next((k for k in srv.data if norm(k) == norm(url)), None), [])
        if not caps:
            return self.reply(404, b'<html>nothing archived</html>')
        target = int((when + '0101000000')[:14]) if len(when) < 14 else int(when)
        c = min(caps, key=lambda c: abs(int(c['ts']) - target))
        if when != c['ts']:
            return self.reply(302, b'', extra=[('Location', f'/web/{c["ts"]}id_/{url}')])
        if c['status'] != 200:
            return self.reply(c['status'], b'<html>x</html>')
        if c['body'] is not None:
            return self.reply(200, c['body'], c['mime'])
        size = c['size']
        self.send_response(200)
        self.send_header('Content-Type', c['mime'])
        if c['length']:
            self.send_header('Content-Length', str(size))
        self.end_headers()
        if self.command == 'HEAD':
            return
        sent = 0
        try:
            while sent < size:
                n_ = min(65536, size - sent)
                self.wfile.write((b'\x1f\x8b' + b'\0' * (n_ - 2)) if sent == 0 else b'\0' * n_)
                sent += n_
                if c['stall'] and sent >= 100_000:
                    time.sleep(60)
                    return
        except OSError:
            with srv.lock:
                srv.aborted.append(url)

    def cdx(self, q):
        srv = self.server
        url = q['url'][0]
        prefix = q.get('matchType', ['exact'])[0] == 'prefix'
        rows = []
        for orig, caps in srv.data.items():
            for c in caps:
                if c['hide']:
                    continue
                if (orig.startswith(url) if prefix else norm(orig) == norm(url)):
                    digest = hashlib.sha1((c['body'] or str(c['size']).encode()) + (b'1' if c['stall'] else b'')).hexdigest()[:12].upper() if c['body'] is not None or c['size'] else 'D'
                    if orig.endswith('v3/DAT.tar.gz') or c['size']:
                        digest = hashlib.sha1(orig.encode()).hexdigest()[:12].upper()
                    rows.append((orig, c['ts'], str(c['status']), c['mime'], digest))
        rows.sort()
        for f in q.get('filter', []):
            neg = f.startswith('!')
            field, regex = f.lstrip('!').split(':', 1)
            idx = {'original': 0, 'statuscode': 2, 'mimetype': 3}[field]
            rows = [r for r in rows if bool(re.fullmatch(regex, r[idx])) != neg]
        if q.get('collapse', [''])[0] == 'urlkey':
            seen, res = set(), []
            for r in rows:
                if r[0] not in seen:
                    seen.add(r[0])
                    res.append(r)
            rows = res
        if q.get('resumeKey'):
            k = q['resumeKey'][0]
            rows = [r for r in rows if f'{r[0]} {r[1]}' > k]
        limit = int(q.get('limit', ['5000'])[0])
        page, more = rows[:limit], len(rows) > limit
        text = ''.join(f'{r[1]} {r[2]} {r[3]} {r[4]} {r[0]}\n' for r in page)
        if more and q.get('showResumeKey'):
            text += f'\n{page[-1][0]} {page[-1][1]}\n'
        self.reply(200, text.encode(), 'text/plain')


def serve(scenario, every429=40, delay=0.0, refuse_filters=False):
    srv = ThreadingHTTPServer(('127.0.0.1', 0), Mock)
    srv.daemon_threads = True
    srv.data, srv.lock, srv.n, srv.hits, srv.aborted, srv.delay, srv.every429, srv.refuse_filters = data(scenario), threading.Lock(), 0, [], [], delay, every429, refuse_filters
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, f'http://127.0.0.1:{srv.server_address[1]}'


def env_for(host, **kw):
    e = dict(os.environ, DALI_ARCHIVE_HOST=host, DALI_DELAY='0', DALI_BACKOFF_SCALE='0.01', DALI_MAX_BYTES='2000000', DALI_CDX_LIMIT='5', DALI_CDX_PAGES='3', DALI_GET_TIME='30')
    e.update(kw)
    return e


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg)
    if not cond:
        FAILS.append(msg)


FAILS = []


def main():
    tmp = tempfile.mkdtemp(prefix='dalisim_')
    # ---- (1) full run
    print('(1) full run against the mock')
    srv, host = serve('main')
    work = os.path.join(tmp, 'w1')
    os.makedirs(work)
    known = os.path.join(work, 'run2_known_404.txt')
    open(known, 'w').write('http://www.bioinformatics.org/dali/DaliLite_1.0.tar.gz\n' + H + 'v3/DAT.tar.gz\n')      # the first: a run-2 pair; the second is also a CDX candidate? (rejected: no)
    out = os.path.join(work, 'out')
    t = time.time()
    r = subprocess.run([sys.executable, '-I', SCRIPT, out], env=env_for(host, DALI_KNOWN_404=known), capture_output=True, text=True, timeout=900)
    print(f'     exit {r.returncode} after {time.time() - t:.0f} s, {srv.n} mock requests')
    log = open(f'{out}/requests.log').read()
    rep = open(f'{out}/report.txt').read()
    paths = [p for (_, p) in srv.hits]
    gets = [p for (m_, p) in srv.hits if m_ == 'GET' and p.startswith('/web/')]
    check(r.returncode == 0, 'exit code 0')
    check(rep.startswith('SUMMARY FINAL') and 'archives_kept=4' in rep.split('\n')[0], 'report.txt starts with SUMMARY FINAL ... archives_kept=4: ' + rep.split('\n')[0][:150])
    check(len(os.listdir(f'{out}/found')) == 4 and not os.path.exists(f'{out}/found/.tmp'), '4 distinct archives in found/, no .tmp left')
    check('TRUNCATED' in log and 'cdx_truncated=1' in rep.split('\n')[0], 'TRUNCATED logged for the 40-row base (5 per page x 3 pages) and counted in the SUMMARY')
    check(any('zzold/DaliLite_2.2' in l and l.startswith('KEPT') for l in log.split('\n')), 'the archive that only CDX page 2 lists was kept (resumeKey pagination works)')
    check(re.search(r'REJECT .*v3/DAT\.tar\.gz', log) is not None and not any('DAT.tar.gz' in p for p in gets), 'DAT.tar.gz rejected by the basename filter and never fetched')
    check(not any('DaliLite_2.4.5' in p and p.startswith('/web/') and m_ == 'GET' for (m_, p) in srv.hits) and 'HEAD content-length' in log, 'huge file with Content-Length: HEAD -> TOO-LARGE, no GET')
    check(any('DaliLite_4.0' in p for p in gets) and re.search(r'RESULT TOO-LARGE .*DaliLite_4\.0\.tar\.gz rc=63', log) is not None and any('DaliLite_4.0' in a for a in srv.aborted), 'huge chunked file: GET started, killed at the cap (rc 63), TOO-LARGE, not retried (1 GET): ' + str(sum('DaliLite_4.0' in p for p in gets)))
    check(sum('DaliLite_4.0' in p for p in gets) == 1, 'the TOO-LARGE GET was not retried')
    v3 = [p for p in gets if 'v3/DaliLite.tar.gz' in p]
    check(len(v3) == 2 and not any('20170709235707' in p for p in v3), 'v3/DaliLite.tar.gz: 3 captures, 2 distinct digests -> exactly 2 GETs, the duplicate capture never fetched: ' + str(v3))
    check(any(l.startswith('KEPT\tnew') and 'DaliLite_2.4.1' in l for l in log.split('\n')), 'archive reachable only through a link on the CDX-listed page download.html was kept')
    pages = os.listdir(f'{out}/pages')
    check(any('dalilicence' in p for p in pages) or any('licence' in open(f'{out}/pages/{p}').read().lower() for p in pages), 'the licence page was saved as text under pages/: ' + str(pages))
    check('BACKOFF' in log and 'RESULT FAILED' not in log, '429 answers retried with back-off, no request FAILED')
    check(not any('bioinformatics.org/dali/DaliLite_1.0.tar.gz' in p for p in paths), 'a (base, name) pair in the run-2 known-404 list was not asked')
    check(not any('DaliLite_1.0.tar.gz' in p and 'www.ebi.ac.uk/dali/DaliLite/DaliLite_1.0' in p for p in paths), 'old bases x old names (run 2) not asked again')
    ins = sorted(os.listdir(f'{out}/inspect'))
    check(len(ins) == 4 and all(os.path.exists(f'{out}/inspect/{d}/members.txt') for d in ins), 'inspection written for every kept archive: ' + str(ins))
    check(rep.count('Query_Sbjct_lines_in_text_members=1') >= 4 and 'licence_file=' in rep and 'CHANGES_versions=3.3' in rep, 'report line per archive: licence file, CHANGES versions, Query/Sbjct line count')
    check(open(f'{out}/found_index.tsv').read().count('\n') == 5, 'found_index.tsv: 4 rows + header')
    check('PROGRESS' in r.stdout and 'PLAN run 3' in r.stdout, 'progress and plan lines on the console (stdout)')
    srv.shutdown()
    # pack_out.sh
    print('(5) pack_out.sh')
    open(f'{out}/found/.tmp', 'wb').write(b'x' * 5_000_000)
    open(f'{out}/big.bin', 'wb').write(b'x' * 4_000_000)
    zf = os.path.join(work, 'o.zip')
    subprocess.run(['bash', os.path.join(HERE, 'pack_out.sh'), out, zf], check=True)
    names = subprocess.run(['unzip', '-Z1', zf], capture_output=True, text=True).stdout.split('\n')
    check(any(n.endswith('requests.log') for n in names) and any(n.endswith('report.txt') for n in names), 'zip holds requests.log and report.txt')
    check(not any('/found/' in n for n in names) and not any('.tmp' in n for n in names) and not any(n.endswith('big.bin') for n in names), 'zip excludes out/found/*, .tmp and the 4 MB file; size %d bytes' % os.path.getsize(zf))
    # ---- (2) DALI_REPROBE
    print('(2) DALI_REPROBE with DALI_QUICK')
    srv, host = serve('main', every429=10 ** 9)
    work2 = os.path.join(tmp, 'w2')
    os.makedirs(work2)
    open(os.path.join(work2, 'run2_known_404.txt'), 'w').write('x\n')
    r = subprocess.run([sys.executable, '-I', SCRIPT, os.path.join(work2, 'out')], env=env_for(host, DALI_QUICK='1', DALI_REPROBE='1'), capture_output=True, text=True, timeout=300)
    check(r.returncode == 0 and any('/web/2008id_/http://www.ebi.ac.uk/~holm/DaliLite/DaliLite_1.0.tar.gz' in p for (_, p) in srv.hits), 'with DALI_REPROBE the run-2 pairs are asked again')
    srv.shutdown()
    # ---- (6) CDX filters refused
    print('(6) CDX filter syntax refused (HTTP 400): fallback to the plain paginated query')
    srv, host = serve('main', every429=10 ** 9, refuse_filters=True)
    work6 = os.path.join(tmp, 'w6')
    os.makedirs(work6)
    r = subprocess.run([sys.executable, '-I', SCRIPT, os.path.join(work6, 'out')], env=env_for(host, DALI_QUICK='1'), capture_output=True, text=True, timeout=300)
    check(r.returncode == 0 and 'fallback: no filters' in r.stdout and any(l.startswith('KEPT\tnew') for l in r.stdout.split('\n')), 'filters refused -> plain query used, archives still found')
    srv.shutdown()
    # ---- (3) SIGTERM while a GET stalls
    print('(3) SIGTERM while a GET stalls')
    srv, host = serve('stall', every429=10 ** 9)
    work3 = os.path.join(tmp, 'w3')
    os.makedirs(work3)
    out3 = os.path.join(work3, 'out')
    p = subprocess.Popen([sys.executable, '-I', SCRIPT, out3], env=env_for(host, DALI_GET_TIME='120'), stdout=open(os.path.join(work3, 'console.txt'), 'w'), stderr=subprocess.STDOUT, text=True)
    deadline = time.time() + 120
    while time.time() < deadline and not any('DaliLite_5.0' in pp and m_ == 'GET' and pp.startswith('/web/') for (m_, pp) in srv.hits):
        time.sleep(0.5)
    time.sleep(1.5)
    t = time.time()
    p.send_signal(signal.SIGTERM)
    try:
        rc = p.wait(timeout=20)
    except subprocess.TimeoutExpired:
        p.kill()
        rc = 'HUNG'
    dt = time.time() - t
    check(rc == 143 and dt < 10, f'SIGTERM during a stalled GET: exit {rc} after {dt:.1f} s')
    rep3 = open(f'{out3}/report.txt').read().split('\n')[0] if os.path.exists(f'{out3}/report.txt') else ''
    check(rep3.startswith('SUMMARY FINAL') and 'STOPPED' in rep3 and os.path.exists(f'{out3}/summary.json'), 'report.txt SUMMARY written on SIGTERM: ' + rep3[:110])
    check(not os.path.exists(f'{out3}/found/.tmp') and 'STOPPED by SIGTERM' in open(f'{out3}/requests.log').read(), 'no .tmp left, STOPPED line logged')
    srv.shutdown()
    # ---- (4) in_container.sh forwards SIGTERM
    print('(4) in_container.sh: SIGTERM is forwarded to the script')
    srv, host = serve('stall', every429=10 ** 9)
    work4 = os.path.join(tmp, 'w4')
    os.makedirs(work4)
    sh = open(os.path.join(HERE, 'in_container.sh')).read()
    sh = re.sub(r'for u in .*?done\n', '', sh, flags=re.S).replace('/work', work4).replace('/opt/fetch_dalilite_snapshots.py', SCRIPT)
    open(os.path.join(work4, 'run.sh'), 'w').write(sh)
    p = subprocess.Popen(['bash', os.path.join(work4, 'run.sh')], env=env_for(host, DALI_GET_TIME='120'), stdout=open(os.path.join(work4, 'console.txt'), 'w'), stderr=subprocess.STDOUT, text=True, start_new_session=True)
    deadline = time.time() + 120
    while time.time() < deadline and not any('DaliLite_5.0' in pp and m_ == 'GET' and pp.startswith('/web/') for (m_, pp) in srv.hits):
        time.sleep(0.5)
    time.sleep(1.5)
    t = time.time()
    p.send_signal(signal.SIGTERM)
    try:
        rc = p.wait(timeout=30)
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        rc = 'HUNG'
    check(rc == 0 and time.time() - t < 15, f'in_container.sh ended with {rc} after {time.time() - t:.1f} s')
    check(os.path.exists(f'{work4}/out/report.txt') and 'STOPPED' in open(f'{work4}/out/report.txt').read().split('\n')[0] and 'STOPPED by SIGTERM' in open(f'{work4}/out/fetch.log').read(), 'out/report.txt and out/fetch.log (tee) written after the signal')
    srv.shutdown()
    print('SELFTEST OK' if not FAILS else f'SELFTEST FAILED ({len(FAILS)}): ' + '; '.join(FAILS))
    sys.exit(1 if FAILS else 0)


if __name__ == '__main__':
    main()
