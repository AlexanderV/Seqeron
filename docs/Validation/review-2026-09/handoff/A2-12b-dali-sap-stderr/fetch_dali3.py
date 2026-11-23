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
        # the wrapper logs every call (arguments, working directory, return code, stdout, stderr) into $DALI_LOG_DIR/call_*.txt and, in call_*.txt.saved, every file that DaliLite's own `rm` removes (the result files:
        # `rm -f ... fort.* subunits.puu units.puu dali.lock list0 domains.puu dali.default puu.default seq.* ...`: a shim `rm` earlier on the PATH copies each file, at most 30 000 bytes, `cat -v`, before it is removed) and the
        # directory listing / the remaining files at the end of the call: the glue's input (`DaliLite -pairwise pdb1.pdb pdb2.pdb > tmp1`) is the standard output, the files show how it was made
        os.makedirs(out + '/shim', exist_ok=True)
        open(out + '/shim/rm', 'w').write('#!/bin/bash\nif [ -n "$DALI_CALL_LOG" ]; then for a in "$@"; do case "$a" in -*) ;; *) if [ -f "$a" ]; then { echo "=== REMOVED $a ($(stat -c %s "$a") bytes)"; head -c 30000 "$a" | cat -v; echo; } >> "$DALI_CALL_LOG"; fi;; esac; done; fi\nexec /bin/rm "$@"\n')
        os.chmod(out + '/shim/rm', 0o755)
        open(out + '/DaliLite.wrap', 'w').write('#!/bin/bash\nD=${DALI_LOG_DIR:-/tmp/dali_log}\nmkdir -p $D\nL=$D/call_$$_$(date +%s%N).txt\n{ echo "ARGS: $*"; echo "CWD: $PWD"; } > $L\n'
                                               f'export DALI_CALL_LOG=$L.saved; PATH={out}/shim:$PATH; : > $DALI_CALL_LOG\n'
                                               f'{exe} "$@" > /tmp/dali_out_$$.txt 2>/tmp/dali_err_$$.txt; rc=$?\n'
                                               'echo "RC: $rc" >> $L; echo "--- STDOUT" >> $L; cat /tmp/dali_out_$$.txt >> $L; echo "--- STDERR" >> $L; cat /tmp/dali_err_$$.txt >> $L\n'
                                               '{ echo "=== LS"; ls -la 2>&1; for f in *; do case "$f" in pdb1.pdb|pdb2.pdb) ;; *) if [ -f "$f" ]; then echo "=== REMAINING $f ($(stat -c %s "$f") bytes)"; head -c 30000 "$f" | cat -v; echo; fi;; esac; done; } >> $DALI_CALL_LOG 2>&1\n'
                                               'cat /tmp/dali_out_$$.txt; cat /tmp/dali_err_$$.txt >&2; /bin/rm -f /tmp/dali_out_$$.txt /tmp/dali_err_$$.txt; exit $rc\n')
        os.chmod(out + '/DaliLite.wrap', 0o755)
        # Run 5 of the package (s4 returned): DaliLite 3.3 has its installation directory compiled into the program / scripts: every call runs
        # `perl /data/backup/zope_data/dali_server/DaliLite_3.3//Bin/forwarder.pl ...` and stops with "Can't open perl script" (354 calls). The program is made to find itself at the path it expects.
        root = os.path.dirname(os.path.dirname(exe))
        # Run 7 / 8 (s4 with full-backbone structures): DSSP runs, but `puu`, `puutos`, `dalicon` (Fortran programs shipped PRECOMPILED in Bin/ for the old g77 runtime) need `libg2c.so.0`.
        # Run 8: (A) `make -B` with my own FFLAGS failed in util.f (labels / `if ... then` lines cut at column 72: the Makefile's own flags, with the g77 line-length switch, were replaced);
        # (B) Ubuntu's libf2c2 linked as libg2c.so.0 loads, but `puu` then stops with `undefined symbol: f_setarg` (g77's libg2c has f_setarg / f_setsig in its main support, the f2c library has not).
        # Now: (A) rebuild keeping going (`make -k -B`) with the Makefile's flags plus `-ffixed-line-length-none -std=legacy -fallow-argument-mismatch`; (B) when a program still does not load,
        # a small shared object named libg2c.so.0 that provides f_setarg / f_setsig (no-ops: the DaliLite programs read their input from files and standard input) and depends on libf2c.
        # Every step is logged (`ldd -r` shows undefined symbols), so that the log decides the next step.
        bindir = os.path.dirname(exe)
        def missing():
            r = sh(f'cd {bindir} && for f in puu puutos dssp dalicon; do [ -e $f ] && (echo "== $f"; file $f | cut -c1-100; ldd -r $f 2>&1 | grep -i "not found\\|undefined"); done; true')
            return r.stdout
        print('LDD BEFORE', missing())
        if 'not found' in missing() or 'undefined' in missing():
            r = sh(f'cd {bindir} && grep -n "FFLAGS\\|^FC\\|^F77\\|g77\\|gfortran\\|^LIB" Makefile | head -30')
            print('MAKEFILE FLAGS', r.stdout)
            r = sh(f'cd {bindir} && make -k -B FC=gfortran F77=gfortran G77=gfortran FFLAGS="-O2 -std=legacy -fallow-argument-mismatch -ffixed-line-length-none -w" 2>&1 | grep -v "^gfortran\\|^gcc" | tail -60')
            print('MAKE -k -B', r.stdout, r.stderr)
            print('LDD AFTER make', missing())
        if 'not found' in missing() or 'undefined' in missing():
            sh('apt-get update -qq; DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends libf2c2 2>&1 | tail -3')
            open(out + '/g2cshim.c', 'w').write('void f_setarg(int argc, char **argv) { (void)argc; (void)argv; }\nvoid f_setsig(void) { }\n')
            r = sh(f'ls /usr/lib/*/libf2c.so.2; D=$(dirname $(ls /usr/lib/*/libf2c.so.2 | head -1)); rm -f $D/libg2c.so.0; gcc -shared -fPIC -o $D/libg2c.so.0 -Wl,-soname,libg2c.so.0 {out}/g2cshim.c -Wl,--no-as-needed -L$D -l:libf2c.so.2 2>&1; ldconfig; ls -l $D/libg2c.so.0')
            print('G2C SHIM', r.stdout, r.stderr)
            print('LDD AFTER shim', missing())
        open(out + '/BUILD_OK', 'w').write(exe + '\n')
        print('BUILD_OK', exe)
        break
else:
    open(out + '/BUILD_FAILED', 'w').write('no DaliLite v3 executable\n')
    print('BUILD_FAILED')
