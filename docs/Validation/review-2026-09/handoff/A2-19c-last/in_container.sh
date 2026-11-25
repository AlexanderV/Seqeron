#!/usr/bin/env bash
# Runs inside the image: for every set of inputs/sets.tsv `mafft <options> in.fa` (WITHOUT --quiet unless the options say it) twice with the logging wrappers of lastdb / lastal; only runs whose two repetitions are equal (return code, standard output,
# standard error, files, every recorded call of LAST) are written. Steps are logged to out/steps.log; nothing aborts.
#   LAST_QUICK=1       only 6 sets (a check of the package, ~ 1 minute)
#   RUN_TIMEOUT=300    seconds per mafft run (a run that exceeds it is recorded as "timeout")
#   WORK, MAFFT_PREFIX directories (defaults /work and /usr/lib/mafft/lib/mafft)
set -uo pipefail
WORK="${WORK:-/work}"
PREFIX="${MAFFT_PREFIX:-/usr/lib/mafft/lib/mafft}"
cd "$WORK"
rm -rf out; mkdir -p out/sets
step() { echo "== STEP $* [$(date -u +%H:%M:%S)]" | tee -a out/steps.log; }
step "versions"
{
  echo "date: $(date -u)"; echo "uname: $(uname -a)"; mafft --version 2>&1 | head -1
  echo "mafft prefix: $PREFIX"; ls -l "$PREFIX" 2>&1 | grep -i -E "last|pairlocalalign" 
  echo "LAST selected: $(readlink -f /opt/last 2>&1)"; "${LAST_REAL_DIR:-/opt/last/bin}/lastal" --version 2>&1 | head -1; sha256sum "${LAST_REAL_DIR:-/opt/last/bin}/lastal" "${LAST_REAL_DIR:-/opt/last/bin}/lastdb" 2>&1
  echo "--- candidates"; cat /opt/last-candidates.txt 2>&1
} > out/versions.txt 2>&1
cp /opt/build_last.log out/ 2>/dev/null; cp /opt/last-candidates.txt out/ 2>/dev/null
step "install the logging wrappers into $PREFIX"
for t in lastdb lastal; do cp "$WORK/last-log" "$PREFIX/$t" && chmod +x "$PREFIX/$t"; done
export LAST_REAL_DIR="${LAST_REAL_DIR:-/opt/last/bin}"
step "the sets"
python3 -I - <<'PY'
import gzip, json, os, re, shutil, subprocess, sys, tempfile, time

QUICK = os.environ.get('LAST_QUICK') == '1'
TIMEOUT = int(os.environ.get('RUN_TIMEOUT', '3600'))          # seconds per mafft run (run 1: 300, but lastal of LAST >= 1447 needs ~26 s per call natively)
TWICE_MAX = int(os.environ.get('TWICE_MAX', '900'))            # a set whose first run took longer is run once only (determinism 'unchecked')
ONLY = os.environ.get('LAST_ONLY', '').split()               # self-test: only these ids
WORK = os.environ.get('WORK', '/work')


def fasta(seqs, prefix='s'):
    return ''.join('>%s%d\n%s\n' % (prefix, i, s) for i, s in enumerate(seqs))


def run_once(td, argv, lastlog):
    shutil.rmtree(lastlog, ignore_errors=True)
    os.makedirs(lastlog)
    for f in os.listdir(td):
        if f not in ('in.fa', 'new.fa') and f != 'lastlog':
            p = os.path.join(td, f)
            if os.path.isfile(p):
                os.remove(p)
    env = dict(os.environ, LASTLOG=lastlog, TMPDIR=td)
    t0 = time.time()
    try:
        p = subprocess.run(argv, capture_output=True, cwd=td, env=env, timeout=TIMEOUT)
        rc, out, err = p.returncode, p.stdout.decode('latin-1'), p.stderr.decode('latin-1')
    except subprocess.TimeoutExpired as e:
        rc, out, err = 'timeout', (e.stdout or b'').decode('latin-1'), (e.stderr or b'').decode('latin-1')
    files = {}
    for f in sorted(os.listdir(td)):
        p = os.path.join(td, f)
        if f in ('in.fa', 'new.fa', 'lastlog') or not os.path.isfile(p):
            continue
        files[f] = open(p, 'rb').read().decode('latin-1')
    groups = {}
    for f in sorted(os.listdir(lastlog)):
        parts = f.split('-')
        key = '-'.join(parts[:3])           # time stamp (ns) - process id - random number: the order of the calls
        rest = '-'.join(parts[3:])          # <program>.<kind>
        prog, kind = rest.split('.', 1)
        g = groups.setdefault(key, {'prog': prog, 'items': {}})
        g['items'][kind] = open(os.path.join(lastlog, f), 'rb').read().decode('latin-1')
    calls = [groups[k] for k in sorted(groups, key=lambda x: int(x.split('-')[0]))]
    return {'rc': rc, 'stdout': out, 'stderr': err, 'files': files, 'calls': calls, 'seconds': round(time.time() - t0, 1)}


def norm(r):
    # run 1 declared every set NONDET because the random names of the temporary directories (mafft.XXXXXXXXXX, /tmp/tmpxxxxxxxx) and the process numbers of the crash message differ from run to run: they are not part of the comparison
    r2 = {k: v for k, v in r.items() if k != 'seconds'}
    r2['calls'] = [{'prog': c['prog'], 'items': {x: y for x, y in c['items'].items() if x != 'ms'}} for c in r['calls']]      # 'ms' = the time of the call: not part of the comparison
    t = json.dumps(r2, sort_keys=True)
    t = re.sub(r'mafft\.[A-Za-z0-9]{6,12}', 'mafft.X', t)
    t = re.sub(r'/tmp/tmp[a-z0-9_]{6,10}', '/tmp/tmpX', t)
    t = re.sub(r'(mafft: line \d+:)\s+\d+ (Segmentation fault)', r'\1 P \2', t)
    t = re.sub(r'\d{10,19}-\d+-\d+', 'T-P-R', t)
    return t


sets = [l.rstrip('\n').split('\t') for l in open(os.path.join(WORK, 'inputs/sets.tsv')) if l.strip()]
if ONLY:
    sets = [x for x in sets if x[0] in ONLY]
elif QUICK:
    sets = sets[:4] + [s for s in sets if s[1] == 'add'][:2]
index = []
nondet = []
for k, mode, opts, old, new in sets:
    old = old.split('|')
    new = new.split('|') if new else []
    with tempfile.TemporaryDirectory() as td:
        open(os.path.join(td, 'new.fa'), 'w').write(fasta(new, 'n') if new else '')
        if mode == 'add':
            q = subprocess.run(['mafft', '--quiet', '-'], input=fasta(old, 'a').encode(), capture_output=True)
            infa = q.stdout.decode('latin-1') if q.returncode == 0 else fasta(old, 'a')
        else:
            infa = fasta(old)
        open(os.path.join(td, 'in.fa'), 'w').write(infa)
        argv = ['mafft'] + opts.split() + ['in.fa']
        lastlog = os.path.join(tempfile.gettempdir(), 'lastlog-%s' % k)
        a = run_once(td, argv, lastlog)
        # the run-to-run difference of the time stamps in the names is not part of the comparison; the calls are compared by their content
        if a['seconds'] > TWICE_MAX:
            b, same = None, None          # too long to repeat: determinism not checked
        else:
            b = run_once(td, argv, lastlog)
            same = norm(a) == norm(b)
        rec = {'id': k, 'mode': mode, 'options': opts, 'in.fa': infa, 'new.fa': open(os.path.join(td, 'new.fa')).read(), 'argv': argv, 'deterministic': same, 'run': a, 'run2': None if same is not False else b}
    with gzip.open(os.path.join(WORK, 'out/sets/%s.json.gz' % k), 'wt', encoding='utf-8') as f:
        json.dump(rec, f)
    label = 'unchecked' if same is None else ('det' if same else 'NONDET')
    index.append('\t'.join([k, mode, opts, str(a['rc']), label, str(len(a['calls'])), str(a['seconds'])]))
    if same is False:
        nondet.append(k)
    print(k, mode, opts, '->', a['rc'], 'calls=%d' % len(a['calls']), label, a['seconds'], 's', flush=True)
    open(os.path.join(WORK, 'out/index.tsv'), 'w').write('\n'.join(index) + '\n')          # after every set: an interrupted run keeps what it has
print('sets', len(index), 'non-deterministic', len(nondet), nondet[:10])
PY
step "valgrind proof for the crashes of pairlocalalign (heap overflow: core/pairlocalalign.c putlocalhom_last is called with localhomtable[i]+j for alg 'R', but row i of the table has njob-i entries)"
mkdir -p out/valgrind
python3 -I - <<'PY'
import gzip, json, os, subprocess, tempfile
WORK = os.environ.get('WORK', '/work')
done = 0
for k in sorted(os.listdir(os.path.join(WORK, 'out/sets'))):
    if done >= 3 or not os.path.exists('/usr/bin/valgrind'):
        break
    rec = json.load(gzip.open(os.path.join(WORK, 'out/sets', k), 'rt', encoding='utf-8'))
    if 'Segmentation fault' not in rec['run']['stderr'] or rec['mode'] != 'plain':
        continue
    with tempfile.TemporaryDirectory() as td:
        open(os.path.join(td, 'in.fa'), 'w').write(rec['in.fa'])
        try:
            p = subprocess.run(['valgrind', '-q', '--trace-children=yes', '--trace-children-skip=*/lastal*,*/lastdb*,*/last-log', '--error-limit=no'] + rec['argv'], capture_output=True, cwd=td, timeout=7200)
            open(os.path.join(WORK, 'out/valgrind', rec['id'] + '.log'), 'wb').write(p.stderr)
        except Exception as e:
            open(os.path.join(WORK, 'out/valgrind', rec['id'] + '.log'), 'w').write('valgrind run failed: %r' % e)
    done += 1
PY
step "done"
[ -n "${HOST_UID:-}" ] && chown -R "$HOST_UID:${HOST_GID:-$HOST_UID}" out || true
exit 0
