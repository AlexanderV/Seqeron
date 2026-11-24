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
import gzip, json, os, shutil, subprocess, sys, tempfile, time

QUICK = os.environ.get('LAST_QUICK') == '1'
TIMEOUT = int(os.environ.get('RUN_TIMEOUT', '300'))
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


sets = [l.rstrip('\n').split('\t') for l in open(os.path.join(WORK, 'inputs/sets.tsv')) if l.strip()]
if QUICK:
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
        b = run_once(td, argv, lastlog)
        strip = lambda r: {x: r[x] for x in r if x != 'seconds'}
        # the run-to-run difference of the time stamps in the names is not part of the comparison; the calls are compared by their content
        same = strip(a) == strip(b)
        rec = {'id': k, 'mode': mode, 'options': opts, 'in.fa': infa, 'new.fa': open(os.path.join(td, 'new.fa')).read(), 'argv': argv, 'deterministic': same, 'run': a, 'run2': None if same else b}
    with gzip.open(os.path.join(WORK, 'out/sets/%s.json.gz' % k), 'wt', encoding='utf-8') as f:
        json.dump(rec, f)
    index.append('\t'.join([k, mode, opts, str(a['rc']), 'det' if same else 'NONDET', str(len(a['calls'])), str(a['seconds'])]))
    if not same:
        nondet.append(k)
    print(k, mode, opts, '->', a['rc'], 'calls=%d' % len(a['calls']), 'det' if same else 'NONDET', flush=True)
open(os.path.join(WORK, 'out/index.tsv'), 'w').write('\n'.join(index) + '\n')
print('sets', len(index), 'non-deterministic', len(nondet), nondet[:10])
PY
step "done"
[ -n "${HOST_UID:-}" ] && chown -R "$HOST_UID:${HOST_GID:-$HOST_UID}" out || true
exit 0
