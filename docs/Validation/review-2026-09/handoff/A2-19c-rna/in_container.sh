#!/usr/bin/env bash
# Runs inside the image: for every set of inputs/sets.tsv `mafft <options> in.fa` twice (the extended build at /opt/mafft-ext when it exists, the Debian build otherwise: the version line says which); every run is recorded
# (return code, standard output, standard error, files, 'deterministic' = both runs equal). Steps are logged to out/steps.log; nothing aborts. RUN_TIMEOUT (default 300) seconds per run. RNA_QUICK=1: 6 sets.
set -uo pipefail
WORK="${WORK:-/work}"
cd "$WORK"
rm -rf out; mkdir -p out
step() { echo "== STEP $* [$(date -u +%H:%M:%S)]" | tee -a out/steps.log; }
step "versions"
{
  echo "date: $(date -u)"; echo "uname: $(uname -m)"
  echo "debian mafft: $(mafft --version 2>&1 | head -1)"; echo "extended mafft: $(/opt/mafft-ext/bin/mafft --version 2>&1 | head -1)"
  echo "RNAalifold: $(command -v RNAalifold) $(RNAalifold --version 2>&1 | head -1)"; echo "mafft on PATH: $(command -v mafft)"; echo "curl: $(command -v curl)"
  echo "contrafold: $(ls -l /opt/mafft-ext/libexec/mafft/contrafold 2>&1)"; echo "dafs: $(ls -l /opt/mafft-ext/libexec/mafft/dafs 2>&1)"
  echo "--- programs of the extended build"; ls -l /opt/mafft-ext/libexec/mafft 2>&1
  echo "--- sha256 of the programs"; sha256sum /opt/mafft-ext/libexec/mafft/* 2>&1 | head -60
} > out/versions.txt 2>&1
cp /opt/build_ext.log /opt/build_vienna.log /opt/build_rnaprogs.log out/ 2>/dev/null
step "the sets"
python3 -I - <<'PY'
import json, os, re, subprocess, tempfile, time, hashlib

WORK = os.environ.get('WORK', '/work')
QUICK = os.environ.get('RNA_QUICK') == '1'
TIMEOUT = int(os.environ.get('RUN_TIMEOUT', '300'))
MAFFT = '/opt/mafft-ext/bin/mafft' if os.path.exists('/opt/mafft-ext/bin/mafft') else 'mafft'


def run_once(td, argv):
    for f in os.listdir(td):
        if f != 'in.fa' and os.path.isfile(os.path.join(td, f)):
            os.remove(os.path.join(td, f))
    env = dict(os.environ, TMPDIR=td)
    t0 = time.time()
    try:
        p = subprocess.run(argv, capture_output=True, cwd=td, env=env, timeout=TIMEOUT)
        rc, out, err = p.returncode, p.stdout.decode('latin-1'), p.stderr.decode('latin-1')
    except subprocess.TimeoutExpired as e:
        rc, out, err = 'timeout', (e.stdout or b'').decode('latin-1'), (e.stderr or b'').decode('latin-1')
    files = {f: open(os.path.join(td, f), 'rb').read().decode('latin-1') for f in sorted(os.listdir(td)) if f != 'in.fa' and os.path.isfile(os.path.join(td, f))}
    return {'rc': rc, 'stdout': out, 'stderr': err, 'files': files, 'seconds': round(time.time() - t0, 1)}


sets = [l.rstrip('\n').split('\t') for l in open(os.path.join(WORK, 'inputs/sets.tsv')) if l.strip()]
if QUICK:
    sets = sets[::24][:6]
os.makedirs(os.path.join(WORK, 'out/sets'), exist_ok=True)
index = []
for k, opts, seqs in sets:
    seqs = seqs.split('|')
    text = ''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs))
    argv = [MAFFT] + opts.split() + ['in.fa']
    with tempfile.TemporaryDirectory() as td:
        open(os.path.join(td, 'in.fa'), 'w').write(text)
        a = run_once(td, argv)
        b = run_once(td, argv)
    def norm(t):
        # the process number in the shell's "mafft: line N: PID Aborted / Segmentation fault" and the random names of the temporary directories are not part of the program's behaviour
        t = re.sub(r'(mafft: line \d+:)\s+\d+ (Aborted|Segmentation fault)', r'\1 P \2', t)
        return re.sub(r'mafft\.[A-Za-z0-9]{6,12}', 'mafft.X', t)
    det = all((norm(a[x]) if x == 'stderr' else a[x]) == (norm(b[x]) if x == 'stderr' else b[x]) for x in ('rc', 'stdout', 'stderr', 'files'))
    rec = {'id': k, 'args': argv[1:], 'input': text, 'deterministic': det, 'run1': a, 'run2': b if not det else None}
    json.dump(rec, open(os.path.join(WORK, 'out/sets', k + '.json'), 'w'))
    index.append((k, opts, a['rc'], det, a['seconds']))
    print(k, opts, 'rc', a['rc'], 'det' if det else 'NON-DETERMINISTIC', a['seconds'], 's', flush=True)
json.dump(index, open(os.path.join(WORK, 'out/index.json'), 'w'))
PY
step "probes (run 4): CONTRAfold and MXSCARNA run directly on inputs/probe_seqs.tsv; their outputs are what the port has to equal exactly"
python3 -I - <<'PY'
import json, os, subprocess, tempfile

WORK = os.environ.get('WORK', '/work')
LIB = '/opt/mafft-ext/libexec/mafft'
TIMEOUT = int(os.environ.get('RUN_TIMEOUT', '300'))


def run(argv, td):
    try:
        p = subprocess.run(argv, capture_output=True, cwd=td, timeout=TIMEOUT, env=dict(os.environ, PATH=LIB + ':' + os.environ.get('PATH', '')))
        return {'rc': p.returncode, 'stdout': p.stdout.decode('latin-1'), 'stderr': p.stderr.decode('latin-1')}
    except subprocess.TimeoutExpired:
        return {'rc': 'timeout', 'stdout': '', 'stderr': ''}


records = []
probes = [l.rstrip('\n').split('\t') for l in open(os.path.join(WORK, 'inputs/probe_seqs.tsv')) if l.strip()]
if os.environ.get('RNA_QUICK') == '1':
    probes = probes[::12]
for k, seq in probes:
    rec = {'id': k, 'seq': seq}
    with tempfile.TemporaryDirectory() as td:
        open(os.path.join(td, 'in.fa'), 'w').write('>in\n' + seq + '\n')
        for name, argv, outfile in (('contrafold', [LIB + '/contrafold', 'predict', 'in.fa', '--posteriors', '0.01', 'out.txt'], 'out.txt'),
                                    ('contrafold_gamma', [LIB + '/contrafold', 'predict', 'in.fa'], None),
                                    ('mxscarnamod', [LIB + '/mxscarnamod', '-m', '-writebpp', 'in.fa'], None)):
            a = run(argv, td)
            b = run(argv, td)
            if outfile and os.path.exists(os.path.join(td, outfile)):
                a['file'] = open(os.path.join(td, outfile), encoding='latin-1').read()
                os.remove(os.path.join(td, outfile))
            a['repeat_equal'] = (a['rc'], a['stdout']) == (b['rc'], b['stdout'])
            rec[name] = a
    records.append(rec)
    print(k, len(seq), rec['contrafold']['rc'], rec['mxscarnamod']['rc'], flush=True)
json.dump(records, open(os.path.join(WORK, 'out/probes.json'), 'w'))
PY
step "done"
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" out 2>/dev/null
exit 0
