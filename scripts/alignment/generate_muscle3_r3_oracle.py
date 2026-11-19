#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/muscle3_r3_oracle.tsv (B09 audit A2-1 / A1-10c part a3) from the ROUND 3 sets.

Round 3 = docs/Validation/review-2026-09/handoff/A1-10c-muscle-r3/inputs/sets.tsv (60 sets R330 ... R350 `-matrix blosum62nostar`, R638 ... R645 `-matrix` with nucleotides, R646 ... R689
one-/two-sequence corners and the score flags -le -sp -sv -spn).

Source of the rows (in this order):
  1. --zip PATH, or the default returned/A1-10c-muscle-r3-out.zip when it exists: the user-run Docker package (muscle3 3.8.1551 built from the Debian source);
  2. otherwise the sets are RUN LOCALLY with the installed Debian `muscle3` binary (`apt-get install muscle3`, /usr/bin/muscle3, 3.8.1551-3 = the same upstream source) through the package's own harness
     run_sets.py (MUSCLE3=<binary>); the version banner is asserted to be MUSCLE v3.8.1551. ENV.txt of that run is kept in handoff/A1-10c-muscle-r3/returned/ENV-local.txt.
Output columns (same as muscle3_r2_oracle.tsv): id kind type options data rc names rows tree1 tree2 dist1 dist2 weights usedtree stdout matrixtext
"""
import os, re, shutil, subprocess, sys, tempfile, zipfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
pkg = os.path.join(root, 'docs', 'Validation', 'review-2026-09', 'handoff', 'A1-10c-muscle-r3')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'muscle3_r3_oracle.tsv')
zpath = os.path.join(pkg, 'returned', 'A1-10c-muscle-r3-out.zip')
if '--zip' in sys.argv:
    zpath = sys.argv[sys.argv.index('--zip') + 1]

sets = {}
with open(os.path.join(pkg, 'inputs', 'sets.tsv')) as f:
    next(f)
    for line in f:
        sid, kind, typ, data, opts = line.rstrip('\n').split('\t')
        sets[sid] = (kind, typ, data, opts)
matrices = {'@' + n[:-4]: open(os.path.join(pkg, 'inputs', 'matrices', n)).read() for n in os.listdir(os.path.join(pkg, 'inputs', 'matrices'))}

if os.path.exists(zpath):
    with zipfile.ZipFile(zpath) as z:
        names = z.namelist()
        env = z.read(next(n for n in names if n.endswith('ENV.txt'))).decode()
        lines = z.read(next(n for n in names if n.endswith('muscle3_r3_oracle.tsv'))).decode().splitlines()
else:
    binary = shutil.which('muscle3') or sys.exit('no returned zip and no muscle3 binary: ' + zpath)
    ver = subprocess.run([binary, '-version'], capture_output=True, text=True).stdout
    env = f'binary: {binary}\nuname -m: {os.uname().machine}\nversion: {ver.strip()}\n'
    with tempfile.TemporaryDirectory() as d:
        # run_sets.py writes ./out/muscle3_r3_oracle.tsv relative to the cwd
        subprocess.run([sys.executable, '-I', os.path.join(pkg, 'run_sets.py')], cwd=d, env={**os.environ, 'MUSCLE3': binary}, check=True, capture_output=True)
        lines = open(os.path.join(d, 'out', 'muscle3_r3_oracle.tsv'), encoding='utf-8').read().splitlines()
    with open(os.path.join(pkg, 'returned', 'ENV-local.txt'), 'w') as fh:
        fh.write(env)
assert 'MUSCLE v3.8.1551' in env, env

esc = lambda s: s.replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n')
rows = []
for line in lines[1:]:
    f = line.split('\t')
    assert len(f) == 14, (len(f), line[:80])
    sid, kind, typ, opts = f[0], f[1], f[2], f[3]
    k, t, data, o = sets[sid]
    assert (k, t, o) == (kind, typ, opts), sid
    assert f[4] != 'RUNNER_ERROR', sid
    mtext = ''
    for key, text in matrices.items():
        if re.search('=' + re.escape(key) + '(,|$)', opts):
            mtext = text
    rows.append([sid, kind, typ, opts, data] + f[4:] + [esc(mtext)])
with open(out_path, 'w', newline='\n') as fh:
    fh.write('id\tkind\ttype\toptions\tdata\trc\tnames\trows\ttree1\ttree2\tdist1\tdist2\tweights\tusedtree\tstdout\tmatrixtext\n')
    for r in rows:
        fh.write('\t'.join(r) + '\n')
print(len(rows), 'rows;', sum(1 for r in rows if not r[5].isdigit()), 'failed/no-output')
