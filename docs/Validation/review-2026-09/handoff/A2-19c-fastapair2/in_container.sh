#!/usr/bin/env bash
# Runs inside the image: for every set of inputs/sets.tsv (1) `mafft <options> in.fa` WITHOUT --quiet with FASTA_4_MAFFT = a logging wrapper (rows, names, return code, standard error, tree / distance files, every raw fasta34 call: arguments, query, library, raw output),
# (2) dndfast7 [-A] run directly on the same sequences in an empty directory: its hat2, hat3 and standard error.
set -uo pipefail
cd /work
rm -rf out/sets out/results.tsv out/progress.log; mkdir -p out/sets
{ mafft --version 2>&1 | head -1; command -v fasta34; sha256sum "$(command -v fasta34)"; uname -m; tail -3 /opt/fasta34.build.log; ls /opt/FASTA34_BUILD_FAILED 2>&1; ls -l /usr/lib/mafft/lib/mafft/dndfast7 2>&1; } > out/versions.txt 2>&1
export FASTA_4_MAFFT=/work/fasta34-log
python3 -I - <<'PY'
import glob, gzip, hashlib, io, json, os, shutil, subprocess, tarfile, tempfile


TIMEOUT = int(os.environ.get('SET_TIMEOUT', '1200'))     # seconds per program run: a hang is recorded (rc 124), it never blocks the run
DIRECT_FROM = int(os.environ.get('DIRECT_FROM', '1'))    # run 3 sets this to 254: the direct dndfast7 runs of sets 1 - 253 are in the archive of run 2
SETS_FROM = int(os.environ.get('SETS_FROM', '1'))


def run_logged(argv, cwd, stdin_path=None):
    d = os.path.join(cwd, 'fastalog')
    shutil.rmtree(d, ignore_errors=True)
    env = dict(os.environ, FASTALOG=d)
    stdin = open(stdin_path, 'rb') if stdin_path else None
    try:
        p = subprocess.run(argv, capture_output=True, cwd=cwd, env=env, stdin=stdin, timeout=TIMEOUT)
    except subprocess.TimeoutExpired as e:
        p = subprocess.CompletedProcess(argv, 124, e.stdout or b'', (e.stderr or b'') + b'\n@@TIMEOUT after %d s@@\n' % TIMEOUT)
    return p, d


def pack_dir(d):
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode='w:gz') as t:
        for f in sorted(os.listdir(d)) if os.path.isdir(d) else []:
            t.add(os.path.join(d, f), arcname=f)
    return buf.getvalue()


# results.tsv is written row by row (run 2 lost it: the run was stopped after set 253 and the file was only written at the end)
results = open('out/results.tsv', 'a')
progress = open('out/progress.log', 'a')
for line in open('inputs/sets.tsv'):
    cols = line.rstrip('\n').split('\t')
    k, opt, seqs = cols[:3]
    if int(k) < SETS_FROM:
        continue
    print('set', k, opt, flush=True)
    progress.write('set %s %s\n' % (k, opt)); progress.flush()
    extra = json.loads(cols[3]) if len(cols) > 3 else {}      # files that the options name (new.fa, seed.fa, table.txt ...)
    seqs = seqs.split('|')
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, 'in.fa')
        open(fa, 'w').write(extra.pop('in.fa', None) or ''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs)))
        for fname, text in extra.items():
            open(os.path.join(td, fname), 'w').write(text)
        p, d = run_logged(['mafft'] + opt.split() + [fa], td)
        tree = open(fa + '.tree', 'rb').read().decode('latin-1') if os.path.exists(fa + '.tree') else ''
        hat2 = open(fa + '.hat2', 'rb').read().decode('latin-1') if os.path.exists(fa + '.hat2') else ''
        calls = pack_dir(d)
        # (2) dndfast7 directly
        wd = os.path.join(td, 'direct'); os.makedirs(wd)
        sw = ['-A'] if '--fastswpair' in opt else []
        exe = next((x for x in ('/usr/lib/mafft/lib/mafft/dndfast7', '/usr/lib/mafft/libexec/mafft/dndfast7', '/usr/libexec/mafft/dndfast7') if os.path.exists(x)), 'dndfast7')
        if int(k) >= DIRECT_FROM:
            q, d2 = run_logged([exe] + sw, wd, fa)
        else:
            q = subprocess.CompletedProcess([], 0, b'', b'')          # not run again: the direct runs of these sets are in the archive of run 2
        h2 = open(os.path.join(wd, 'hat2'), 'rb').read().decode('latin-1') if os.path.exists(os.path.join(wd, 'hat2')) else ''
        h3 = open(os.path.join(wd, 'hat3'), 'rb').read().decode('latin-1') if os.path.exists(os.path.join(wd, 'hat3')) else ''
    out, names = [], []
    for l in p.stdout.decode('latin-1').splitlines():
        if l.startswith('>'):
            names.append(l[1:]); out.append('')
        elif l.strip():
            out[-1] += l.strip()
    with gzip.open('out/sets/%s.stderr.gz' % k, 'wb') as f:
        f.write(p.stderr)
    with open('out/sets/%s.fastalog.tar.gz' % k, 'wb') as f:
        f.write(calls)
    with gzip.open('out/sets/%s.direct.json.gz' % k, 'wt', encoding='utf-8') as f:
        json.dump({'rc': q.returncode, 'stderr': q.stderr.decode('latin-1'), 'hat2': h2, 'hat3': h3}, f)
    esc = lambda s: s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')
    results.write('\t'.join([k, opt, '|'.join(seqs), str(p.returncode), '|'.join(out), '|'.join(names), hashlib.sha256(p.stderr).hexdigest(), esc(tree), esc(hat2)]) + '\n')
    results.flush()
print('done')
PY
[ -n "${HOST_UID:-}" ] && chown -R "$HOST_UID:${HOST_GID:-$HOST_UID}" out || true
