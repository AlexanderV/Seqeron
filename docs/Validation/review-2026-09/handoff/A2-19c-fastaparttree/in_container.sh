#!/usr/bin/env bash
# Runs inside the image: for every set of inputs/sets.tsv `mafft <options> in.fa` (the program's alignment rows, return code, SHA-256 of standard error) -> out/results.tsv, plus the versions.
set -uo pipefail
cd /work
mkdir -p out
{ mafft --version 2>&1 | head -1; command -v fasta34; sha256sum "$(command -v fasta34)"; cat /opt/fasta34.build.log; ls /opt/FASTA34_BUILD_FAILED 2>&1; } > out/versions.txt 2>&1
python3 -I - <<'PY'
import hashlib, os, subprocess, tempfile
rows = []
for line in open('inputs/sets.tsv'):
    k, opt, seqs = line.rstrip('\n').split('\t')
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, 'in.fa')
        open(fa, 'w').write(''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs.split('|'))))
        p = subprocess.run(['mafft', '--quiet'] + opt.split() + [fa], capture_output=True)
    out = []
    for l in p.stdout.decode('latin-1').splitlines():
        if l.startswith('>'):
            out.append('')
        elif l.strip():
            out[-1] += l.strip()
    rows.append('\t'.join([k, opt, seqs, str(p.returncode), '|'.join(out), hashlib.sha256(p.stderr).hexdigest(), p.stderr.decode('latin-1')[:200].replace('\n', '\\n')]))
open('out/results.tsv', 'w').write('\n'.join(rows) + '\n')
print('wrote', len(rows), 'rows')
PY
[ -n "${HOST_UID:-}" ] && chown -R "$HOST_UID:${HOST_GID:-$HOST_UID}" out || true
