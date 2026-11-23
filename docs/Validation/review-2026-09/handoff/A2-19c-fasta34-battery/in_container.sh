#!/usr/bin/env bash
# Runs inside the image: (1) the raw fasta34 battery (inputs/battery.tsv + the 20 probes of getfastascoremtx): every call as splittbfast.c makes it; (2) `mafft <options> in.fa` (WITHOUT --quiet) for every set of inputs/sets.tsv.
set -uo pipefail
cd /work
rm -rf out; mkdir -p out/battery out/sets
{ mafft --version 2>&1 | head -1; command -v fasta34; sha256sum "$(command -v fasta34)"; uname -m; cat /opt/fasta34.build.log | tail -3; ls /opt/FASTA34_BUILD_FAILED 2>&1; } > out/versions.txt 2>&1
python3 -I - <<'PY'
import gzip, hashlib, os, subprocess, tempfile


def fasta_text(items):
    # writeData_pointer: ">" + name (the first character of the name is skipped) + newline, 60 letters per line
    return ''.join('>%s\n' % n + ''.join(s[j:j + 60] + '\n' for j in range(0, len(s), 60)) for n, s in items)


def call(kind, query, lib, tag):
    with tempfile.TemporaryDirectory() as td:
        d, q, r = (os.path.join(td, x) for x in ('data', 'query', 'result'))
        open(d, 'w').write(fasta_text(('+===========+%d                      ' % i, s) for i, s in enumerate(lib)))
        open(q, 'w').write(fasta_text([('+==========+0                      ', query)]))
        n = len(lib)
        if kind == 'p':
            cmd = 'fasta34 -z3 -m10 -p -Q -H -b%d -E%d -d0 %s %s 2 > %s' % (n, n, q, d, r)
        else:
            cmd = 'fasta34 -z3 -m10 -n -Q -H -b%d -E%d -d0 %s %s 6 > %s' % (n, n, q, d, r)
        rc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        out = open(r, 'rb').read()
    with gzip.open('out/battery/%s.out.gz' % tag, 'wb') as f:
        f.write(out)
    return rc.returncode, hashlib.sha256(out).hexdigest(), rc.stderr.decode('latin-1')[:200]


res = []
# the 20 probes of getfastascoremtx (protein): AAAAAA ... against 'AAAAAAXXXXXXCCCCCCXXXXXX...'
lib = ''.join(c * 6 + 'X' * 6 for c in 'ACDEFGHIKLMNPQRSTVWY')
for c in 'ACDEFGHIKLMNPQRSTVWY':
    rc, h, e = call('p', c * 6, [lib], 'probe-' + c)
    res.append('\t'.join(['probe-' + c, 'p', str(rc), h, e]))
for line in open('inputs/battery.tsv'):
    k, kind, q, lib = line.rstrip('\n').split('\t')
    rc, h, e = call(kind, q, lib.split('|'), 'b' + k)
    res.append('\t'.join(['b' + k, kind, str(rc), h, e]))
open('out/battery.tsv', 'w').write('\n'.join(res) + '\n')
print('battery', len(res))

rows = []
for line in open('inputs/sets.tsv'):
    k, opt, seqs = line.rstrip('\n').split('\t')
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, 'in.fa')
        open(fa, 'w').write(''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs.split('|'))))
        p = subprocess.run(['mafft'] + opt.split() + [fa], capture_output=True, cwd=td)
        tf = os.path.join(td, 'in.fa.tree')
        tree = open(tf, 'rb').read().decode('latin-1') if os.path.exists(tf) else ''
    out, names = [], []
    for l in p.stdout.decode('latin-1').splitlines():
        if l.startswith('>'):
            names.append(l[1:])
            out.append('')
        elif l.strip():
            out[-1] += l.strip()
    with gzip.open('out/sets/%s.stderr.gz' % k, 'wb') as f:
        f.write(p.stderr)
    esc = lambda s: s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t')
    rows.append('\t'.join([k, opt, seqs, str(p.returncode), '|'.join(out), '|'.join(names), hashlib.sha256(p.stderr).hexdigest(), esc(tree)]))
open('out/results.tsv', 'w').write('\n'.join(rows) + '\n')
print('sets', len(rows))
PY
[ -n "${HOST_UID:-}" ] && chown -R "$HOST_UID:${HOST_GID:-$HOST_UID}" out || true
