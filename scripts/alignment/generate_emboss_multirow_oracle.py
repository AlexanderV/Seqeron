#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/emboss_multirow_oracle.tsv (B09 audit A2-27d): the EMBOSS 6.6.0 alignment output formats of three or more rows, written by the installed `aligncopy`.

Every row = one run `aligncopy -sequence a.fa -aformat FORMAT [-awidth W] -outfile stdout -auto` on the aligned FASTA of the row. Columns: id, format, width (0 = default), the aligned FASTA, the output (run date replaced by the
fixed date `Sat 21 Nov 2026 19:00:00`) or `ERROR:` + standard error when the program refuses (the pair-markup formats of three or more rows). Text columns escaped (backslash, newline, tab, carriage return).
Run: python3 scripts/alignment/generate_emboss_multirow_oracle.py
"""
import os
import random
import re
import subprocess
import tempfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'emboss_multirow_oracle.tsv')
rng = random.Random(20261322)
FORMATS = ['unknown', 'multiple', 'simple', 'srs', 'pair', 'srspair', 'markx0', 'markx3', 'match', 'score', 'sam']
PROT = 'ACDEFGHIKLMNPQRSTVWY'


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def make(prot):
    n = rng.choice([3, 3, 4, 5, 6, 8])
    cols = rng.choice([10, 25, 51, 60, 130])
    alpha = PROT if prot else 'ACGT'
    base = ''.join(rng.choice(alpha) for _ in range(cols))
    rows = []
    for _ in range(n):
        r = []
        for c in base:
            x = rng.random()
            if x < 0.12:
                r.append('-')
            elif x < 0.30:
                r.append(rng.choice(alpha))
            else:
                r.append(c)
        s = ''.join(r)
        if rng.random() < 0.3:
            s = ''.join(ch.lower() if rng.random() < 0.3 else ch for ch in s)
        rows.append(s)
    # no column of only gaps (the program drops nothing, but keep the data sane)
    for i in range(cols):
        if all(r[i] == '-' for r in rows):
            rows[0] = rows[0][:i] + base[i] + rows[0][i + 1:]
    names = []
    for k in range(n):
        names.append(rng.choice(['seq%d' % (k + 1), 'sp|P%05d|X_HUMAN' % rng.randrange(99999), 'a%d' % k, 'longer_name_%d_abcdefghijklmnop' % k]))
    names = [nm + ('' if names.index(nm) == i else '_%d' % i) for i, nm in enumerate(names)]
    return ''.join('>%s\n%s\n' % (names[i], rows[i]) for i in range(n))


rows = []
with tempfile.TemporaryDirectory() as d:
    for i in range(1200):
        prot = rng.random() < 0.5
        fa = make(prot)
        fmt = ['unknown', 'multiple', 'simple', 'srs'][i % 4] if i % 3 else FORMATS[4 + (i // 3) % 7]
        width = rng.choice([0, 0, 0, 20, 33, 50, 80, 200])
        open(os.path.join(d, 'a.fa'), 'w').write(fa)
        cmd = ['aligncopy', '-sequence', 'a.fa', '-aformat', fmt, '-outfile', 'stdout', '-auto']
        if width:
            cmd += ['-awidth', str(width)]
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
        if r.returncode != 0:
            out = 'ERROR:' + r.stderr
        else:
            out = re.sub(r'# Rundate: .*', '# Rundate: Sat 21 Nov 2026 19:00:00', r.stdout)
        rows.append((fmt, str(width), fa, out))

with open(out_path, 'w', newline='\n') as fh:
    fh.write('id\tformat\twidth\tfasta\toutput\n')
    for i, r in enumerate(rows):
        fh.write('M%04d\t%s\n' % (i + 1, '\t'.join(esc(c) for c in r)))
print(len(rows), 'rows;', sum(1 for r in rows if r[3].startswith('ERROR')), 'refusals')
