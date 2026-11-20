#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/muscle3_default_iters_oracle.tsv (B09 audit A3-7) with the installed muscle3 3.8.1551 (`muscle3`, Debian package muscle3).

Question: what is the default of `-maxiters`? Source (muscle3 3.8.1551-2, params.cpp line 109): `unsigned g_uMaxIters = 8;`; the program's own `-help` text says "default 16" (stale).
The binary decides it: for the protein families below `muscle3 -in x.fa` (no -maxiters) writes the rows of `-maxiters 8` and NOT those of `-maxiters 16` (the two differ on the first family).
Rows of every run are recorded in INPUT order (the program writes tree order; the names map them back). Columns: id, variant (default | 8 | 16 | 2), sequences (';'), rows (';').
Run: python3 scripts/alignment/generate_muscle3_default_iters_oracle.py
"""
import os, random, subprocess, tempfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'muscle3_default_iters_oracle.tsv')
P = 'ACDEFGHIKLMNPQRSTVWY'


def family(rng, n, length, div):
    base = ''.join(rng.choice(P) for _ in range(length))
    seqs = []
    for _ in range(n):
        s = []
        for ch in base:
            x = rng.random()
            if x < 0.05:
                continue
            if x < 0.10:
                s.extend([ch] + [rng.choice(P) for _ in range(rng.randint(1, 3))])
                continue
            s.append(rng.choice(P) if rng.random() < div else ch)
        seqs.append(''.join(s))
    return seqs


def run(seqs, variant):
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, 'x.fa'), 'w') as f:
            for i, s in enumerate(seqs):
                f.write('>s%d\n%s\n' % (i, s))
        cmd = ['muscle3', '-in', 'x.fa', '-quiet'] + ([] if variant == 'default' else ['-maxiters', variant])
        out = subprocess.run(cmd, cwd=d, capture_output=True, text=True).stdout
    rows = {}
    name = None
    for line in out.splitlines():
        if line.startswith('>'):
            name = line[1:]
            rows[name] = ''
        else:
            rows[name] += line.strip()
    return [rows['s%d' % i] for i in range(len(seqs))]


rng = random.Random(77)
# the search that found the first family: 40 sequences, 300 residues, 50 % substitutions (seed 77, 10th draw); here the sets are drawn directly
sets = []
for k in range(40):
    n = rng.choice([20, 30, 40]); L = rng.choice([120, 200, 300]); div = rng.choice([0.3, 0.5, 0.65])
    seqs = family(rng, n, L, div)
    if k in (9, 27):
        sets.append(('D%02d' % k, seqs))
with open(out_path, 'w', newline='\n') as fh:
    fh.write('# ' + subprocess.run(['muscle3', '-version'], capture_output=True, text=True).stdout.strip() + '; uname -m ' + os.uname().machine + '\n')
    fh.write('id\tvariant\tsequences\trows\n')
    for sid, seqs in sets:
        for variant in ('default', '8', '16', '2'):
            fh.write('\t'.join([sid, variant, ';'.join(seqs), ';'.join(run(seqs, variant))]) + '\n')
print(len(sets), 'sets')
