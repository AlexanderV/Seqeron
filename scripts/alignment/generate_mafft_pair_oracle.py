#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/mafft_pair_oracle.tsv (B09 audit A2-19, stage 1: two sequences) with the installed MAFFT 7.505 (`mafft`, Debian 7.505-1).

Every row = one run `mafft --quiet --thread 1 [--anysymbol] [options] in.fa` on a deterministic pair (python `random.Random`, fixed seed):
  * amino-acid pairs (related with substitutions / indels, unrelated, 1 - 400 residues), nucleotide pairs (DNA / RNA, IUPAC codes, lower case), mixed alphabets (a protein against a nucleotide sequence),
    unusual symbols (digits, '*', '#', '.', letters outside the alphabets) with --anysymbol, '-' characters in the input;
  * the options of the alignment of a pair: --op, --ep, --bl 30 / 45 / 50 / 62 / 80 / 90 / 100 / 0, --nuc, --amino, and the accepted ones (--retree 1 / 2, --nofft, --fft, --reorder).
Columns: id, options, seqs ('|' separated; the two input sequences), rc, rows ('|' separated).  Run: python3 scripts/alignment/generate_mafft_pair_oracle.py
"""
import os, random, subprocess, tempfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'mafft_pair_oracle.tsv')
rng = random.Random(7505)
P = 'ACDEFGHIKLMNPQRSTVWY'
N = 'ACGT'


def mut(s, alpha, r, ind):
    out = []
    for c in s:
        x = rng.random()
        if x < ind / 2:
            continue
        if x < ind:
            out.append(c); out.append(rng.choice(alpha)); continue
        out.append(rng.choice(alpha) if rng.random() < r else c)
    return ''.join(out) or s[:1]


def rand(alpha, n):
    return ''.join(rng.choice(alpha) for _ in range(n))


def related(alpha, n, r=None, ind=None):
    b = rand(alpha, n)
    return b, mut(b, alpha, rng.choice([0.05, 0.2, 0.4]) if r is None else r, rng.choice([0.02, 0.08, 0.15]) if ind is None else ind)


cases = []  # (options, [a, b])
lens_p = [1, 2, 3, 5, 8, 12, 20, 30, 45, 70, 100, 150, 220, 300]
lens_n = [1, 2, 4, 8, 15, 25, 40, 60, 100, 160, 250, 400]
for k in range(220):
    n = rng.choice(lens_p)
    a, b = related(P, n) if k % 5 else (rand(P, n), rand(P, rng.choice(lens_p)))
    cases.append(([], [a, b]))
for k in range(200):
    n = rng.choice(lens_n)
    a, b = related(N, n) if k % 5 else (rand(N, n), rand(N, rng.choice(lens_n)))
    if k % 4 == 3:
        a = a.replace('T', 'U'); b = b.replace('T', 'U')
    if k % 6 == 2:
        a = a.lower()
    cases.append(([], [a, b]))
# IUPAC / unusual symbols (--anysymbol is on in every run of this file)
for k in range(80):
    n = rng.choice([6, 10, 20, 40, 80])
    a, b = related(N, n)
    a = list(a); b = list(b)
    for s in (a, b):
        for _ in range(rng.randrange(1, 4)):
            s[rng.randrange(len(s))] = rng.choice('NRYKMSWBDHVXnryx')
    cases.append(([], [''.join(a), ''.join(b)]))
for k in range(80):
    n = rng.choice([6, 10, 20, 40, 80])
    a, b = related(P, n)
    a = list(a); b = list(b)
    for s in (a, b):
        for _ in range(rng.randrange(1, 4)):
            s[rng.randrange(len(s))] = rng.choice('XBZJOU*#1.@x')
    cases.append(([], [''.join(a), ''.join(b)]))
# mixed alphabets, '-' in the input, lower / upper mixtures
for k in range(30):
    cases.append(([], [rand(P, rng.choice([5, 12, 30])), rand(N, rng.choice([5, 12, 30]))]))
for k in range(30):
    a, b = related(P, rng.choice([10, 30, 60]))
    a = list(a); a.insert(rng.randrange(len(a)), '-'); a.insert(rng.randrange(len(a)), '-')
    cases.append(([], [''.join(a), b.lower() if k % 2 else b]))
# options
for k in range(160):
    prot = k % 2 == 0
    n = rng.choice([10, 30, 60, 120] if prot else [15, 40, 90, 200])
    a, b = related(P if prot else N, n)
    opt = []
    r = k % 8
    if r == 0:
        opt = ['--op', rng.choice(['0.5', '1.0', '2.0', '3.0', '5.0'])]
    elif r == 1:
        opt = ['--ep', rng.choice(['0.0', '0.1', '0.5', '1.0', '2.0', '-0.5'])]
    elif r == 2:
        opt = ['--op', rng.choice(['0.8', '1.5', '2.5']), '--ep', rng.choice(['0.0', '0.2', '0.7'])]
    elif r == 3:
        opt = ['--bl', str(rng.choice([30, 45, 50, 62, 80, 90, 100, 0]))]
    elif r == 4:
        opt = ['--nuc'] if prot else ['--amino']
        if not prot:
            a = a.replace('A', 'M'); b = b.replace('T', 'V')
    elif r == 5:
        opt = ['--retree', rng.choice(['1', '2'])]
    elif r == 6:
        opt = [rng.choice(['--nofft', '--fft', '--reorder'])]
    else:
        opt = ['--bl', str(rng.choice([30, 45, 80, 100])), '--op', rng.choice(['1.0', '2.0']), '--ep', rng.choice(['0.1', '0.3'])]
    cases.append((opt, [a, b]))
# the plain command (no --anysymbol): standard letters only
plain = []
for k in range(60):
    prot = k % 2 == 0
    n = rng.choice([8, 20, 50, 100])
    a, b = related(P if prot else N, n)
    plain.append((['PLAIN'], [a, b]))
cases += plain
# longer pairs (the FFT stage of the default FFT-NS-2 would start to matter here)
for k in range(24):
    prot = k % 2 == 0
    n = rng.choice([400, 600, 900, 1500])
    a, b = related(P if prot else N, n, r=rng.choice([0.05, 0.2]), ind=0.04)
    cases.append(([], [a, b]))

rows = []
fails = 0
with tempfile.TemporaryDirectory() as d:
    for k, (opts, seqs) in enumerate(cases):
        fin = os.path.join(d, 'in.fa')
        with open(fin, 'w', newline='\n') as f:
            for i, s in enumerate(seqs):
                f.write('>s%d\n%s\n' % (i + 1, s))
        plainrun = opts == ['PLAIN']
        cmd = ['mafft', '--quiet', '--thread', '1'] + ([] if plainrun else ['--anysymbol']) + ([] if plainrun else opts) + [fin]
        p = subprocess.run(cmd, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=300)
        out_rows = []
        if p.returncode == 0:
            for line in p.stdout.split('\n'):
                if line.startswith('>'):
                    out_rows.append('')
                elif line:
                    out_rows[-1] += line
        else:
            fails += 1
        rows.append(['M%04d' % (k + 1), 'plain' if plainrun else (' '.join(opts) or '-'), '|'.join(seqs), str(p.returncode), '|'.join(out_rows)])
ver = subprocess.run(['mafft', '--version'], capture_output=True, text=True).stderr.strip() or subprocess.run(['mafft', '--version'], capture_output=True, text=True).stdout.strip()
with open(out_path, 'w', newline='\n') as fh:
    fh.write('# ' + ver + '; uname -m ' + os.uname().machine + '\n')
    fh.write('id\toptions\tseqs\trc\trows\n')
    for r in rows:
        fh.write('\t'.join(r) + '\n')
print(len(rows), 'rows;', fails, 'failed;', ver)
