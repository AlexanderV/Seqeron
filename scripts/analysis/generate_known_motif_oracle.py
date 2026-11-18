#!/usr/bin/env python3
"""GENOMIC-MOTIFS-001 oracle: EMBOSS fuzznuc 6.6.0 vs independent Python regex, for known-motif search.

Random DNA (alphabets ACGT / AC / A, tandem-rich), exact and IUPAC patterns (1-8 nt). For every
(sequence, pattern) the reference positions are computed twice and asserted equal:
  * fuzznuc -pattern P -pmismatch 0 -complement Y (GFF): '+' hits = forward window starts, '-' hits =
    windows matching the reverse complement (forward-strand coordinates, 1-based in the report -> 0-based here);
  * Python `re` with a zero-width lookahead (overlapping), IUPAC classes expanded from the NC-IUB table;
    reverse hits = forward hits of the IUPAC reverse-complement pattern.
Output TSV, embedded as tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/known_motif_oracle.tsv:
  sequence  pattern  forward(0-based starts, comma | -)  reverse(0-based forward-strand window starts | -)
Requires the `fuzznuc` binary (Debian package emboss). Usage: python3 generate_known_motif_oracle.py > known_motif_oracle.tsv
"""
import os, random, re, subprocess, sys, tempfile

rnd = random.Random(20261118)
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT',
         'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
        'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}


def rc(p):
    return ''.join(COMP[c] for c in reversed(p))


def regex(p):
    return re.compile('(?=(' + ''.join('[' + IUPAC[c] + ']' for c in p) + '))')


def find(seq, p):
    return [m.start() for m in regex(p).finditer(seq)]


def rand(n, alpha):
    return ''.join(rnd.choice(alpha) for _ in range(n))


seqs = []
for _ in range(40):
    seqs.append(rand(rnd.randint(8, 70), 'ACGT'))
for _ in range(10):
    seqs.append(rand(rnd.randint(8, 40), 'AC'))
seqs += ['AAAAA', 'AAAAAAAAAA', 'GAATTCAAAGAATTCTTTATATAAAGTTCAAGAAC', 'ACACACACACAC', 'GAATTC', 'ACGT', 'TTTT', 'GGGATCGGG']
for _ in range(10):
    unit = rand(rnd.randint(1, 4), 'ACGT')
    seqs.append(unit * rnd.randint(3, 12) + rand(rnd.randint(0, 6), 'ACGT'))

exact = ['GAATTC', 'AAA', 'AA', 'A', 'ACGT', 'TATAAA', 'AATAAA', 'GGATCC', 'CCCGGG', 'AT', 'TTT', 'ACAC', 'GCCACC',
         'CAAT', 'GGGCGG', 'TTCAAG']
degenerate = ['TATAWAW', 'GCCRCC', 'WGATAR', 'RGATCY', 'YGGCCR', 'NNN', 'ANA', 'KMS', 'BDHV', 'GGNCC', 'RYRY', 'SWSW',
              'AGCTNN', 'CCWGG', 'GCNNGC', 'RRR', 'YYY', 'N', 'R', 'ACNGT']
patterns = exact + degenerate

with tempfile.TemporaryDirectory() as d:
    fa = os.path.join(d, 's.fa')
    with open(fa, 'w') as f:
        for i, s in enumerate(seqs):
            f.write('>s%d\n%s\n' % (i, s))
    rows = []
    for p in patterns:
        out = subprocess.run(['fuzznuc', '-sequence', fa, '-pattern', p, '-pmismatch', '0', '-complement', 'Y',
                              '-outfile', 'stdout', '-rformat', 'gff'], capture_output=True, text=True, check=True).stdout
        fz = {i: ([], []) for i in range(len(seqs))}
        for line in out.splitlines():
            if line.startswith('#') or not line.strip() or '\t' not in line:
                continue
            c = line.split('\t')
            i = int(c[0][1:])
            start0 = int(c[3]) - 1
            fz[i][0 if c[6] == '+' else 1].append(start0)
        for i, s in enumerate(seqs):
            fwd, rev = sorted(fz[i][0]), sorted(fz[i][1])
            assert fwd == find(s, p), (s, p, fwd, find(s, p))
            assert rev == find(s, rc(p)), (s, p, rev, find(s, rc(p)))
            rows.append((s, p, fwd, rev))

for s, p, fwd, rev in rows:
    print('\t'.join([s, p, ','.join(map(str, fwd)) or '-', ','.join(map(str, rev)) or '-']))
print('rows', len(rows), 'with hits fwd', sum(1 for r in rows if r[2]), 'rev', sum(1 for r in rows if r[3]), file=sys.stderr)
