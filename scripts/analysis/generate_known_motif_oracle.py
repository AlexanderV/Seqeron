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

Mismatch mode (A1-20, GenomicAnalyzer.FindKnownMotifHits `maxMismatches`):
  python3 generate_known_motif_oracle.py --mismatches > known_motif_mismatch_oracle.tsv
For every (sequence, pattern, k) with 1 <= k < len(pattern) and len(sequence) >= len(pattern) the hits of
`fuzznuc -pmismatch k -complement Y` (position + the `*mismatch n` note, none = 0) are asserted equal to an independent
Python Hamming scan (a window position agrees when the base is in the pattern code's NC-IUB base set; reverse hits =
the same scan with the IUPAC reverse-complement pattern). Sequences include copies of the pattern with 0-3 planted
substitutions on both strands. Output TSV (embedded in the test project):
  sequence  pattern  k  forward(`pos:mismatches` comma list | -)  reverse(`pos:mismatches` | -)
Not part of the oracle: fuzznuc reports hits that overhang the end of a sequence SHORTER than the pattern (e.g.
pattern AAA, -pmismatch 2 reports one hit on the sequence `AC`); every other fuzznuc hit equals the Hamming scan.
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


if '--mismatches' in sys.argv:
    rnd2 = random.Random(20261120)

    def ham(w, p):
        return sum(1 for a, b in zip(w, p) if a not in IUPAC[b])

    def find_k(seq, p, k):
        return [(i, ham(seq[i:i + len(p)], p)) for i in range(len(seq) - len(p) + 1) if ham(seq[i:i + len(p)], p) <= k]

    def concrete(p):
        return ''.join(rnd2.choice(IUPAC[c]) for c in p)

    def mutate(w, j):
        w = list(w)
        for i in rnd2.sample(range(len(w)), min(j, len(w))):
            w[i] = rnd2.choice([b for b in 'ACGT' if b != w[i]])
        return ''.join(w)

    pats = ['GAATTC', 'TATAAA', 'GGATCC', 'CCCGGG', 'ACGT', 'TTGACA', 'TATAAT', 'GCCACC', 'AAAAA', 'ACACAC',
            'TATAWAW', 'GCCRCC', 'WGATAR', 'RGATCY', 'GGNCC', 'RYRY', 'AGCTNN', 'CCWGG', 'TGASTCA', 'GGGRNWYYCC']
    seqs = []
    for _ in range(36):
        seqs.append(rand(rnd2.randint(12, 80), 'ACGT'))
    for _ in range(8):
        seqs.append(rand(rnd2.randint(10, 40), 'AC'))
    seqs += ['AAAAAAAAAA', 'GAATTC', 'GAATTCGAATTC', 'ACGTACGTACGT', 'TTTTTTTT', 'GGGATCGGG']
    for p in pats:
        for _ in range(3):
            w = mutate(concrete(p), rnd2.randint(0, 3))
            if rnd2.random() < 0.5:
                w = rc(w)
            seqs.append(rand(rnd2.randint(0, 15), 'ACGT') + w + rand(rnd2.randint(0, 15), 'ACGT'))
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 's.fa')
        with open(fa, 'w') as f:
            for i, sq in enumerate(seqs):
                f.write('>s%d\n%s\n' % (i, sq))
        n_rows = n_fwd = n_rev = n_mm = 0
        for p in pats:
            for k in (1, 2, 3):
                if k >= len(p):
                    continue
                out = subprocess.run(['fuzznuc', '-sequence', fa, '-pattern', p, '-pmismatch', str(k), '-complement', 'Y',
                                      '-outfile', 'stdout', '-rformat', 'gff'], capture_output=True, text=True, check=True).stdout
                fz = {i: ([], []) for i in range(len(seqs))}
                for line in out.splitlines():
                    if line.startswith('#') or not line.strip() or '\t' not in line:
                        continue
                    c = line.split('\t')
                    m = re.search(r'\*mismatch (\d+)', c[8])
                    fz[int(c[0][1:])][0 if c[6] == '+' else 1].append((int(c[3]) - 1, int(m.group(1)) if m else 0))
                for i, sq in enumerate(seqs):
                    if len(sq) < len(p):
                        continue
                    fwd, rev = sorted(fz[i][0]), sorted(fz[i][1])
                    assert fwd == find_k(sq, p, k), (sq, p, k, fwd, find_k(sq, p, k))
                    assert rev == find_k(sq, rc(p), k), (sq, p, k, rev, find_k(sq, rc(p), k))
                    n_rows += 1
                    n_fwd += bool(fwd)
                    n_rev += bool(rev)
                    n_mm += any(x[1] > 0 for x in fwd + rev)
                    f2 = ','.join('%d:%d' % x for x in fwd) or '-'
                    r2 = ','.join('%d:%d' % x for x in rev) or '-'
                    print('\t'.join([sq, p, str(k), f2, r2]))
    print('rows', n_rows, 'with fwd hits', n_fwd, 'rev hits', n_rev, 'with mismatching hits', n_mm, file=sys.stderr)
    sys.exit(0)

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
