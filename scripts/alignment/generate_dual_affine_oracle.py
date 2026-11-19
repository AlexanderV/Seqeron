#!/usr/bin/env python3
"""ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit A2-29, B09 F123) oracle: ksw2-style two-piece affine gap alignment, global and local, full matrix and banded.

Seqeron API: SequenceAligner.GlobalAlignDualAffine / LocalAlignDualAffine(seq1, seq2, DualAffineScoring, band).
Reference: (1) Biopython 1.88 PairwiseAligner with the callable gap score g(k) = max(o1 + k e1, o2 + k e2) (Waterman-Smith-Beyer / generic-gap algorithms, mode global / local) for the full matrix; (2) an independent
Python Waterman-Smith-Beyer O(m n (m + n)) dynamic program written from the definition for every case, restricted to the cells |i - j| <= band for the banded cases (it must equal Biopython on the full-matrix cases:
asserted here). Output TSV (tests/.../TestData/Alignment/dual_affine_oracle.tsv):  id  mode  band  seq1  seq2  match  mismatch  open1  ext1  open2  ext2  score
Usage: python3 scripts/alignment/generate_dual_affine_oracle.py OUTFILE [casesPerVariant] [seed]
"""
import random, sys
from Bio.Align import PairwiseAligner

NEG = -10 ** 9


def wsb(a, b, M, X, g, local, band):
    m, n = len(a), len(b)
    ok = lambda i, j: band < 0 or abs(i - j) <= band
    H = [[NEG] * (n + 1) for _ in range(m + 1)]
    H[0][0] = 0
    best = 0
    for i in range(m + 1):
        for j in range(n + 1):
            if (i == 0 and j == 0) or not ok(i, j):
                continue
            v = NEG
            if i and j and H[i - 1][j - 1] > NEG:
                v = H[i - 1][j - 1] + (M if a[i - 1] == b[j - 1] else X)
            for k in range(1, j + 1):
                if H[i][j - k] > NEG: v = max(v, H[i][j - k] + g(k))
            for k in range(1, i + 1):
                if H[i - k][j] > NEG: v = max(v, H[i - k][j] + g(k))
            if local and v < 0: v = 0
            H[i][j] = v
            if local and v > best: best = v
    return best if local else H[m][n]


def mutate(rnd, s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate: out.append(rnd.choice('ACGT'))
        elif r < rate * 1.3: continue
        elif r < rate * 1.6: out.append(c); out.extend(rnd.choice('ACGT') for _ in range(rnd.randint(1, 8)))
        else: out.append(c)
    return ''.join(out)


def main(outfile, per, seed):
    rnd = random.Random(seed)
    rows = []
    for mode in ['global', 'local']:
        for banded in [False, True]:
            for _ in range(per):
                n = rnd.choice([4, 10, 20, 35, 50])
                a = ''.join(rnd.choice('ACGT') for _ in range(n))
                b = mutate(rnd, a, rnd.choice([0.05, 0.15, 0.3])) or 'A'
                if mode == 'local' and rnd.random() < 0.5: b = ''.join(rnd.choice('ACGT') for _ in range(rnd.choice([5, 20, 40])))
                M = rnd.randint(1, 5); X = -rnd.randint(1, 6)
                o1 = -rnd.randint(0, 8); e1 = -rnd.randint(0, 4); o2 = -rnd.randint(0, 30); e2 = -rnd.randint(0, 3)
                g = lambda k, o1=o1, e1=e1, o2=o2, e2=e2: max(o1 + k * e1, o2 + k * e2)
                band = -1
                if banded:
                    band = rnd.choice([0, 1, 2, 3, 5, 8, 15])
                    if mode == 'global': band = max(band, abs(len(a) - len(b)))
                s = wsb(a, b, M, X, g, mode == 'local', band)
                if not banded:
                    al = PairwiseAligner(); al.mode = mode
                    al.match_score = M; al.mismatch_score = X; al.gap_score = lambda pos, k, g=g: g(k)
                    bs = al.score(a, b)
                    assert abs(bs - s) < 1e-9, (a, b, M, X, o1, e1, o2, e2, bs, s)
                rows.append('\t'.join(map(str, [len(rows), mode, band, a, b, M, X, o1, e1, o2, e2, s])))
    with open(outfile, 'w') as f:
        f.write('id\tmode\tband\tseq1\tseq2\tmatch\tmismatch\topen1\text1\topen2\text2\tscore\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 400, int(sys.argv[3]) if len(sys.argv) > 3 else 20261124)
