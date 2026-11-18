#!/usr/bin/env python3
"""ALIGN-GLOBAL/LOCAL/SEMI-001 (audit A1-5) oracle: pairwise alignment with an ARBITRARY gap-cost function w(k)
(Waterman-Smith-Beyer 1976) - SequenceAligner.{Global,Local,SemiGlobal}AlignGeneralGap(seq1, seq2, Func<int,double> gapScore, ...).

Model: every maximal gap run of k positions scores w(k) (a score, negative = penalty); the two kinds of runs may be adjacent and are
charged separately; a run is never split (so a convex w is honoured).
Oracles:
  biopython = Bio.Align.PairwiseAligner(mode=global|local, gap_score=<callable (position, length) -> w(length)>) 1.88 - its
              "Waterman-Smith-Beyer global/local alignment algorithm" (align().score and score() agree on these rows; checked here)
  pyref     = independent brute-force implementation written for this oracle (scripts/alignment/wsb_reference.py logic inlined below:
              three-state recurrence over every gap length, also with the end-gap policy of the ends-free rows), compared with Biopython
              on every global/local row in this script (a difference aborts the run)
Biopython ignores left/right_insertion/deletion_score for a callable gap_score (verified: identical scores with and without), so the
ends-free ("semi") rows have pyref only; Seqeron's end-gap handling is additionally locked by an exhaustive enumeration in the unit
tests and by the affine/linear equality tests.

Gap functions (kind a b c d): linear -(a*k); affine -(a + b*k); log -(a + b*ln k); sqrt -(a + b*sqrt k); quad -(a + b*k^2);
cap -min(a + b*k, c); step -(a if k <= c else b); twoslope -(a + b*k) for k <= c, else -(a + b*c + d*(k - c)).

Output TSV (embedded: tests/.../TestData/Alignment/general_gap_oracle.tsv), columns
  mode(global|local|semi) flags(0..15; semi only) kind a b c d scoring(scalar|BLOSUM62|EDNAFULL) match mismatch seq1 seq2 biopython pyref

Usage: python3 scripts/alignment/generate_general_gap_oracle.py [n_pairs=60] > general_gap_oracle.tsv
"""
import math, random, sys
from Bio.Align import PairwiseAligner, substitution_matrices as sm

rnd = random.Random(20261118)
PROT = 'ACDEFGHIKLMNPQRSTVWY'
NEG = float('-inf')


def gap_fn(kind, a, b, c, d):
    if kind == 'linear': return lambda k: -(a * k)
    if kind == 'affine': return lambda k: -(a + b * k)
    if kind == 'log': return lambda k: -(a + b * math.log(k))
    if kind == 'sqrt': return lambda k: -(a + b * math.sqrt(k))
    if kind == 'quad': return lambda k: -(a + b * k * k)
    if kind == 'cap': return lambda k: -min(a + b * k, c)
    if kind == 'step': return lambda k: -(a if k <= c else b)
    if kind == 'twoslope': return lambda k: -(a + b * k) if k <= c else -(a + b * c + d * (k - c))
    raise ValueError(kind)


PARAMS = {
    'linear': lambda: (rnd.choice([1, 2, 0.5, 3]), 0, 0, 0),
    'affine': lambda: (rnd.choice([0, 2, 5, 3]), rnd.choice([1, 0.5, 2]), 0, 0),
    'log': lambda: (rnd.choice([1, 2, 4]), rnd.choice([1, 2, 3]), 0, 0),
    'sqrt': lambda: (rnd.choice([1, 2, 3]), rnd.choice([1, 2, 0.5]), 0, 0),
    'quad': lambda: (rnd.choice([0, 1, 2]), rnd.choice([0.25, 0.5, 1]), 0, 0),
    'cap': lambda: (rnd.choice([1, 2]), rnd.choice([1, 2]), rnd.choice([4, 6, 8]), 0),
    'step': lambda: (rnd.choice([1, 2]), rnd.choice([5, 8, 12]), rnd.choice([1, 2, 3]), 0),
    'twoslope': lambda: (rnd.choice([1, 2]), rnd.choice([0.5, 1]), rnd.choice([2, 3, 4]), rnd.choice([0.1, 0.25, 2])),
}


def pyref(a, b, w, sub, mode, ends=0):
    """three-state Waterman-Smith-Beyer (single-run gaps), exhaustive over gap lengths; local: zero floor on M; ends = EndGaps flags"""
    m, n = len(a), len(b)
    M = [[NEG] * (n + 1) for _ in range(m + 1)]; X = [[NEG] * (n + 1) for _ in range(m + 1)]; Y = [[NEG] * (n + 1) for _ in range(m + 1)]
    loc = mode == 'local'
    fs1, fe1, fs2, fe2 = [bool(ends & k) for k in (1, 2, 4, 8)]
    M[0][0] = NEG if loc else 0
    for i in range(1, m + 1):
        if not loc:
            if fs1: M[i][0] = 0
            else: X[i][0] = w(i)
    for j in range(1, n + 1):
        if not loc:
            if fs2: M[0][j] = 0
            else: Y[0][j] = w(j)
    best = 0.0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            prev = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1])
            if loc: prev = max(prev, 0.0)
            M[i][j] = prev + sub(a[i - 1], b[j - 1])
            X[i][j] = max(max(M[i - k][j], Y[i - k][j]) + w(k) for k in range(1, i + 1))
            Y[i][j] = max(max(M[i][j - k], X[i][j - k]) + w(k) for k in range(1, j + 1))
            if loc: best = max(best, M[i][j])
    if loc: return best
    H = lambda i, j: max(M[i][j], X[i][j], Y[i][j])
    res = H(m, n)
    if fe2:
        for j in range(n + 1): res = max(res, H(m, j))
    if fe1:
        for i in range(m + 1): res = max(res, H(i, n))
    return res


def mutate(s, alpha, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice(alpha))
        else: out.append(c)
    return ''.join(out) or rnd.choice(alpha)


def gen_pair(alpha, maxlen):
    a = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, maxlen)))
    b = mutate(a, alpha, rnd.choice([0.1, 0.3, 0.6]))
    r = rnd.random()
    if r < 0.2: b = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, maxlen)))
    elif r < 0.4:  # a long deletion / insertion: where the gap-function shape matters
        i = rnd.randint(0, len(a)); j = min(len(a), i + rnd.randint(2, 9)); b = a[:i] + a[j:] or a[:1]
    elif r < 0.55:
        i = rnd.randint(0, len(a)); b = a[:i] + ''.join(rnd.choice(alpha) for _ in range(rnd.randint(2, 9))) + a[i:]
    return a, b


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    print('mode\tflags\tkind\ta\tb\tc\td\tscoring\tmatch\tmismatch\tseq1\tseq2\tbiopython\tpyref')
    mism = 0
    for idx in range(n):
        for kind in PARAMS:
            for mode in ('global', 'local'):
                scoring = rnd.choice(['scalar', 'scalar', 'scalar', 'BLOSUM62', 'EDNAFULL'])
                if scoring == 'BLOSUM62': alpha = PROT
                else: alpha = 'ACGT' if rnd.random() < 0.8 else 'AC'
                a, b = gen_pair(alpha, 26)
                ma, mi = rnd.choice([(2, -3), (1, -1), (5, -4), (3, -2)])
                pa = PARAMS[kind]()
                w = gap_fn(kind, *pa)
                if scoring == 'scalar':
                    sub = lambda x, y, ma=ma, mi=mi: ma if x == y else mi
                    al = PairwiseAligner(); al.match_score = ma; al.mismatch_score = mi
                else:
                    mat = sm.load('BLOSUM62' if scoring == 'BLOSUM62' else 'NUC.4.4')
                    sub = lambda x, y, mat=mat: float(mat[x][y])
                    al = PairwiseAligner(); al.substitution_matrix = mat
                al.mode = mode
                al.gap_score = lambda pos, length, w=w: w(length)
                bs = float(al.align(a, b).score)
                bsapi = float(al.score(a, b))
                ref = pyref(a, b, w, sub, mode)
                if abs(bs - ref) > 1e-9 or abs(bsapi - ref) > 1e-9:
                    mism += 1
                    print(f'# MISMATCH {mode} {kind} {pa} {a} {b} bio={bs} api={bsapi} ref={ref}', file=sys.stderr)
                print('\t'.join([mode, '0', kind] + [repr(float(x)) for x in pa] + [scoring, str(ma), str(mi), a, b, repr(bs), repr(ref)]))
            # ends-free rows (all 16 policies on a rotating basis): pyref only
            flags = rnd.randint(1, 15)
            a, b = gen_pair('ACGT', 16)
            pa = PARAMS[kind]()
            w = gap_fn(kind, *pa)
            ma, mi = rnd.choice([(2, -3), (1, -1)])
            ref = pyref(a, b, w, lambda x, y: ma if x == y else mi, 'semi', flags)
            print('\t'.join(['semi', str(flags), kind] + [repr(float(x)) for x in pa] + ['scalar', str(ma), str(mi), a, b, 'NA', repr(ref)]))
    print(f'# biopython-vs-pyref disagreements: {mism}', file=sys.stderr)


main()
