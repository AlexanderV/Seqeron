#!/usr/bin/env python3
"""ALIGN-SEMI-001 oracle (A3-10): Biopython 1.88 PairwiseAligner with the twelve per-direction / per-end gap-score properties
({target,query}_{internal,left,right}_{open,extend}_gap_score).

Seqeron sequence1 = Biopython target, sequence2 = Biopython query (a gap in sequence 1 is a target gap: target_*; a gap in sequence 2 is a query gap: query_*).
Rows: random DNA pairs (alphabets A / AC / ACGT, lengths 1-16) x mode (global | local) x linear (open == extend for all twelve) | affine (independent open / extend,
open >= -6, extend >= -3, both <= 0, including open weaker than extend) x scalar (match / mismatch integers) | matrix (EDNAFULL / BLOSUM62, fractional scores);
roughly a third of the global rows have free (0) left / right ends (semi-global family), a third have all ends priced differently, a third a mixture.
`score` = PairwiseAligner.score(target, query) (Biopython 1.88).
TSV (tests/.../TestData/Alignment/gap_profile_oracle.tsv):
  seq1 \t seq2 \t mode(global|local) \t gaps(linear|affine) \t scoring(scalar:M:X | EDNAFULL | BLOSUM62) \t twelve scores (target_internal_open, target_internal_extend,
  target_left_open, target_left_extend, target_right_open, target_right_extend, query_internal_open, ..., query_right_extend; comma separated) \t biopython score \t model score
`model score` = an independent Python reference (three-state Gotoh DP, gap scores a function of the matrix position: a vertical run in column 0 is a start run, in the last
column an end run; a horizontal run in row 0 / the last row likewise; local: interior scores only) = the Seqeron semantics. It equals Biopython on every row except one documented
Biopython 1.88 corner (below); the generator aborts a row only if it differs for another reason.
Biopython 1.88 corner (pairwisealigner.c): the LINEAR Needleman-Wunsch score (all open == extend) initialises the first column of the last row as
`row[0] = nA * right_gap_extend_B` (the right / end score of the deletion direction), i.e. a leading run that consumes ALL of sequence 1 is charged like an end gap, whereas the
affine Gotoh routine of the same file initialises it with the LEFT score (`Ix_row[0] = left_gap_open_B + left_gap_extend_B * (i - 1)`). Seqeron follows the position rule (start run),
so these rows (an all-gap alignment whose first run spans the whole of sequence 1, only reachable when the end scores are cheaper than the start scores) differ from Biopython's linear
score; the reference marks them: model score != biopython score, and the `corner` variant of the model (that cell charged with the right score) reproduces Biopython.
Usage: python3 scripts/alignment/generate_gap_profile_oracle.py [rows] > gap_profile_oracle.tsv
"""
import random, sys
from Bio import Align
from Bio.Align import substitution_matrices

rnd = random.Random(20261121)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 900
PROPS = [f'{d}_{p}_{k}_gap_score' for d in ('target', 'query') for p in ('internal', 'left', 'right') for k in ('open', 'extend')]
PROTEIN = 'ACDEFGHIKLMNPQRSTVWY'
GAPVALS = [0, -0.5, -1, -1.5, -2, -2.5, -3, -4, -5, -6]


def seq(alpha, lo=1, hi=16):
    return ''.join(rnd.choice(alpha) for _ in range(rnd.randint(lo, hi)))


def pick_scores(affine, style, fractional):
    def ext():
        v = rnd.choice([0, -1, -2, -3] if not fractional else [0, -0.5, -1, -1.5, -2, -0.25, -0.1])
        return v
    def opn():
        return rnd.choice([0, -1, -2, -3, -4, -5, -6] if not fractional else [0, -0.5, -1, -2.5, -10, -3.5, -0.1])
    vals = {}
    for d in ('target', 'query'):
        for p in ('internal', 'left', 'right'):
            e = ext(); o = opn() if affine else e
            if not affine: o = e
            vals[f'{d}_{p}_open_gap_score'] = o
            vals[f'{d}_{p}_extend_gap_score'] = e
    for d in ('target', 'query'):
        for p in ('left', 'right'):
            if style == 'free' and rnd.random() < 0.7:
                vals[f'{d}_{p}_open_gap_score'] = 0
                vals[f'{d}_{p}_extend_gap_score'] = 0
            elif style == 'same':
                vals[f'{d}_{p}_open_gap_score'] = vals[f'{d}_internal_open_gap_score']
                vals[f'{d}_{p}_extend_gap_score'] = vals[f'{d}_internal_extend_gap_score']
    return vals


NEG = -1e18
NAMES = PROPS


def model(a, b, sub, vals, mode, corner=False):
    g = dict(zip(NAMES, vals))
    m, n = len(a), len(b)
    local = mode == 'local'

    def X(j, i=None):  # vertical run (gap in the query = sequence-1 residue), column j
        pos = 'internal' if local else 'left' if j == 0 else 'right' if j == n else 'internal'
        return g['query_%s_open_gap_score' % pos], g['query_%s_extend_gap_score' % pos]

    def Y(i):  # horizontal run (gap in the target), row i
        pos = 'internal' if local else 'left' if i == 0 else 'right' if i == m else 'internal'
        return g['target_%s_open_gap_score' % pos], g['target_%s_extend_gap_score' % pos]

    floor = 0 if local else NEG
    M = [[NEG] * (n + 1) for _ in range(m + 1)]
    Xs = [[NEG] * (n + 1) for _ in range(m + 1)]
    Ys = [[NEG] * (n + 1) for _ in range(m + 1)]
    M[0][0] = 0
    best = 0
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 and j == 0: continue
            if i > 0 and j > 0:
                M[i][j] = max(floor, M[i - 1][j - 1], Xs[i - 1][j - 1], Ys[i - 1][j - 1]) + sub(a[i - 1], b[j - 1])
                if local: best = max(best, M[i][j])
            if i > 0 and not (local and j == 0):
                o, e = X(j, i)
                Xs[i][j] = max(M[i - 1][j] + o, Ys[i - 1][j] + o, Xs[i - 1][j] + e)
                if corner and not local and j == 0 and i == m:  # NW linear corner: a column-0 run consuming ALL of sequence 1 is charged with the end (right) score
                    Xs[i][j] = g['query_right_open_gap_score'] + (m - 1) * g['query_right_extend_gap_score']
            if j > 0 and not (local and i == 0):
                o, e = Y(i)
                Ys[i][j] = max(M[i][j - 1] + o, Xs[i][j - 1] + o, Ys[i][j - 1] + e)
    if local: return best
    return max(M[m][n], Xs[m][n], Ys[m][n])


def main():
    edna = substitution_matrices.load('NUC.4.4')
    blosum = substitution_matrices.load('BLOSUM62')
    out = 0
    for i in range(N):
        mode = 'global' if rnd.random() < 0.8 else 'local'
        affine = rnd.random() < 0.6
        kind = rnd.random()
        fractional = False
        if kind < 0.7:
            alpha = rnd.choice(['A', 'AC', 'ACGT', 'ACGT'])
            a, b = seq(alpha), seq(alpha)
            m, x = rnd.choice([1, 2, 3, 5]), rnd.choice([0, -1, -2, -3])
            sc = 'scalar:%d:%d' % (m, x)
        elif kind < 0.85:
            a, b = seq('ACGT'), seq('ACGT')
            sc = 'EDNAFULL'; fractional = True
        else:
            a, b = seq(PROTEIN, 1, 20), seq(PROTEIN, 1, 20)
            sc = 'BLOSUM62'; fractional = True
        style = rnd.choice(['free', 'priced', 'priced', 'same'])
        vals = pick_scores(affine, style, fractional)
        al = Align.PairwiseAligner()
        al.mode = mode
        if sc.startswith('scalar'):
            _, m, x = sc.split(':')
            al.match_score = int(m); al.mismatch_score = int(x)
        else:
            al.substitution_matrix = edna if sc == 'EDNAFULL' else blosum
        for k, v in vals.items():
            d, p_, kind, _, _ = k.split('_')
            setattr(al, '%s_%s_%s_score' % (kind, p_, 'insertion' if d == 'target' else 'deletion'), v)  # Biopython 1.88 names (the target_* names are deprecated aliases)
        try:
            score = al.score(a, b)
        except Exception as ex:
            sys.stderr.write('FAIL %s %s %s: %s\n' % (a, b, vals, ex))
            continue
        if sc.startswith('scalar'):
            _, mm, xx = sc.split(':'); sub = lambda p, q, mm=int(mm), xx=int(xx): mm if p == q else xx
        else:
            mat = edna if sc == 'EDNAFULL' else blosum
            sub = lambda p, q, mat=mat: float(mat[p][q])
        pv = [float(vals[p]) for p in PROPS]
        ms = model(a, b, sub, pv, mode)
        if abs(ms - score) > 1e-9:
            mc = model(a, b, sub, pv, mode, corner=True)
            if affine or abs(mc - score) > 1e-9:
                sys.stderr.write('UNEXPLAINED %s %s %s %s %s: biopython %s model %s corner %s\n' % (a, b, mode, affine, pv, score, ms, mc))
                continue
            sys.stderr.write('corner row %s %s\n' % (a, b))
        print('\t'.join([a, b, mode, 'affine' if affine else 'linear', sc, ','.join(repr(v) for v in pv), repr(float(score)), repr(float(ms))]))
        out += 1
    sys.stderr.write('rows: %d\n' % out)


if __name__ == '__main__':
    main()
