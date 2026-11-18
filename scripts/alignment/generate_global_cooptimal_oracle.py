#!/usr/bin/env python3
"""ALIGN-GLOBAL-001 (audit A1-7b) oracle: ALL co-optimal global / ends-free alignments.

Seqeron API: SequenceAligner.GlobalAlignCoOptimal / GlobalAlignAffineCoOptimal(seq1, seq2, scoring, maxAlignments) and
SemiGlobalAlignCoOptimal / SemiGlobalAlignAffineCoOptimal(seq1, seq2, scoring, EndGaps, maxAlignments); EndGaps flag set
{FreeStart1=1, FreeEnd1=2, FreeStart2=4, FreeEnd2=8} (see generate_endgaps_oracle.py).

Definition (checked here against two independent references): an alignment is a gapped column sequence; its score = substitution
columns + gap cost (linear: GapExtend per gap position; affine: GapOpen + k*GapExtend per maximal run of gap columns of the SAME row,
a gap in row 1 directly followed by a gap in row 2 are two runs); the leading / trailing gap columns of a FREE end cost 0 (the first
column decides which side's leading flank it is; the last column decides the trailing flank). The co-optimal set = all column
sequences of maximum score, every distinct gapped-row pair once.
References:
  * exhaustive(): enumeration of EVERY column sequence of tiny inputs (len <= 5 both) with the score definition above (written from
    the definition, not from the C# code) -> asserted equal to Biopython's set on every tiny row.
  * Biopython 1.88 PairwiseAligner(mode='global').align(a, b): the set of all co-optimal alignments (rows rebuilt from the
    coordinates; left/right_insertion/deletion_score = 0 for the free ends; open_gap_score = GapOpen + GapExtend, extend_gap_score =
    GapExtend; linear: both = GapExtend) -> column `alignments` ("row1|row2;..." sorted; NA when more than CAP alignments).
Output TSV (embedded: tests/.../TestData/Alignment/global_cooptimal_oracle.tsv), columns
  kind(scalar|matrix) model(linear|affine) flags(0..15) matrix match mismatch gapOpen gapExtend seq1 seq2 biopython count alignments
Usage: python3 scripts/alignment/generate_global_cooptimal_oracle.py [n_pairs=120] > global_cooptimal_oracle.tsv
"""
import random, sys
from Bio.Align import PairwiseAligner, substitution_matrices as sm

rnd = random.Random(20261120)
CAP = 150
PROT = 'ACDEFGHIKLMNPQRSTVWY'
SCHEMES = [(1, -1, -1, -1), (2, -3, 0, -2), (1, -1, -2, -1), (1, -4, -1, -1), (0, 0, 0, 0), (1, 0, 0, -1), (3, -2, -1, 0), (2, -1, -3, -1),
           (1, -1, 0, 0), (5, -4, -10, -1), (2, -3, -3, -1), (1, -2, -5, -2), (1, 0, -1, -1), (4, -5, 0, -6)]
MATRIX_GAPS = [(-9.5, -0.5), (-4.5, -0.5), (-2.0, -1.0), (0.0, -4.0), (-3.25, -0.25), (-7.5, -1.5), (-0.5, -0.5), (0.0, 0.0), (-1.0, 0.0)]
BLOSUM = sm.load('BLOSUM62')


def rows_of(target, query, alignment):
    c = alignment.coordinates
    r1, r2 = [], []
    for k in range(c.shape[1] - 1):
        dt = c[0][k + 1] - c[0][k]; dq = c[1][k + 1] - c[1][k]
        if dt and dq:
            assert dt == dq
            r1.append(target[c[0][k]:c[0][k + 1]]); r2.append(query[c[1][k]:c[1][k + 1]])
        elif dt: r1.append(target[c[0][k]:c[0][k + 1]]); r2.append('-' * dt)
        elif dq: r1.append('-' * dq); r2.append(query[c[1][k]:c[1][k + 1]])
    return ''.join(r1), ''.join(r2)


def exhaustive(a, b, sub, go, ge, affine, flags):
    """(best, set of (row1,row2)) over every column sequence; sub(x, y) = substitution score."""
    m, n = len(a), len(b)
    fs1, fe1, fs2, fe2 = flags & 1, flags & 2, flags & 4, flags & 8
    res = {}

    def score(cols):
        L = len(cols); s = 0; lead = 0
        if cols and cols[0][1] == '-' and fs1:
            while lead < L and cols[lead][1] == '-': lead += 1
        elif cols and cols[0][0] == '-' and fs2:
            while lead < L and cols[lead][0] == '-': lead += 1
        trail = 0
        if cols and cols[-1][1] == '-' and fe1:
            while trail < L - lead and cols[L - 1 - trail][1] == '-': trail += 1
        elif cols and cols[-1][0] == '-' and fe2:
            while trail < L - lead and cols[L - 1 - trail][0] == '-': trail += 1
        last = None
        for k in range(lead, L - trail):
            x, y = cols[k]
            if x != '-' and y != '-': s += sub(x, y); last = 'M'
            elif y == '-': s += ge if (last == 'X' or not affine) else go + ge; last = 'X'
            else: s += ge if (last == 'Y' or not affine) else go + ge; last = 'Y'
        return s

    def rec(i, j, cols):
        if i == m and j == n:
            res[(''.join(c[0] for c in cols), ''.join(c[1] for c in cols))] = score(cols); return
        if i < m and j < n: rec(i + 1, j + 1, cols + [(a[i], b[j])])
        if i < m: rec(i + 1, j, cols + [(a[i], '-')])
        if j < n: rec(i, j + 1, cols + [('-', b[j])])
    rec(0, 0, [])
    best = max(res.values())
    return best, {k for k, v in res.items() if abs(v - best) < 1e-9}


def biopython(a, b, matrix, ma, mi, go, ge, affine, flags):
    al = PairwiseAligner(); al.mode = 'global'
    if matrix != '-': al.substitution_matrix = sm.load(matrix)
    else: al.match_score = ma; al.mismatch_score = mi
    al.open_gap_score = go + ge if affine else ge
    al.extend_gap_score = ge
    if flags & 1: al.left_deletion_score = 0
    if flags & 2: al.right_deletion_score = 0
    if flags & 4: al.left_insertion_score = 0
    if flags & 8: al.right_insertion_score = 0
    alns = al.align(a, b)
    score = float(alns.score)
    try:
        cnt = len(alns)
    except OverflowError:
        return score, None, None
    if cnt > CAP: return score, cnt, None
    return score, cnt, {rows_of(a, b, x) for x in alns}


def mutate(s, alpha, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice(alpha))
        else: out.append(c)
    return ''.join(out) or rnd.choice(alpha)


def gen_pair(alpha, maxlen, tiny):
    if tiny:
        return ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 5))), ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 5)))
    r = rnd.random()
    if r < 0.3:   # long gap: b = a without a block (several co-optimal gap placements inside repeats), or a with an inserted block
        a = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(6, maxlen)))
        i = rnd.randint(0, len(a) - 3); k = rnd.randint(2, min(6, len(a) - i))
        b = a[:i] + a[i + k:]
        if rnd.random() < 0.5: a, b = b, a
        return a or alpha[0], b or alpha[0]
    if r < 0.55:  # tandem repeats
        unit = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 3)))
        return unit * rnd.randint(2, 5), mutate(unit * rnd.randint(1, 5), alpha, 0.2)
    if r < 0.8:
        a = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(2, maxlen)))
        return a, mutate(a, alpha, rnd.choice([0.1, 0.3, 0.6]))
    return ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, maxlen))), ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, maxlen)))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    print('kind\tmodel\tflags\tmatrix\tmatch\tmismatch\tgapOpen\tgapExtend\tseq1\tseq2\tbiopython\tcount\talignments')
    stats = {'rows': 0, 'exh': 0, 'na': 0, 'multi': 0, 'dupgap': 0}

    def emit(kind, model, flags, matrix, ma, mi, go, ge, a, b, sub, tiny):
        affine = model == 'affine'
        bs, cnt, rows = biopython(a, b, matrix, ma, mi, go, ge, affine, flags)
        if tiny and kind == 'scalar':
            best, ex = exhaustive(a, b, sub, go, ge, affine, flags)
            assert abs(best - bs) < 1e-9, ('score', a, b, flags, model, best, bs)
            if rows is not None: assert ex == rows, ('set mismatch', a, b, flags, model, (ma, mi, go, ge), sorted(ex ^ rows)[:4])
            stats['exh'] += 1
        if rows is None: stats['na'] += 1
        elif len(rows) > 1: stats['multi'] += 1
        al = 'NA' if rows is None else ';'.join('|'.join(r) for r in sorted(rows))
        print('\t'.join([kind, model, str(flags), matrix, str(ma), str(mi), repr(float(go)), repr(float(ge)), a, b, repr(bs), str(cnt if cnt is not None else 'NA'), al]))
        stats['rows'] += 1

    for idx in range(n):
        tiny = idx % 3 == 0
        alpha = rnd.choice(['AC', 'ACG', 'ACGT', 'A'] if tiny else ['AC', 'ACGT', 'ACGT', 'AC'])
        a, b = gen_pair(alpha, 12, tiny)
        for model in ('linear', 'affine'):
            for flags in (rnd.sample(range(16), 4) + [0]):
                ma, mi, go, ge = rnd.choice(SCHEMES)
                if model == 'linear': go = 0
                s = lambda x, y, ma=ma, mi=mi: ma if x == y else mi
                emit('scalar', model, flags, '-', ma, mi, go, ge, a, b, s, tiny and len(a) <= 5 and len(b) <= 5)
    for idx in range(n):     # protein, BLOSUM62, fractional gap scores
        a, b = gen_pair(rnd.choice([PROT, 'LIVM', 'AGST', 'DE']), 9, False)
        for model in ('linear', 'affine'):
            for flags in (rnd.sample(range(16), 2) + [0]):
                go, ge = rnd.choice(MATRIX_GAPS)
                if model == 'linear': go = 0.0
                emit('matrix', model, flags, 'BLOSUM62', 0, 0, go, ge, a, b, None, False)
    print(f'# rows {stats["rows"]}, exhaustive-verified {stats["exh"]}, multi {stats["multi"]}, NA {stats["na"]}', file=sys.stderr)


main()
