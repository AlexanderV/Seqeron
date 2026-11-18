#!/usr/bin/env python3
"""ALIGN-LOCAL-001 (audit A1-7) oracle: Waterman-Eggert (1987) non-intersecting local alignments + co-optimal local alignments.

Seqeron API: SequenceAligner.LocalAlignTopN / LocalAlignAffineTopN(seq1, seq2, scoring, maxAlignments, minScore) and
LocalAlignCoOptimal / LocalAlignAffineCoOptimal(seq1, seq2, scoring, maxAlignments).

Definition (Waterman & Eggert 1987, J Mol Biol 197:723; the one EMBOSS matcher -alternatives (Huang-Miller SIM) implements):
an aligned pair is a diagonal column (a_i opposite b_j); alignment k is a maximum-score local alignment among those that contain
none of the aligned pairs of alignments 1..k-1. Tie rule of Seqeron (so that the greedy sequence is deterministic): first maximal
cell in row-major order, traceback diagonal > up > left (linear) / M > X > Y (affine), stop at the first zero predecessor.

Independent references written here from the definition (not from the C# code):
  * we_ref.we_linear / we_affine: textbook full recomputation of the Smith-Waterman matrix with the forbidden pairs after every
    round (the C# implementation recomputes only the affected cells, as the 1987 paper does) -> column `topN`.
  * brute(): exhaustive enumeration of ALL local alignments (every gapped column sequence) of tiny inputs; every reported alignment
    must score the maximum among the alignments that avoid the earlier pairs (asserted here for every tiny row, and again in the C# test).
  * EMBOSS matcher 6.6.0 `-alternatives k` (installed): its alignments are re-scored and verified round by round to be optimal for
    the definition above (tie-independent check: round k score == max of the matrix with the pairs of matcher's rounds 1..k-1
    forbidden, no pair reused) -> column `emboss` = matcher's score list (NA when matcher cannot run: non-integer gap penalties,
    custom scores).
  * Biopython 1.88 PairwiseAligner(mode='local').align(a, b): the set of all co-optimal alignments -> column `coOptimal`
    (NA when more than 60) and `biopython` = optimum.
Gap convention: gap of length k = GapOpen + k*GapExtend (linear: GapExtend per position). Biopython open_gap_score = GapOpen+GapExtend,
extend_gap_score = GapExtend. EMBOSS: -gapopen G -gapextend E = GapOpen -(G - E), GapExtend -E.

Output TSV (embedded: tests/.../TestData/Alignment/topn_oracle.tsv), columns
  kind(scalar|matrix) model(linear|affine) matrix match mismatch gapOpen gapExtend seq1 seq2 maxAlignments biopython topN emboss coOptimal
topN / coOptimal entries: "score|start1|start2|end1|end2|row1|row2" separated by ';' (0-based inclusive starts/ends of the aligned substrings).

Usage: python3 scripts/alignment/generate_topn_oracle.py [n_pairs=150] > topn_oracle.tsv     (needs /tmp scratch, matcher on PATH)
"""
import os, random, re, subprocess, sys, tempfile
from itertools import product
from Bio.Align import PairwiseAligner, substitution_matrices as sm

NEG = -10**9
rnd = random.Random(20261119)
PROT = 'ACDEFGHIKLMNPQRSTVWY'
SCALAR_SCHEMES = [(2, -3, 0, -2), (2, -3, -3, -1), (1, -1, -2, -1), (5, -4, -10, -1), (3, -2, 0, -4), (1, -2, -5, -2),
                  (5, -4, -12, -4), (5, -4, 0, -8), (5, -4, -6, -2), (1, 0, -1, -1)]
EMBOSS_DNA = [(5, -4, -12, -4), (5, -4, 0, -8), (5, -4, -6, -2), (5, -4, -4, -4), (5, -4, -15, -1)]   # EDNAFULL ACGT = +5/-4
MATRIX_GAPS = [(-9.5, -0.5), (-4.5, -0.5), (-2.0, -1.0), (0.0, -4.0), (-3.25, -0.25), (-7.5, -1.5)]
EMBOSS_PROT_GAPS = [(-2.0, -1.0), (0.0, -4.0), (-9.0, -1.0), (-12.0, -2.0), (-6.0, -2.0)]
BLOSUM = sm.load('BLOSUM62')
BTAB = {(a, b): int(BLOSUM[a][b]) for a in PROT for b in PROT}


# ---------------------------------------------------------------- independent Waterman-Eggert reference (full recompute)
def dp_linear(a, b, s, d, forb):
    m, n = len(a), len(b)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = NEG if (i, j) in forb else H[i - 1][j - 1] + s(a[i - 1], b[j - 1])
            H[i][j] = max(0, diag, H[i - 1][j] + d, H[i][j - 1] + d)
    return H


def best_cell(H, m, n):
    best, bi, bj = 0, 0, 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if H[i][j] > best: best, bi, bj = H[i][j], i, j
    return best, bi, bj


def tb_linear(a, b, H, i, j, s, d, forb):
    r1, r2, pairs = [], [], []
    ei, ej = i, j
    while i > 0 and j > 0 and H[i][j] > 0:
        if (i, j) not in forb and H[i][j] == H[i - 1][j - 1] + s(a[i - 1], b[j - 1]):
            r1.append(a[i - 1]); r2.append(b[j - 1]); pairs.append((i, j)); i -= 1; j -= 1
        elif H[i][j] == H[i - 1][j] + d:
            r1.append(a[i - 1]); r2.append('-'); i -= 1
        else:
            r1.append('-'); r2.append(b[j - 1]); j -= 1
    return (''.join(reversed(r1)), ''.join(reversed(r2)), H[ei][ej], i, j, ei - 1, ej - 1, pairs)


def dp_affine(a, b, s, o, e, forb):
    m, n = len(a), len(b); oe = o + e
    M = [[NEG] * (n + 1) for _ in range(m + 1)]; X = [[NEG] * (n + 1) for _ in range(m + 1)]; Y = [[NEG] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            M[i][j] = NEG if (i, j) in forb else max(0, M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]) + s(a[i - 1], b[j - 1])
            X[i][j] = max(NEG, M[i - 1][j] + oe, X[i - 1][j] + e, Y[i - 1][j] + oe)
            Y[i][j] = max(NEG, M[i][j - 1] + oe, Y[i][j - 1] + e, X[i][j - 1] + oe)
    return M, X, Y


def tb_affine(a, b, M, X, Y, i, j, s, o, e):
    oe = o + e
    r1, r2, pairs = [], [], []
    ei, ej = i, j; best = M[i][j]; state = 'M'

    def pick(target, cands):
        for val, st in cands:
            if val == target: return st
        raise RuntimeError('traceback')
    while True:
        if state == 'M':
            prev = M[i][j] - s(a[i - 1], b[j - 1])
            r1.append(a[i - 1]); r2.append(b[j - 1]); pairs.append((i, j)); i -= 1; j -= 1
            if prev == 0: break
            state = pick(prev, [(M[i][j], 'M'), (X[i][j], 'X'), (Y[i][j], 'Y')])
        elif state == 'X':
            cur = X[i][j]; r1.append(a[i - 1]); r2.append('-'); i -= 1
            state = pick(cur, [(M[i][j] + oe, 'M'), (X[i][j] + e, 'X'), (Y[i][j] + oe, 'Y')])
        else:
            cur = Y[i][j]; r1.append('-'); r2.append(b[j - 1]); j -= 1
            state = pick(cur, [(M[i][j] + oe, 'M'), (Y[i][j] + e, 'Y'), (X[i][j] + oe, 'X')])
    return (''.join(reversed(r1)), ''.join(reversed(r2)), best, i, j, ei - 1, ej - 1, pairs)


def we(a, b, s, o, e, k, affine):
    """Greedy Waterman-Eggert; returns [(row1, row2, scaledScore, start1, start2, end1, end2)] and the final forbidden set."""
    forb = set(); out = []
    while len(out) < k:
        if affine:
            M, X, Y = dp_affine(a, b, s, o, e, forb); H = M
        else:
            H = dp_linear(a, b, s, e, forb)
        best, bi, bj = best_cell(H, len(a), len(b))
        if best <= 0: break
        r = tb_affine(a, b, M, X, Y, bi, bj, s, o, e) if affine else tb_linear(a, b, H, bi, bj, s, e, forb)
        out.append(r[:7]); forb.update(r[7])
    return out


# ---------------------------------------------------------------- exhaustive definition check for tiny inputs
def brute_best(a, b, s, o, e, forb):
    """Maximum score over ALL local alignments (any gapped column sequence beginning and ending with an aligned pair) that use no forbidden pair."""
    best = 0
    m, n = len(a), len(b)

    def rec(i, j, score, last):   # i, j = residues consumed so far; last = last column kind
        nonlocal best
        if last == 'P' and score > best: best = score
        if i < m and j < n and (i + 1, j + 1) not in forb:
            rec(i + 1, j + 1, score + s(a[i], b[j]), 'P')
        if last != 'S':   # no leading gap: a local alignment begins with a pair
            if i < m: rec(i + 1, j, score + (e if last == 'X' else o + e), 'X')
            if j < n: rec(i, j + 1, score + (e if last == 'Y' else o + e), 'Y')
    for i0 in range(m):
        for j0 in range(n):
            if (i0 + 1, j0 + 1) not in forb: rec(i0 + 1, j0 + 1, s(a[i0], b[j0]), 'P')
    return best


# ---------------------------------------------------------------- EMBOSS matcher, verified against the definition
def matcher_scores(a, b, k, G, E, matrix, kind, s, o, e, d):
    open(d + '/a.fa', 'w').write('>a\n' + a + '\n'); open(d + '/b.fa', 'w').write('>b\n' + b + '\n')
    typ = ['-snucleotide1', '-snucleotide2'] if kind == 'dna' else ['-sprotein1', '-sprotein2']
    r = subprocess.run(['matcher', d + '/a.fa', d + '/b.fa', '-gapopen', str(G), '-gapextend', str(E), '-datafile', matrix,
                        '-alternatives', str(k), '-aformat3', 'srspair', '-outfile', d + '/o.txt', '-auto'] + typ, capture_output=True, text=True)
    if r.returncode != 0: return None
    txt = open(d + '/o.txt').read()
    forb = set(); scores = []
    for part in re.split(r'#=======================================\n#\n# Aligned_sequences', txt)[1:]:
        sc = int(re.search(r'# Score: (\d+)', part).group(1))
        ra = rb = ''; sa = sb = None
        for line in part.split('\n'):
            m_ = re.match(r'^a\s+(\d+) (\S+)\s+(\d+)$', line)
            if m_: sa = sa or int(m_.group(1)); ra += m_.group(2); continue
            m_ = re.match(r'^b\s+(\d+) (\S+)\s+(\d+)$', line)
            if m_: sb = sb or int(m_.group(1)); rb += m_.group(2)
        i, j = sa - 1, sb - 1; tot = 0; pairs = []; gap = None
        for x, y in zip(ra, rb):
            if x != '-' and y != '-': i += 1; j += 1; tot += s(x, y); pairs.append((i, j)); gap = None
            elif x == '-': j += 1; tot += e if gap == 'y' else o + e; gap = 'y'
            else: i += 1; tot += e if gap == 'x' else o + e; gap = 'x'
        M, X, Y = dp_affine(a, b, s, o, e, forb)
        opt = max(max(row) for row in M)
        assert tot == sc and opt == sc and not any(p in forb for p in pairs), ('matcher deviates from the Waterman-Eggert definition', a, b, k, G, E)
        forb.update(pairs); scores.append(sc)
    return scores


# ---------------------------------------------------------------- Biopython co-optimal set
def rows_of(target, query, alignment):
    c = alignment.coordinates
    r1, r2 = [], []
    for k in range(c.shape[1] - 1):
        dt = c[0][k + 1] - c[0][k]; dq = c[1][k + 1] - c[1][k]
        if dt and dq: r1.append(target[c[0][k]:c[0][k + 1]]); r2.append(query[c[1][k]:c[1][k + 1]])
        elif dt: r1.append(target[c[0][k]:c[0][k + 1]]); r2.append('-' * dt)
        elif dq: r1.append('-' * dq); r2.append(query[c[1][k]:c[1][k + 1]])
    return ''.join(r1), ''.join(r2)


def biopython(a, b, kind, matrix, ma, mi, go, ge, model, cap=60):
    al = PairwiseAligner(); al.mode = 'local'
    if matrix: al.substitution_matrix = sm.load(matrix)
    else: al.match_score = ma; al.mismatch_score = mi
    al.open_gap_score = ge if model == 'linear' else go + ge
    al.extend_gap_score = ge
    alns = al.align(a, b)
    try:
        score = float(alns.score)
    except Exception:
        return 0.0, 'NA'
    if score <= 0: return score, ''
    try:
        if len(alns) > cap: return score, 'NA'
    except OverflowError:
        return score, 'NA'
    out = set()
    for x in alns:
        r1, r2 = rows_of(a, b, x)
        c = x.coordinates
        out.add((r1, r2, int(c[0][0]), int(c[1][0]), int(c[0][-1]) - 1, int(c[1][-1]) - 1))
    return score, out


# ---------------------------------------------------------------- sequence generation (tandem-rich + random)
def mutate(s, alpha, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice(alpha))
        else: out.append(c)
    return ''.join(out) or rnd.choice(alpha)


def gen_pair(alpha, maxlen, tiny=False):
    if tiny:
        return ''.join(rnd.choice(alpha) for _ in range(rnd.randint(2, 5))), ''.join(rnd.choice(alpha) for _ in range(rnd.randint(2, 5)))
    r = rnd.random()
    unit = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(2, 6)))
    flank = lambda: ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, maxlen // 4)))
    if r < 0.55:    # tandem repeats on both sides -> many alignments
        a = flank() + unit * rnd.randint(2, 5) + flank(); b = flank() + mutate(unit * rnd.randint(1, 5), alpha, 0.15) + flank()
    elif r < 0.8:   # one shared domain repeated at two places in a
        dom = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(5, 12)))
        a = flank() + dom + flank() + mutate(dom, alpha, 0.2) + flank(); b = flank() + mutate(dom, alpha, 0.1) + flank()
    else:
        a = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, maxlen))); b = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, maxlen)))
    return a[:maxlen * 2] or 'A', b[:maxlen * 2] or 'A'


def fmt(alns, scale):
    return ';'.join(f'{sc / scale!r}|{s1}|{s2}|{e1}|{e2}|{r1}|{r2}' for (r1, r2, sc, s1, s2, e1, e2) in alns)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    print('kind\tmodel\tmatrix\tmatch\tmismatch\tgapOpen\tgapExtend\tseq1\tseq2\tmaxAlignments\tbiopython\ttopN\temboss\tcoOptimal')
    stats = {'rows': 0, 'matcher': 0, 'brute': 0}
    with tempfile.TemporaryDirectory() as d:
        def emit(kind, model, matrix, ma, mi, go, ge, a, b, k, scale, s, bio_matrix, emboss):
            o_s, e_s = int(round(go * scale)), int(round(ge * scale))
            res = we(a, b, s, o_s, e_s, k, model == 'affine')
            for (_, _, sc, *_r) in res: assert sc > 0
            # exhaustive definition check on tiny rows
            if len(a) <= 5 and len(b) <= 5:
                forb = set()
                for idx, (r1, r2, sc, s1, s2, e1, e2) in enumerate(res):
                    gopen = o_s if model == 'affine' else 0
                    assert brute_best(a, b, s, gopen, e_s, forb) == sc, ('WE round is not the brute-force optimum', a, b, idx)
                    i, j = s1, s2
                    for x, y in zip(r1, r2):
                        if x != '-': i += 1
                        if y != '-': j += 1
                        if x != '-' and y != '-': forb.add((i, j))
                if res:
                    gopen = o_s if model == 'affine' else 0
                    assert brute_best(a, b, s, gopen, e_s, forb) == 0 or len(res) == k, ('WE stopped early', a, b)
                stats['brute'] += 1
            bs, co = biopython(a, b, kind, bio_matrix, ma, mi, go, ge, model)
            emb = 'NA'
            if emboss:
                G, E = emboss
                ms = matcher_scores(a, b, k, G, E, 'EDNAFULL' if kind == 'dna' else 'EBLOSUM62', kind, s, o_s, e_s, d)
                if ms is not None:
                    emb = ','.join(map(str, ms)); stats['matcher'] += 1
            co_s = co if isinstance(co, str) else ';'.join(f'{bs!r}|{s1}|{s2}|{e1}|{e2}|{r1}|{r2}' for (r1, r2, s1, s2, e1, e2) in sorted(co))
            if not co_s: co_s = ''
            print('\t'.join([('matrix' if matrix != '-' else 'scalar'), model, matrix, str(ma), str(mi), repr(float(go)), repr(float(ge)), a, b, str(k),
                             repr(bs), fmt(res, scale), emb, co_s]))
            stats['rows'] += 1

        for idx in range(n):
            tiny = idx % 5 == 0
            alpha = rnd.choice(['ACGT', 'ACGT', 'AC'])
            a, b = gen_pair(alpha, 18, tiny)
            for model in ('linear', 'affine'):
                ma, mi, go, ge = rnd.choice(SCALAR_SCHEMES)
                s = lambda x, y, ma=ma, mi=mi: ma if x == y else mi
                emit('dna', model, '-', ma, mi, go, ge, a, b, rnd.randint(1, 8), 1, s, None, None)
            # EMBOSS matcher rows: EDNAFULL ACGT (+5/-4), integer gaps
            ma, mi, go, ge = rnd.choice(EMBOSS_DNA)
            model = 'linear' if go == 0 else 'affine'
            s = lambda x, y: 5 if x == y else -4
            emit('dna', model, '-', 5, -4, go, ge, a, b, rnd.randint(2, 7), 1, s, None, (int(-(go + ge)), int(-ge)))
        for idx in range(n // 2):   # protein, BLOSUM62, fractional gap scores (scale 4) and matcher rows with integer gaps
            a, b = gen_pair(rnd.choice([PROT, 'ACDEF', 'LIVM']), 25)
            for model in ('linear', 'affine'):
                go, ge = rnd.choice(MATRIX_GAPS)
                if model == 'linear': go = 0.0
                s = lambda x, y: 4 * BTAB[(x, y)]
                emit('prot', model, 'BLOSUM62', 0, 0, go, ge, a, b, rnd.randint(1, 6), 4, s, 'BLOSUM62', None)
            go, ge = rnd.choice(EMBOSS_PROT_GAPS)
            s = lambda x, y: BTAB[(x, y)]
            emit('prot', 'linear' if go == 0 else 'affine', 'BLOSUM62', 0, 0, go, ge, a, b, rnd.randint(2, 6), 1, s, 'BLOSUM62', (int(-(go + ge)), int(-ge)))
    print(f'# rows {stats["rows"]}, matcher-verified {stats["matcher"]}, brute-force-verified {stats["brute"]}', file=sys.stderr)


main()
