#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-10) oracle: affine-gap multiple sequence alignment.

Seqeron API: MsaOptions { GapModel = MsaGapModel.Affine } on SequenceAligner.MultipleAlign / MultipleAlignAnchored /
MultipleAlignProgressive / MultipleAlignIterative / MultipleAlignConsistency and SequenceAligner.SumOfPairsScore(rows, scoring, gapModel).

Objective (exact, Altschul 1989 "natural" gap costs): for every pair of rows delete the columns that are gaps in both rows; residue pairs
score match/mismatch; every maximal run of columns of the projection in which the same row has a gap costs open + length * extend
(gap in row i directly followed by gap in row j = two runs). Gap convention of Seqeron: gap of length k = GapOpen + k * GapExtend (negative
scores); Biopython open_gap_score = GapOpen + GapExtend, extend_gap_score = GapExtend.

Files written to OUTDIR (tests/.../TestData/Alignment/):
  msa_affine_sp_oracle.tsv     random MSAs (gapped rows): match mismatch open ext rows(;) linearSP affineSP pairBound
        linearSP / affineSP: independent Python re-score (groupby over the pair projection); pairBound = sum over pairs of the Biopython 1.88
        PairwiseAligner(mode='global') OPTIMUM of the two ungapped sequences (an upper bound of any MSA's affine SP).
  msa_affine_exact_oracle.tsv  tiny sets (k = 2, 3; length <= 7): match mismatch open ext seqs optimal optimalRows(;) numOptimalRowsChecked
        optimal = EXACT optimum of the natural affine SP by dynamic programming over (i, j, k, run state of each of the 3 pairs), validated
        against exhaustive enumeration of every alignment (all column sequences) for lengths <= 4 (asserted here).
  msa_affine_oracle.tsv        match mismatch open ext seqs star progressive iterative consistency
        (consistency: T-Coffee primary library of Gotoh global + Smith-Waterman-Gotoh local pair alignments, triplet extension, zero-gap-penalty
        progressive DP along the same guide tree)
        independent Python implementations written from the algorithm description of the C# region comment: Gotoh pair alignment
        (Flouri border, ties M > own gap > other gap), Gusfield center (best total optimal affine pair score, ties lowest index) + once-a-gap-
        always-a-gap merge, guide tree = UPGMA (exact rationals) over 1 - identities/min(len) of the Gotoh alignments, affine profile DP
        (integer Gotoh over columns with occupancy-scaled open/extend), MUSCLE-style edge-partition refinement accepted on a strict
        improvement of the natural affine SP.
Usage: python3 scripts/alignment/generate_msa_affine_oracle.py OUTDIR
"""
import itertools, os, random, sys
from fractions import Fraction
from Bio.Align import PairwiseAligner

NEG = -(1 << 60)
rnd = random.Random(20261121)


# ---------------------------------------------------------------- objective
def pair_affine(a, b, sc):
    M, X, O, E = sc
    cols = [(x, y) for x, y in zip(a, b) if not (x == '-' and y == '-')]
    score = 0
    for kind, grp in itertools.groupby(cols, key=lambda c: 'a' if c[0] == '-' else ('b' if c[1] == '-' else 'm')):
        grp = list(grp)
        if kind == 'm': score += sum(M if x == y else X for x, y in grp)
        else: score += O + E * len(grp)
    return score


def sp_affine(rows, sc):
    return sum(pair_affine(rows[i], rows[j], sc) for i, j in itertools.combinations(range(len(rows)), 2))


def sp_linear(rows, sc):
    M, X, O, E = sc; tot = 0
    for i, j in itertools.combinations(range(len(rows)), 2):
        for x, y in zip(rows[i], rows[j]):
            if x == '-' and y == '-': continue
            tot += E if (x == '-' or y == '-') else (M if x == y else X)
    return tot


def bio_optimum(a, b, sc):
    M, X, O, E = sc
    if not a or not b: return O + E * max(len(a), len(b)) if (a or b) else 0
    al = PairwiseAligner(mode='global', match_score=M, mismatch_score=X, open_gap_score=O + E, extend_gap_score=E)
    return int(round(al.score(a, b)))


# ---------------------------------------------------------------- Gotoh pair alignment (Flouri border, ties M > own > other)
def gotoh(a, b, sc):
    M, X, O, E = sc; m, n = len(a), len(b); oe = O + E
    mm = [[NEG] * (n + 1) for _ in range(m + 1)]; xx = [[NEG] * (n + 1) for _ in range(m + 1)]; yy = [[NEG] * (n + 1) for _ in range(m + 1)]
    mm[0][0] = 0
    for i in range(1, m + 1): xx[i][0] = max(NEG, mm[i-1][0] + oe, xx[i-1][0] + E, yy[i-1][0] + oe)
    for j in range(1, n + 1): yy[0][j] = max(NEG, mm[0][j-1] + oe, yy[0][j-1] + E, xx[0][j-1] + oe)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = M if a[i-1] == b[j-1] else X
            mm[i][j] = max(mm[i-1][j-1], xx[i-1][j-1], yy[i-1][j-1]) + s
            xx[i][j] = max(NEG, mm[i-1][j] + oe, xx[i-1][j] + E, yy[i-1][j] + oe)
            yy[i][j] = max(NEG, mm[i][j-1] + oe, yy[i][j-1] + E, xx[i][j-1] + oe)
            if mm[i][j] < NEG: mm[i][j] = NEG
    best = max(mm[m][n], xx[m][n], yy[m][n])
    st = 'M' if mm[m][n] == best else ('X' if xx[m][n] == best else 'Y')
    i, j = m, n; A = []; B = []
    while i > 0 or j > 0:
        if st == 'M':
            prev = mm[i][j] - (M if a[i-1] == b[j-1] else X)
            A.append(a[i-1]); B.append(b[j-1]); i -= 1; j -= 1
            st = 'M' if mm[i][j] == prev else ('X' if xx[i][j] == prev else 'Y')
        elif st == 'X':
            cur = xx[i][j]; A.append(a[i-1]); B.append('-'); i -= 1
            if i == 0 and j == 0: break
            st = 'M' if mm[i][j] + oe == cur else ('X' if xx[i][j] + E == cur else 'Y')
        else:
            cur = yy[i][j]; A.append('-'); B.append(b[j-1]); j -= 1
            if i == 0 and j == 0: break
            st = 'M' if mm[i][j] + oe == cur else ('Y' if yy[i][j] + E == cur else 'X')
    return ''.join(reversed(A)), ''.join(reversed(B)), best


# ---------------------------------------------------------------- star
def star(seqs, sc):
    k = len(seqs)
    S = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            S[i][j] = S[j][i] = gotoh(seqs[i], seqs[j], sc)[2]
    tot = [sum(S[c][j] for j in range(k) if j != c) for c in range(k)]
    c = max(range(k), key=lambda i: (tot[i], -i))
    cen = seqs[c]; al = {}; gaps = [0] * (len(cen) + 1)
    for i in range(k):
        if i == c: continue
        A, B, _ = gotoh(cen, seqs[i], sc); al[i] = (A, B)
        p = 0; g = 0
        for ch in A:
            if ch == '-': g += 1
            else: gaps[p] = max(gaps[p], g); p += 1; g = 0
        gaps[len(cen)] = max(gaps[len(cen)], g)
    al[c] = (cen, cen)
    out = []
    for i in range(k):
        A, B = al[i]; row = []; pos = 0
        for p in range(len(cen) + 1):
            t = 0
            while pos + t < len(A) and A[pos + t] == '-': t += 1
            row.append(B[pos:pos + t]); pos += t
            row.append('-' * (gaps[p] - t))
            if p < len(cen): row.append(B[pos]); pos += 1
        out.append(''.join(row))
    return out


# ---------------------------------------------------------------- progressive / iterative
def aff_dist(a, b, sc):
    if not a and not b: return Fraction(0)
    if not a or not b: return Fraction(1)
    A, B, _ = gotoh(a, b, sc)
    idn = sum(1 for x, y in zip(A, B) if x != '-' and x == y)
    return 1 - Fraction(idn, min(len(a), len(b)))


def upgma(k, D):
    nodes = {i: i for i in range(k)}; size = {i: 1 for i in range(k)}; d = {}
    for i in range(k):
        for j in range(k):
            if i != j: d[(i, j)] = D[i][j]
    active = list(range(k))
    while len(active) > 1:
        best = None
        for ii in range(len(active)):
            for jj in range(ii + 1, len(active)):
                a, b = active[ii], active[jj]; v = d[(a, b)]
                if best is None or v < best[0]: best = (v, a, b)
        _, a, b = best
        new = (nodes[a], nodes[b]); ns = size[a] + size[b]
        for c in active:
            if c in (a, b): continue
            nd = (d[(a, c)] * size[a] + d[(b, c)] * size[b]) / ns
            d[(a, c)] = nd; d[(c, a)] = nd
        nodes[a] = new; size[a] = ns; active.remove(b)
    return nodes[active[0]]


def alignprof_affine(p1, p2, sc):
    """Affine profile DP: M/X/Y over columns, occupancy-scaled open/extend (see the C# region comment)."""
    M, Xs, O, E = sc; oe = O + E
    r1, i1 = p1; r2, i2 = p2
    m, n = len(r1[0]), len(r2[0]); na, nb = len(r1), len(r2)
    res1 = [sum(1 for r in r1 if r[i] != '-') for i in range(m)]
    res2 = [sum(1 for r in r2 if r[j] != '-') for j in range(n)]

    def sub(i, j):
        t = 0
        for x in r1:
            for y in r2:
                a, b = x[i], y[j]
                if a == '-' and b == '-': continue
                t += E if (a == '-' or b == '-') else (M if a == b else Xs)
        return t
    def ea(i): return E * res1[i - 1] * nb          # A column i (1-based) against a new gap column of B
    def oea(i): return oe * res1[i - 1] * nb
    def eb(j): return E * res2[j - 1] * na
    def oeb(j): return oe * res2[j - 1] * na
    mm = [[NEG] * (n + 1) for _ in range(m + 1)]; xx = [[NEG] * (n + 1) for _ in range(m + 1)]; yy = [[NEG] * (n + 1) for _ in range(m + 1)]
    S = {}
    mm[0][0] = 0
    for i in range(0, m + 1):
        for j in range(0, n + 1):
            if i == 0 and j == 0: continue
            if i > 0 and j > 0:
                S[(i, j)] = sub(i - 1, j - 1)
                v = max(mm[i-1][j-1], xx[i-1][j-1], yy[i-1][j-1])
                mm[i][j] = NEG if v <= NEG else v + S[(i, j)]
            if i > 0:
                xx[i][j] = max(NEG, mm[i-1][j] + oea(i), xx[i-1][j] + ea(i), yy[i-1][j] + oea(i))
            if j > 0:
                yy[i][j] = max(NEG, mm[i][j-1] + oeb(j), yy[i][j-1] + eb(j), xx[i][j-1] + oeb(j))
    best = max(mm[m][n], xx[m][n], yy[m][n])
    st = 'M' if mm[m][n] == best else ('X' if xx[m][n] == best else 'Y')
    i, j = m, n; c1 = []; c2 = []
    while i > 0 or j > 0:
        if st == 'M':
            prev = mm[i][j] - S[(i, j)]
            c1.append(i - 1); c2.append(j - 1); i -= 1; j -= 1
            st = 'M' if mm[i][j] == prev else ('X' if xx[i][j] == prev else 'Y')
        elif st == 'X':
            cur = xx[i][j]; c1.append(i - 1); c2.append(-1)
            st = 'M' if mm[i-1][j] + oea(i) == cur else ('X' if xx[i-1][j] + ea(i) == cur else 'Y')
            i -= 1
        else:
            cur = yy[i][j]; c1.append(-1); c2.append(j - 1)
            st = 'M' if mm[i][j-1] + oeb(j) == cur else ('Y' if yy[i][j-1] + eb(j) == cur else 'X')
            j -= 1
    c1.reverse(); c2.reverse()
    rows = [''.join('-' if c < 0 else r[c] for c in c1) for r in r1] + [''.join('-' if c < 0 else r[c] for c in c2) for r in r2]
    return rows, i1 + i2


def build(t, seqs, sc):
    if isinstance(t, int): return [seqs[t]], [t]
    return alignprof_affine(build(t[0], seqs, sc), build(t[1], seqs, sc), sc)


def leaves(t): return [t] if isinstance(t, int) else leaves(t[0]) + leaves(t[1])


def parts(t):
    res = []
    def vis(n, depth, root):
        if isinstance(n, int): return
        if not root: res.append((depth, sorted(leaves(n))))
        vis(n[0], depth + 1, False); vis(n[1], depth + 1, False)
    vis(t, 0, True)
    res.sort(key=lambda x: (-x[0], x[1]))
    return [set(r[1]) for r in res]


def split(rows, idx, keep):
    sel = [r for r, i in enumerate(idx) if keep(i)]
    if not sel: return [], []
    cols = [c for c in range(len(rows[0])) if any(rows[r][c] != '-' for r in sel)]
    return [''.join(rows[r][c] for c in cols) for r in sel], [idx[r] for r in sel]


def progressive_iterative(seqs, sc):
    k = len(seqs)
    D = [[aff_dist(seqs[i], seqs[j], sc) if i != j else 0 for j in range(k)] for i in range(k)]
    t = upgma(k, D)
    rows, idx = build(t, seqs, sc)
    prog = [None] * k
    for r, i in zip(rows, idx): prog[i] = r
    best = sp_affine(rows, sc)
    for _ in range(16):
        imp = False
        for g in parts(t):
            a = split(rows, idx, lambda i: i in g); b = split(rows, idx, lambda i: i not in g)
            if not a[0] or not b[0]: continue
            cand, cidx = alignprof_affine(a, b, sc)
            s = sp_affine(cand, sc)
            if s > best: rows, idx, best, imp = cand, cidx, s, True
        if not imp: break
    it = [None] * k
    for r, i in zip(rows, idx): it[i] = r
    return prog, it



# ---------------------------------------------------------------- consistency (T-Coffee): affine pairwise library, zero-gap-penalty DP
def sw_affine(a, b, sc):
    """Smith-Waterman-Gotoh: end = first row-major maximum of M, traceback ties M > own gap > other gap, stops when the predecessor is 0."""
    M, Xs, O, E = sc; m, n = len(a), len(b); oe = O + E
    mm = [[NEG] * (n + 1) for _ in range(m + 1)]; xx = [[NEG] * (n + 1) for _ in range(m + 1)]; yy = [[NEG] * (n + 1) for _ in range(m + 1)]
    best = 0; bi = bj = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s_ = M if a[i-1] == b[j-1] else Xs
            mm[i][j] = max(0, mm[i-1][j-1], xx[i-1][j-1], yy[i-1][j-1]) + s_
            xx[i][j] = max(NEG, mm[i-1][j] + oe, xx[i-1][j] + E, yy[i-1][j] + oe)
            yy[i][j] = max(NEG, mm[i][j-1] + oe, yy[i][j-1] + E, xx[i][j-1] + oe)
            if mm[i][j] > best: best = mm[i][j]; bi, bj = i, j
    if best == 0: return '', '', -1, -1
    i, j = bi, bj; st = 'M'; A = []; B = []
    while True:
        if st == 'M':
            prev = mm[i][j] - (M if a[i-1] == b[j-1] else Xs)
            A.append(a[i-1]); B.append(b[j-1]); i -= 1; j -= 1
            if prev == 0: break
            st = 'M' if mm[i][j] == prev else ('X' if xx[i][j] == prev else 'Y')
        elif st == 'X':
            cur = xx[i][j]; A.append(a[i-1]); B.append('-'); i -= 1
            if i == 0 and j == 0: break
            st = 'M' if mm[i][j] + oe == cur else ('X' if xx[i][j] + E == cur else 'Y')
        else:
            cur = yy[i][j]; A.append('-'); B.append(b[j-1]); j -= 1
            if i == 0 and j == 0: break
            st = 'M' if mm[i][j] + oe == cur else ('Y' if yy[i][j] + E == cur else 'X')
    return ''.join(reversed(A)), ''.join(reversed(B)), i, j


def pid(A, B):
    al = sum(1 for x, y in zip(A, B) if x != '-' and y != '-')
    idn = sum(1 for x, y in zip(A, B) if x != '-' and x == y)
    return 0 if al == 0 else (100 * idn) // al


def add_pairs(lib, i, j, A, B, pi, pj, w):
    for x, y in zip(A, B):
        if x != '-' and y != '-':
            key = (i, pi, j, pj); lib[key] = lib.get(key, 0) + w; pi += 1; pj += 1
        elif x != '-': pi += 1
        elif y != '-': pj += 1


def consistency(seqs, sc):
    k = len(seqs); lib = {}
    for i in range(k):
        for j in range(i + 1, k):
            A, B, _ = gotoh(seqs[i], seqs[j], sc); w = pid(A, B)
            if w > 0: add_pairs(lib, i, j, A, B, 0, 0, w)
            if seqs[i] and seqs[j]:
                A, B, si, sj = sw_affine(seqs[i], seqs[j], sc); w = pid(A, B) if A else 0
                if w > 0: add_pairs(lib, i, j, A, B, si, sj, w)
    nb = {}
    for (i, p, j, q), w in lib.items():
        nb.setdefault((i, p), {})[(j, q)] = w; nb.setdefault((j, q), {})[(i, p)] = w
    def ext(s1, p1, s2, p2):
        if s1 > s2: s1, p1, s2, p2 = s2, p2, s1, p1
        t = lib.get((s1, p1, s2, p2), 0)
        for (c, r), w1 in nb.get((s1, p1), {}).items():
            if c in (s1, s2): continue
            w2 = nb.get((s2, p2), {}).get((c, r))
            if w2 is not None: t += min(w1, w2)
        return t
    D = [[aff_dist(seqs[i], seqs[j], sc) if i != j else 0 for j in range(k)] for i in range(k)]
    tree = upgma(k, D)
    def build_c(t):
        if isinstance(t, int): return [seqs[t]], [t], [list(range(len(seqs[t])))]
        r1, i1, m1 = build_c(t[0]); r2, i2, m2 = build_c(t[1])
        m = len(r1[0]); n = len(r2[0])
        def cs(c1, c2):
            tot = 0
            for a in range(len(r1)):
                pa = m1[a][c1]
                if pa < 0: continue
                for b in range(len(r2)):
                    pb = m2[b][c2]
                    if pb < 0: continue
                    tot += ext(i1[a], pa, i2[b], pb)
            return tot
        F = [[0] * (n + 1) for _ in range(m + 1)]
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                F[i][j] = max(F[i-1][j-1] + cs(i-1, j-1), F[i-1][j], F[i][j-1])
        i, j = m, n; c1 = []; c2 = []
        while i > 0 or j > 0:
            if i > 0 and j > 0 and F[i][j] == F[i-1][j-1] + cs(i-1, j-1): c1.append(i-1); c2.append(j-1); i -= 1; j -= 1
            elif i > 0 and (j == 0 or F[i][j] == F[i-1][j]): c1.append(i-1); c2.append(-1); i -= 1
            else: c1.append(-1); c2.append(j-1); j -= 1
        c1.reverse(); c2.reverse()
        rows = [''.join('-' if c < 0 else r[c] for c in c1) for r in r1] + [''.join('-' if c < 0 else r[c] for c in c2) for r in r2]
        maps = [[-1 if c < 0 else mm_[c] for c in c1] for mm_ in m1] + [[-1 if c < 0 else mm_[c] for c in c2] for mm_ in m2]
        return rows, i1 + i2, maps
    rows, idx, _ = build_c(tree)
    out = [None] * k
    for r, i in zip(rows, idx): out[i] = r
    return out

# ---------------------------------------------------------------- exact optimum (k = 2, 3)
def exact_optimum(seqs, sc):
    """DP over (i, j, l, run state of every pair); returns (optimum, rows of one optimal alignment)."""
    M, X, O, E = sc; k = len(seqs); L = [len(s) for s in seqs]
    pairs = list(itertools.combinations(range(k), 2))
    from functools import lru_cache
    import sys as _s
    _s.setrecursionlimit(100000)
    moves = [mv for mv in itertools.product((0, 1), repeat=k) if any(mv)]
    @lru_cache(maxsize=None)
    def f(pos, st):            # best score of the remaining alignment, st = run state of every pair
        if pos == tuple(L): return 0, ()
        best = None
        for mv in moves:
            npos = tuple(p + d for p, d in zip(pos, mv))
            if any(np_ > l for np_, l in zip(npos, L)): continue
            col = tuple(seqs[r][pos[r]] if mv[r] else '-' for r in range(k))
            cost = 0; nst = list(st)
            for pi, (a, b) in enumerate(pairs):
                x, y = col[a], col[b]
                if x == '-' and y == '-': continue
                if x == '-' or y == '-':
                    kind = 1 if x == '-' else 2
                    if kind != st[pi]: cost += O
                    cost += E; nst[pi] = kind
                else:
                    cost += M if x == y else X; nst[pi] = 0
            sub, tail = f(npos, tuple(nst))
            v = cost + sub
            if best is None or v > best[0]: best = (v, (col,) + tail)
        return best
    v, cols = f(tuple([0] * k), tuple([0] * len(pairs)))
    rows = [''.join(c[r] for c in cols) for r in range(k)]
    assert sp_affine(rows, sc) == v
    return v, rows


def brute_optimum(seqs, sc):
    """Exhaustive enumeration of every alignment (every sequence of columns)."""
    k = len(seqs); L = [len(s) for s in seqs]
    moves = [mv for mv in itertools.product((0, 1), repeat=k) if any(mv)]
    best = [None]
    def rec(pos, cols):
        if pos == tuple(L):
            rows = [''.join(c[r] for c in cols) for r in range(k)]
            v = sp_affine(rows, sc)
            if best[0] is None or v > best[0]: best[0] = v
            return
        for mv in moves:
            npos = tuple(p + d for p, d in zip(pos, mv))
            if any(a > l for a, l in zip(npos, L)): continue
            rec(npos, cols + [tuple(seqs[r][pos[r]] if mv[r] else '-' for r in range(k))])
    rec(tuple([0] * k), [])
    return best[0]


# ---------------------------------------------------------------- generation
def rseq(n, alpha='ACGT'): return ''.join(rnd.choice(alpha) for _ in range(n))


def mutate(s, rate):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.4: out.append(rnd.choice('ACGT'))
        elif r < rate * 0.7: continue
        elif r < rate: out.append(ch); out.extend(rnd.choice('ACGT') for _ in range(rnd.randint(1, 3)))
        else: out.append(ch)
    return ''.join(out) or 'A'


def make_set(maxlen, kmin=3, kmax=6):
    k = rnd.randint(kmin, kmax)
    base = rseq(rnd.randint(4, maxlen), rnd.choice(['ACGT', 'ACGT', 'AC', 'AAAC']))
    return [mutate(base, rnd.choice([0.05, 0.15, 0.3])) if rnd.random() < 0.7 else rseq(rnd.randint(2, maxlen)) for _ in range(k)]


SCHEMES = [(1, -1, -2, -1), (2, -3, -5, -2), (5, -4, -10, -1), (1, -2, -3, -1), (3, -2, 0, -4), (2, -1, -6, -2), (1, -1, -4, 0)]


def random_msa(k, L):
    """Random gapped rows of equal length (no all-gap column)."""
    while True:
        rows = [''.join(rnd.choice('ACGT') if rnd.random() > 0.35 else '-' for _ in range(L)) for _ in range(k)]
        if all(any(r[c] != '-' for r in rows) or rnd.random() < 0.3 for c in range(L)):   # all-gap columns allowed rarely (they are neutral)
            return rows


def main(outdir):
    # 1. objective on random MSAs
    with open(os.path.join(outdir, 'msa_affine_sp_oracle.tsv'), 'w') as f:
        f.write('match\tmismatch\topen\text\trows\tlinearSP\taffineSP\tpairBound\n')
        n = 0
        for _ in range(450):
            sc = rnd.choice(SCHEMES); k = rnd.randint(2, 6); rows = random_msa(k, rnd.randint(1, 24))
            seqs = [r.replace('-', '') for r in rows]
            bound = sum(bio_optimum(seqs[i], seqs[j], sc) for i, j in itertools.combinations(range(k), 2))
            a = sp_affine(rows, sc)
            assert a <= bound, (rows, sc, a, bound)
            f.write('\t'.join(map(str, [*sc, ';'.join(rows), sp_linear(rows, sc), a, bound])) + '\n'); n += 1
        # all-gap columns / empty rows are neutral
        for rows in (['A-C', 'A-C'], ['ACG', '---', 'AC-'], ['---', '---']):
            sc = (1, -1, -2, -1); k = len(rows); seqs = [r.replace('-', '') for r in rows]
            bound = sum(bio_optimum(seqs[i], seqs[j], sc) for i, j in itertools.combinations(range(k), 2))
            assert sp_affine(rows, sc) <= bound
            f.write('\t'.join(map(str, [*sc, ';'.join(rows), sp_linear(rows, sc), sp_affine(rows, sc), bound])) + '\n')
    print('sp rows', n + 3)

    # 2. exact optimum on tiny sets
    with open(os.path.join(outdir, 'msa_affine_exact_oracle.tsv'), 'w') as f:
        f.write('match\tmismatch\topen\text\tseqs\toptimal\toptimalRows\n')
        n = 0
        for t in range(260):
            sc = rnd.choice(SCHEMES); k = 2 if t % 8 == 0 else 3
            seqs = [rseq(rnd.randint(1, 6 if k == 3 else 8), rnd.choice(['ACGT', 'AC'])) for _ in range(k)]
            v, rows = exact_optimum(seqs, sc)
            if k == 2: assert v == bio_optimum(seqs[0], seqs[1], sc), (seqs, sc)
            if sum(map(len, seqs)) <= 9: assert v == brute_optimum(seqs, sc), (seqs, sc, v)
            f.write('\t'.join(map(str, [*sc, ';'.join(seqs), v, ';'.join(rows)])) + '\n'); n += 1
    print('exact rows', n)

    # 3. star / progressive / iterative rows
    with open(os.path.join(outdir, 'msa_affine_oracle.tsv'), 'w') as f:
        f.write('match\tmismatch\topen\text\tseqs\tstar\tprogressive\titerative\tconsistency\n')
        n = 0
        for t in range(150):
            sc = SCHEMES[t % len(SCHEMES)]
            seqs = make_set(24)
            st = star(seqs, sc)
            pr, it = progressive_iterative(seqs, sc)
            f.write('\t'.join(map(str, [*sc, ';'.join(seqs), ';'.join(st), ';'.join(pr), ';'.join(it), ';'.join(consistency(seqs, sc))])) + '\n'); n += 1
    print('msa rows', n)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
