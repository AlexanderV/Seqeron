#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-14) oracle: protein / substitution-matrix multiple sequence alignment.

Seqeron API: SequenceAligner.MultipleAlign / MultipleAlignAnchored / MultipleAlignProgressive / MultipleAlignIterative / MultipleAlignConsistency /
MultipleAlignExact (+ MultipleAlignExactWithStatistics) and MatrixSumOfPairsScore with (IEnumerable<ProteinSequence> | IEnumerable<string>,
MatrixScoring, MsaOptions) - BLOSUM62 / PAM250 / EDNAFULL residue-pair scores, real-valued (fractional) gap scores, linear or affine gaps.

Everything below is an independent Python implementation written from the algorithm descriptions (nothing of the C# code is reused);
the DNA-analogue structure (Gusfield star + once-a-gap-always-a-gap merge, UPGMA over 1 - identities / min(len) of the optimal pairwise
alignments, profile-profile dynamic program, MUSCLE-style edge-partition refinement, T-Coffee primary library + triplet extension, anchored
star = suffix-tree matching-statistics anchors + colinear chain + Needleman-Wunsch between anchors) is that of generate_msa_affine_oracle.py,
generalised to a residue-pair score sub(a, b) and an optional linear gap model. All arithmetic is exact: every score is multiplied by the
smallest integer `scale` that makes the gap scores integral (2 for 0.5, 10 for 0.1, ...); profile column averages are Fractions.

Pairwise reference: Biopython 1.88 PairwiseAligner(substitution_matrix=..., open_gap_score = GapOpen + GapExtend, extend_gap_score = GapExtend)
- the optimum of every k = 2 set (star / progressive / iterative / exact), the Carrillo-Lipman style upper bound of the SP score of any MSA
(sum over pairs of the pairwise optima) and the exhaustive enumeration of tiny sets are asserted here.

Files written to OUTDIR (tests/.../TestData/Alignment/):
  msa_matrix_sp_oracle.tsv     random gapped MSAs: matrix open ext rows(;) linearSP affineSP pairBoundLinear pairBoundAffine (exact decimal strings)
  msa_matrix_exact_oracle.tsv  model matrix open ext seqs(;) optimum optimalRows(;)   (exact DP over the k-dimensional lattice, k = 2..4; k = 2 = Biopython;
                               tiny sets also brute-forced over every column sequence)
  msa_matrix_oracle.tsv        model matrix open ext seqs(;) star anchored progressive iterative consistency (rows ;-separated, '' = skipped)
Usage: python3 scripts/alignment/generate_msa_matrix_oracle.py OUTDIR
"""
import importlib.util, itertools, os, random, sys
from fractions import Fraction
from functools import lru_cache
from math import gcd, sqrt
from multiprocessing import Pool

from Bio.Align import PairwiseAligner, substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('aff', os.path.join(HERE, 'generate_msa_affine_oracle.py'))
aff = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aff)

NEG = -(1 << 60)
sys.setrecursionlimit(1000000)

AA = 'ARNDCQEGHILKMFPSTWYV'
DNA = 'ACGT'
IUPAC = 'ACGTNRYSWKM'


# ---------------------------------------------------------------- scoring model (exact integers, scaled)
class Sc:
    def __init__(self, matrix, gap_open, gap_ext):
        self.name = matrix
        self.go, self.ge = Fraction(str(gap_open)), Fraction(str(gap_ext))          # real scores (<= 0), Seqeron convention: gap of k = go + k*ge
        self.bio = substitution_matrices.load({'EDNAFULL': 'NUC.4.4'}.get(matrix, matrix))
        den = 1
        for v in (self.go, self.ge): den = den * v.denominator // gcd(den, v.denominator)
        self.scale = den
        self.O, self.E = int(self.go * den), int(self.ge * den)           # scaled integer gap scores
        self._t = {}

    def sub(self, a, b):                                                   # scaled integer residue-pair score
        k = (a, b)
        v = self._t.get(k)
        if v is None: v = self._t[k] = int(self.bio[a][b]) * self.scale
        return v

    def real(self, scaled):
        return Fraction(scaled, self.scale)


def fmt(fr):
    """Exact decimal text of a Fraction with a terminating expansion."""
    fr = Fraction(fr)
    if fr.denominator == 1: return str(fr.numerator)
    d = fr.denominator
    while d % 2 == 0: d //= 2
    while d % 5 == 0: d //= 5
    assert d == 1, fr
    s = f'{float(fr):.12f}'.rstrip('0').rstrip('.')
    assert Fraction(s) == fr, (s, fr)
    return s


# ---------------------------------------------------------------- objectives
def pair_affine(a, b, sc):
    cols = [(x, y) for x, y in zip(a, b) if not (x == '-' and y == '-')]
    score = 0
    for kind, grp in itertools.groupby(cols, key=lambda c: 'a' if c[0] == '-' else ('b' if c[1] == '-' else 'm')):
        grp = list(grp)
        if kind == 'm': score += sum(sc.sub(x, y) for x, y in grp)
        else: score += sc.O + sc.E * len(grp)
    return score


def sp_affine(rows, sc):
    return sum(pair_affine(rows[i], rows[j], sc) for i, j in itertools.combinations(range(len(rows)), 2))


def sp_linear(rows, sc):
    tot = 0
    for i, j in itertools.combinations(range(len(rows)), 2):
        for x, y in zip(rows[i], rows[j]):
            if x == '-' and y == '-': continue
            tot += sc.E if (x == '-' or y == '-') else sc.sub(x, y)
    return tot


def sp(rows, sc, model): return sp_affine(rows, sc) if model == 'affine' else sp_linear(rows, sc)


# ---------------------------------------------------------------- Biopython pairwise optimum (real score)
def bio_aligner(sc, model, mode='global'):
    og = float(sc.go + sc.ge) if model == 'affine' else float(sc.ge)
    al = PairwiseAligner(mode=mode, substitution_matrix=sc.bio, open_gap_score=og, extend_gap_score=float(sc.ge))
    return al


def bio_optimum(a, b, sc, model):
    """Optimum of the ungapped sequences, as a scaled integer."""
    if not a or not b:
        n = max(len(a), len(b))
        return (sc.O * (1 if model == 'affine' else 0) + sc.E * n) if n else 0
    return int(round(bio_aligner(sc, model).score(a, b) * sc.scale))


# ---------------------------------------------------------------- pairwise alignments
def nw_linear(a, b, sc):
    """Needleman-Wunsch, linear gaps; traceback ties diagonal > up (gap in b) > left."""
    m, n, E = len(a), len(b), sc.E
    F = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1): F[i][0] = i * E
    for j in range(n + 1): F[0][j] = j * E
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            F[i][j] = max(F[i-1][j-1] + sc.sub(a[i-1], b[j-1]), F[i-1][j] + E, F[i][j-1] + E)
    i, j, A, B = m, n, [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i-1][j-1] + sc.sub(a[i-1], b[j-1]):
            A.append(a[i-1]); B.append(b[j-1]); i -= 1; j -= 1
        elif i > 0 and (j == 0 or F[i][j] == F[i-1][j] + E):
            A.append(a[i-1]); B.append('-'); i -= 1
        else:
            A.append('-'); B.append(b[j-1]); j -= 1
    return ''.join(reversed(A)), ''.join(reversed(B)), F[m][n]


def sw_linear(a, b, sc):
    """Smith-Waterman, linear gaps: end = first row-major maximum, traceback until the cell value is 0, ties diagonal > up > left."""
    m, n, E = len(a), len(b), sc.E
    H = [[0] * (n + 1) for _ in range(m + 1)]
    best = 0; bi = bj = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            H[i][j] = max(0, H[i-1][j-1] + sc.sub(a[i-1], b[j-1]), H[i-1][j] + E, H[i][j-1] + E)
            if H[i][j] > best: best, bi, bj = H[i][j], i, j
    if best == 0: return '', '', -1, -1
    i, j, A, B = bi, bj, [], []
    while i > 0 and j > 0 and H[i][j] > 0:
        if H[i][j] == H[i-1][j-1] + sc.sub(a[i-1], b[j-1]):
            A.append(a[i-1]); B.append(b[j-1]); i -= 1; j -= 1
        elif H[i][j] == H[i-1][j] + E:
            A.append(a[i-1]); B.append('-'); i -= 1
        else:
            A.append('-'); B.append(b[j-1]); j -= 1
    return ''.join(reversed(A)), ''.join(reversed(B)), i, j


def gotoh(a, b, sc):
    """Gotoh affine global alignment (Flouri border; ties M > own gap > other gap)."""
    E, O = sc.E, sc.O; m, n = len(a), len(b); oe = O + E
    mm = [[NEG] * (n + 1) for _ in range(m + 1)]; xx = [[NEG] * (n + 1) for _ in range(m + 1)]; yy = [[NEG] * (n + 1) for _ in range(m + 1)]
    mm[0][0] = 0
    for i in range(1, m + 1): xx[i][0] = max(NEG, mm[i-1][0] + oe, xx[i-1][0] + E, yy[i-1][0] + oe)
    for j in range(1, n + 1): yy[0][j] = max(NEG, mm[0][j-1] + oe, yy[0][j-1] + E, xx[0][j-1] + oe)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = sc.sub(a[i-1], b[j-1])
            mm[i][j] = max(mm[i-1][j-1], xx[i-1][j-1], yy[i-1][j-1]) + s
            xx[i][j] = max(NEG, mm[i-1][j] + oe, xx[i-1][j] + E, yy[i-1][j] + oe)
            yy[i][j] = max(NEG, mm[i][j-1] + oe, yy[i][j-1] + E, xx[i][j-1] + oe)
            if mm[i][j] < NEG: mm[i][j] = NEG
    best = max(mm[m][n], xx[m][n], yy[m][n])
    st = 'M' if mm[m][n] == best else ('X' if xx[m][n] == best else 'Y')
    i, j = m, n; A = []; B = []
    while i > 0 or j > 0:
        if st == 'M':
            prev = mm[i][j] - sc.sub(a[i-1], b[j-1])
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


def sw_affine(a, b, sc):
    E, O = sc.E, sc.O; m, n = len(a), len(b); oe = O + E
    mm = [[NEG] * (n + 1) for _ in range(m + 1)]; xx = [[NEG] * (n + 1) for _ in range(m + 1)]; yy = [[NEG] * (n + 1) for _ in range(m + 1)]
    best = 0; bi = bj = 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s_ = sc.sub(a[i-1], b[j-1])
            mm[i][j] = max(0, mm[i-1][j-1], xx[i-1][j-1], yy[i-1][j-1]) + s_
            xx[i][j] = max(NEG, mm[i-1][j] + oe, xx[i-1][j] + E, yy[i-1][j] + oe)
            yy[i][j] = max(NEG, mm[i][j-1] + oe, yy[i][j-1] + E, xx[i][j-1] + oe)
            if mm[i][j] > best: best = mm[i][j]; bi, bj = i, j
    if best == 0: return '', '', -1, -1
    i, j = bi, bj; st = 'M'; A = []; B = []
    while True:
        if st == 'M':
            prev = mm[i][j] - sc.sub(a[i-1], b[j-1])
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


def galign(a, b, sc, model):
    return gotoh(a, b, sc) if model == 'affine' else nw_linear(a, b, sc)


def lalign(a, b, sc, model):
    return sw_affine(a, b, sc) if model == 'affine' else sw_linear(a, b, sc)


# ---------------------------------------------------------------- star (Gusfield center + once a gap always a gap)
def merge_star(cen, al, c, k):
    gaps = [0] * (len(cen) + 1)
    for i in range(k):
        if i == c: continue
        A = al[i][0]; p = 0; g = 0
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


def star(seqs, sc, model):
    k = len(seqs)
    S = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            S[i][j] = S[j][i] = galign(seqs[i], seqs[j], sc, model)[2]
    tot = [sum(S[c][j] for j in range(k) if j != c) for c in range(k)]
    c = max(range(k), key=lambda i: (tot[i], -i))
    al = {i: galign(seqs[c], seqs[i], sc, model)[:2] for i in range(k) if i != c}
    return merge_star(seqs[c], al, c, k)


# ---------------------------------------------------------------- anchored star
def kmer_center(seqs):
    """4-mer count cosine similarity center (same float operations and summation order as the C# code)."""
    k = len(seqs)
    first = next((i for i in range(k) if seqs[i]), -1)
    if first < 0: return 0
    if k <= 2: return first
    prof = []; norms = []
    for s in seqs:
        d = {}
        for p in range(len(s) - 3): d[s[p:p+4]] = d.get(s[p:p+4], 0) + 1
        prof.append(d)
        nn = 0.0
        for v in d.values(): nn += float(v) * v
        norms.append(sqrt(nn))
    tot = [0.0] * k
    for i in range(k):
        for j in range(i + 1, k):
            small, large = (prof[i], prof[j]) if len(prof[i]) <= len(prof[j]) else (prof[j], prof[i])
            dot = 0.0
            for kk, cnt in small.items():
                if kk in large: dot += float(cnt) * large[kk]
            sim = dot / (norms[i] * norms[j]) if norms[i] > 0 and norms[j] > 0 else 0.0
            tot[i] += sim; tot[j] += sim
    best = first
    for i in range(first + 1, k):
        if not seqs[i]: continue
        if tot[i] > tot[best]: best = i
    return best


def anchors_of(ref, query, min_len=8):
    """One anchor per maximal run of query end positions with matching statistic >= min_len: the first peak of the statistic.
    Returns None when an anchor has more than one text occurrence (the C# tree reports an arbitrary one)."""
    ms = []
    for i in range(len(query)):
        L = 0
        # longest suffix of query[0..i] occurring in ref
        for l in range(min(i + 1, len(ref)), 0, -1):
            if query[i - l + 1:i + 1] in ref: L = l; break
        ms.append(L)
    res = []; peak = 0; peak_end = -1
    for i, L in enumerate(ms):
        if L >= min_len:
            if L > peak: peak, peak_end = L, i
        elif peak >= min_len:
            res.append((peak_end, peak)); peak = 0; peak_end = -1
    if peak >= min_len: res.append((peak_end, peak))
    out = []
    for end, L in res:
        q = end - L + 1; txt = query[q:q + L]
        occ = [p for p in range(len(ref) - L + 1) if ref[p:p + L] == txt]
        if len(occ) != 1: return None
        out.append((occ[0], q, L))
    return out


def chain(anchors):
    srt = sorted(anchors, key=lambda a: (a[0], a[1], -a[2]))
    n = len(srt)
    if n == 0: return []
    best = [0] * n; prev = [-1] * n; trim_of = [0] * n
    for i in range(n):
        rs, qs, L = srt[i]
        best[i] = L
        bg, bj = 0, -1
        for j in range(i):
            jr, jq, jl = srt[j]
            re_j, qe_j = jr + jl, jq + jl
            if not (re_j < rs + L and qe_j < qs + L): continue
            trim = max(0, re_j - rs, qe_j - qs)
            gain = best[j] - trim
            if gain > 0 and (gain > bg or (gain == bg and j < bj)): bg, bj = gain, j
        if bj >= 0:
            best[i] += bg; prev[i] = bj
            jr, jq, jl = srt[bj]
            trim_of[i] = max(0, jr + jl - rs, jq + jl - qs)
    bi = 0
    for i in range(1, n):
        if best[i] > best[bi]: bi = i
    out = []; idx = bi
    while idx >= 0:
        rs, qs, L = srt[idx]; t = trim_of[idx] if prev[idx] >= 0 else 0
        out.append((rs + t, qs + t, L - t)); idx = prev[idx]
    out.reverse()
    return out


def gap_fill(rg, qg, sc, model):
    if not rg and not qg: return '', ''
    if not rg: return '-' * len(qg), qg
    if not qg: return rg, '-' * len(rg)
    if rg == qg: return rg, qg
    A, B, _ = galign(rg, qg, sc, model)
    return A, B


def anchored_pair(ref, query, sc, model):
    anc = anchors_of(ref, query)
    if anc is None: return None
    ch = chain(anc)
    A, B = [], []; rp = qp = 0
    for rs, qs, L in ch:
        a, b = gap_fill(ref[rp:rs], query[qp:qs], sc, model); A.append(a); B.append(b)
        A.append(ref[rs:rs + L]); B.append(ref[rs:rs + L]); rp, qp = rs + L, qs + L
    a, b = gap_fill(ref[rp:], query[qp:], sc, model); A.append(a); B.append(b)
    return ''.join(A), ''.join(B)


def anchored(seqs, sc, model):
    k = len(seqs); c = kmer_center(seqs)
    al = {}
    for i in range(k):
        if i == c: continue
        r = anchored_pair(seqs[c], seqs[i], sc, model)
        if r is None: return None
        al[i] = r
    return merge_star(seqs[c], al, c, k)


# ---------------------------------------------------------------- progressive / iterative
def dist(a, b, sc, model):
    if not a and not b: return Fraction(0)
    if not a or not b: return Fraction(1)
    A, B, _ = galign(a, b, sc, model)
    idn = sum(1 for x, y in zip(A, B) if x != '-' and x == y)
    return 1 - Fraction(idn, min(len(a), len(b)))


def alignprof_linear(p1, p2, sc):
    """Linear profile DP: column score = average over the cross pairs (exact Fractions), traceback diagonal > up > left."""
    E = sc.E
    r1, i1 = p1; r2, i2 = p2
    m, n = len(r1[0]), len(r2[0]); na, nb = len(r1), len(r2)
    pairs = na * nb

    def col(i, j):
        t = 0
        for x in r1:
            for y in r2:
                a, b = x[i], y[j]
                if a == '-' and b == '-': continue
                t += E if (a == '-' or b == '-') else sc.sub(a, b)
        return Fraction(t, pairs)
    def g1(i): return Fraction(E * sum(1 for r in r1 if r[i] != '-') * nb, pairs)
    def g2(j): return Fraction(E * sum(1 for r in r2 if r[j] != '-') * na, pairs)
    f = [[Fraction(0)] * (n + 1) for _ in range(m + 1)]
    cs = {}
    for i in range(1, m + 1): f[i][0] = f[i-1][0] + g1(i - 1)
    for j in range(1, n + 1): f[0][j] = f[0][j-1] + g2(j - 1)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cs[(i, j)] = col(i - 1, j - 1)
            f[i][j] = max(f[i-1][j-1] + cs[(i, j)], f[i-1][j] + g1(i - 1), f[i][j-1] + g2(j - 1))
    i, j = m, n; c1 = []; c2 = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and f[i][j] == f[i-1][j-1] + cs[(i, j)]:
            c1.append(i - 1); c2.append(j - 1); i -= 1; j -= 1; continue
        if i > 0 and (j == 0 or f[i][j] == f[i-1][j] + g1(i - 1)):
            c1.append(i - 1); c2.append(-1); i -= 1; continue
        c1.append(-1); c2.append(j - 1); j -= 1
    c1.reverse(); c2.reverse()
    rows = [''.join('-' if c < 0 else r[c] for c in c1) for r in r1] + [''.join('-' if c < 0 else r[c] for c in c2) for r in r2]
    return rows, i1 + i2


def alignprof_affine(p1, p2, sc):
    """Affine profile DP: M/X/Y over columns, occupancy-scaled open/extend (the C# region comment)."""
    E, O = sc.E, sc.O; oe = O + E
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
                t += E if (a == '-' or b == '-') else sc.sub(a, b)
        return t
    def ea(i): return E * res1[i - 1] * nb
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


def alignprof(p1, p2, sc, model):
    return alignprof_affine(p1, p2, sc) if model == 'affine' else alignprof_linear(p1, p2, sc)


def build(t, seqs, sc, model):
    if isinstance(t, int): return [seqs[t]], [t]
    return alignprof(build(t[0], seqs, sc, model), build(t[1], seqs, sc, model), sc, model)


def progressive_iterative(seqs, sc, model):
    k = len(seqs)
    D = [[dist(seqs[i], seqs[j], sc, model) if i != j else 0 for j in range(k)] for i in range(k)]
    t = aff.upgma(k, D)
    rows, idx = build(t, seqs, sc, model)
    prog = [None] * k
    for r, i in zip(rows, idx): prog[i] = r
    best = sp(rows, sc, model)
    for _ in range(16):
        imp = False
        for g in aff.parts(t):
            a = aff.split(rows, idx, lambda i: i in g); b = aff.split(rows, idx, lambda i: i not in g)
            if not a[0] or not b[0]: continue
            cand, cidx = alignprof(a, b, sc, model)
            s = sp(cand, sc, model)
            if s > best: rows, idx, best, imp = cand, cidx, s, True
        if not imp: break
    it = [None] * k
    for r, i in zip(rows, idx): it[i] = r
    return prog, it


# ---------------------------------------------------------------- consistency (T-Coffee paper extension)
def consistency(seqs, sc, model):
    k = len(seqs); lib = {}
    for i in range(k):
        for j in range(i + 1, k):
            A, B, _ = galign(seqs[i], seqs[j], sc, model); w = aff.pid(A, B)
            if w > 0: aff.add_pairs(lib, i, j, A, B, 0, 0, w)
            if seqs[i] and seqs[j]:
                A, B, si, sj = lalign(seqs[i], seqs[j], sc, model); w = aff.pid(A, B) if A else 0
                if w > 0: aff.add_pairs(lib, i, j, A, B, si, sj, w)
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
    D = [[dist(seqs[i], seqs[j], sc, model) if i != j else 0 for j in range(k)] for i in range(k)]
    tree = aff.upgma(k, D)
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


# ---------------------------------------------------------------- exact optimum (k-dimensional lattice)
def exact(seqs, sc, model):
    k = len(seqs); L = tuple(len(s) for s in seqs)
    pairs = list(itertools.combinations(range(k), 2))
    moves = [mv for mv in itertools.product((0, 1), repeat=k) if any(mv)]

    @lru_cache(maxsize=None)
    def f(pos, st):
        if pos == L: return 0, ()
        best = None
        for mv in moves:
            npos = tuple(p + d for p, d in zip(pos, mv))
            if any(a > l for a, l in zip(npos, L)): continue
            col = tuple(seqs[r][pos[r]] if mv[r] else '-' for r in range(k))
            cost = 0; nst = list(st)
            for pi, (a, b) in enumerate(pairs):
                x, y = col[a], col[b]
                if x == '-' and y == '-': continue
                if x == '-' or y == '-':
                    kind = 1 if x == '-' else 2
                    if model == 'affine' and kind != st[pi]: cost += sc.O
                    cost += sc.E; nst[pi] = kind
                else:
                    cost += sc.sub(x, y); nst[pi] = 0
            if model == 'linear': nst = [0] * len(pairs)
            sub, tail = f(npos, tuple(nst))
            v = cost + sub
            if best is None or v > best[0]: best = (v, (col,) + tail)
        return best
    v, cols = f(tuple([0] * k), tuple([0] * len(pairs)))
    rows = [''.join(c[r] for c in cols) for r in range(k)]
    assert sp(rows, sc, model) == v
    return v, rows


def brute(seqs, sc, model):
    k = len(seqs); L = tuple(len(s) for s in seqs)
    moves = [mv for mv in itertools.product((0, 1), repeat=k) if any(mv)]
    best = [None]
    def rec(pos, cols):
        if pos == L:
            rows = [''.join(c[r] for c in cols) for r in range(k)]
            v = sp(rows, sc, model)
            if best[0] is None or v > best[0]: best[0] = v
            return
        for mv in moves:
            npos = tuple(p + d for p, d in zip(pos, mv))
            if any(a > l for a, l in zip(npos, L)): continue
            rec(npos, cols + [tuple(seqs[r][pos[r]] if mv[r] else '-' for r in range(k))])
    rec(tuple([0] * k), [])
    return best[0]


# ---------------------------------------------------------------- generation
SCHEMES = [  # (matrix, gap open, gap extend)  Seqeron convention: gap of length k = open + k * extend (EMBOSS gapopen G / gapextend E = (-(G - E), -E))
    ('BLOSUM62', -9.5, -0.5), ('BLOSUM62', -10, -1), ('BLOSUM62', -11, -1), ('BLOSUM62', 0, -4), ('BLOSUM62', -2.5, -0.5),
    ('PAM250', -9.5, -0.5), ('PAM250', -8, -2), ('PAM250', -4.9, -0.1), ('PAM250', 0, -8),
    ('BLOSUM62', -0.25, -0.75), ('BLOSUM62', -5, 0),
]
DNA_SCHEMES = [('EDNAFULL', -9.5, -0.5), ('EDNAFULL', -10, -1), ('EDNAFULL', -5, -2), ('EDNAFULL', 0, -3)]


def rseq(rnd, n, alpha): return ''.join(rnd.choice(alpha) for _ in range(n))


def mutate(rnd, s, rate, alpha):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.4: out.append(rnd.choice(alpha))
        elif r < rate * 0.7: continue
        elif r < rate: out.append(ch); out.extend(rnd.choice(alpha) for _ in range(rnd.randint(1, 3)))
        else: out.append(ch)
    return ''.join(out) or rnd.choice(alpha)


def make_set(rnd, maxlen, alpha, kmin=3, kmax=6):
    k = rnd.randint(kmin, kmax)
    base = rseq(rnd, rnd.randint(8, maxlen), alpha)
    return [mutate(rnd, base, rnd.choice([0.05, 0.15, 0.3]), alpha) if rnd.random() < 0.75 else rseq(rnd, rnd.randint(4, maxlen), alpha) for _ in range(k)]


def one_set(job):
    t, model, scheme, seqs = job
    sc = Sc(*scheme)
    k = len(seqs)
    st = star(seqs, sc, model)
    an = anchored(seqs, sc, model)
    pr, it = progressive_iterative(seqs, sc, model)
    co = consistency(seqs, sc, model)
    for rows in (st, pr, it, co) + ((an,) if an else ()):
        assert [r.replace('-', '') for r in rows] == seqs
    if k == 2:                                  # k = 2: every method is an optimal pairwise alignment, = Biopython
        opt = bio_optimum(seqs[0], seqs[1], sc, model)
        for rows in (st, pr, it):
            assert sp(rows, sc, model) == opt, (model, scheme, seqs, rows, opt)
    return '\t'.join([model, scheme[0], str(scheme[1]), str(scheme[2]), ';'.join(seqs), ';'.join(st), ';'.join(an) if an else '', ';'.join(pr), ';'.join(it), ';'.join(co)])


def msa_jobs(rnd):
    jobs = []
    for t in range(190):
        model = 'linear' if t % 2 == 0 else 'affine'
        scheme = SCHEMES[t % len(SCHEMES)]
        jobs.append((t, model, scheme, make_set(rnd, 36, AA)))
    for t in range(24):                         # IUPAC DNA + EDNAFULL
        jobs.append((t, 'linear' if t % 2 == 0 else 'affine', DNA_SCHEMES[t % len(DNA_SCHEMES)], make_set(rnd, 30, IUPAC)))
    for t in range(24):                         # k = 2 (Biopython optimum asserted)
        a = rseq(rnd, rnd.randint(8, 40), AA)
        jobs.append((t, 'linear' if t % 2 == 0 else 'affine', SCHEMES[t % len(SCHEMES)], [a, mutate(rnd, a, 0.2, AA) if t % 3 else rseq(rnd, rnd.randint(8, 40), AA)]))
    return jobs


def exact_job(job):
    model, scheme, seqs = job
    sc = Sc(*scheme)
    k = len(seqs)
    v, rows = exact(seqs, sc, model)
    if k == 2: assert v == bio_optimum(seqs[0], seqs[1], sc, model), (model, scheme, seqs)
    if sum(map(len, seqs)) <= (9 if k == 3 else 8): assert v == brute(seqs, sc, model), (model, scheme, seqs, v)
    return '\t'.join([model, scheme[0], str(scheme[1]), str(scheme[2]), ';'.join(seqs), fmt(sc.real(v)), ';'.join(rows)])


def exact_jobs(rnd):
    jobs = []
    for t in range(170):                        # k = 3, lengths <= 9
        model = 'linear' if t % 2 == 0 else 'affine'
        base = rseq(rnd, rnd.randint(3, 8), AA)
        seqs = [mutate(rnd, base, rnd.choice([0.1, 0.3, 0.5]), AA)[:9] for _ in range(3)] if rnd.random() < 0.6 else [rseq(rnd, rnd.randint(1, 8), AA) for _ in range(3)]
        jobs.append((model, SCHEMES[t % len(SCHEMES)], seqs))
    for t in range(40):                         # k = 2 (Biopython)
        jobs.append(('linear' if t % 2 else 'affine', SCHEMES[t % len(SCHEMES)], [rseq(rnd, rnd.randint(0, 14), AA) for _ in range(2)]))
    for t in range(50):                         # k = 4, tiny (brute force)
        jobs.append(('linear' if t % 2 else 'affine', SCHEMES[t % len(SCHEMES)], [rseq(rnd, rnd.randint(1, 3), AA) for _ in range(4)]))
    for t in range(24):                         # k = 3 DNA / EDNAFULL
        jobs.append(('linear' if t % 2 else 'affine', DNA_SCHEMES[t % len(DNA_SCHEMES)], [rseq(rnd, rnd.randint(1, 7), IUPAC) for _ in range(3)]))
    for t in range(24):                         # k = 3, longer related (pruning)
        base = rseq(rnd, rnd.randint(14, 20), AA)
        jobs.append(('linear' if t % 2 else 'affine', SCHEMES[t % len(SCHEMES)], [mutate(rnd, base, 0.15, AA)[:22] for _ in range(3)]))
    return jobs


def random_msa(rnd, k, L, alpha):
    while True:
        rows = [''.join(rnd.choice(alpha) if rnd.random() > 0.35 else '-' for _ in range(L)) for _ in range(k)]
        if all(any(r[c] != '-' for r in rows) or rnd.random() < 0.3 for c in range(L)): return rows


def sp_rows(rnd):
    out = []
    for _ in range(420):
        scheme = rnd.choice(SCHEMES + DNA_SCHEMES); sc = Sc(*scheme)
        alpha = IUPAC if scheme[0] == 'EDNAFULL' else AA
        k = rnd.randint(2, 6); rows = random_msa(rnd, k, rnd.randint(1, 24), alpha)
        seqs = [r.replace('-', '') for r in rows]
        bl = sum(bio_optimum(seqs[i], seqs[j], sc, 'linear') for i, j in itertools.combinations(range(k), 2))
        ba = sum(bio_optimum(seqs[i], seqs[j], sc, 'affine') for i, j in itertools.combinations(range(k), 2))
        a, l = sp_affine(rows, sc), sp_linear(rows, sc)
        assert a <= ba and l <= bl, (rows, scheme, a, ba, l, bl)
        out.append('\t'.join([scheme[0], str(scheme[1]), str(scheme[2]), ';'.join(rows), fmt(sc.real(l)), fmt(sc.real(a)), fmt(sc.real(bl)), fmt(sc.real(ba))]))
    return out


def main(outdir):
    rnd = random.Random(20261123)
    rows = sp_rows(rnd)
    with open(os.path.join(outdir, 'msa_matrix_sp_oracle.tsv'), 'w') as f:
        f.write('matrix\topen\text\trows\tlinearSP\taffineSP\tpairBoundLinear\tpairBoundAffine\n' + '\n'.join(rows) + '\n')
    print('sp rows', len(rows))

    with Pool(os.cpu_count()) as pool:
        ex = pool.map(exact_job, exact_jobs(rnd), chunksize=4)
        with open(os.path.join(outdir, 'msa_matrix_exact_oracle.tsv'), 'w') as f:
            f.write('model\tmatrix\topen\text\tseqs\toptimum\toptimalRows\n' + '\n'.join(ex) + '\n')
        print('exact rows', len(ex))

        jobs = msa_jobs(rnd)
        res = pool.map(one_set, jobs, chunksize=2)
    with open(os.path.join(outdir, 'msa_matrix_oracle.tsv'), 'w') as f:
        f.write('model\tmatrix\topen\text\tseqs\tstar\tanchored\tprogressive\titerative\tconsistency\n' + '\n'.join(res) + '\n')
    print('msa rows', len(res), 'anchored present:', sum(1 for r in res if r.split('\t')[6]))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
