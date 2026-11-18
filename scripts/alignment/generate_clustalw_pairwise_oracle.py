#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-15) oracle: ClustalW 2.1 pairwise stage (percent identity / guide-tree distance) and the
progressive / iterative MSA built on it.

Seqeron API: SequenceAligner.ClustalWPairwiseAlign(a, b) and MultipleAlignProgressive / MultipleAlignIterative(seqs, scoring,
new MsaOptions { Distance = MsaDistanceMode.ClustalWAffine }).

References (all run here):
  * `clustalw` 2.1 (Debian clustalw 2.1+lgpl-7, installed) `-TYPE=DNA -ALIGN`: for every pair it prints
    "Sequences (i:j) Aligned. Score:  N" with N = (int)(100 * identities / MIN(len1, len2)) (single precision; FullPairwiseAlign.cpp).
    For min(len1, len2) <= 100 the integer N determines the identity count uniquely, for longer pairs it is a floor check.
  * ClustalW's pairwise stage is NOT a global alignment: forwardPass/reversePass/diff of FullPairwiseAlign.cpp (affine
    Smith-Waterman-Gotoh, end = first row-major maximum, start from the reverse pass, Myers-Miller midpoint rule for the path).
    `CW` below is a line-by-line Python port of those members (source: sources.debian.org clustalw 2.1+lgpl-7
    src/pairwise/FullPairwiseAlign.cpp, src/substitutionMatrix/SubMatrix.cpp + matrices.h). It is the generator's model of the
    binary; the binary decides (every row is asserted equal here before it is written).
  * The MSA rows come from an independent Python implementation (NW + sum-of-pairs profile DP + UPGMA with exact rationals +
    MUSCLE-style tree-partition refinement, written from the algorithm descriptions, not from the C# code) driven by the
    ClustalW distance of the port.

DNA pairwise parameters of ClustalW 2.1: IUB matrix x 100 (matrices.h `swgapdnamt`), gap open 2*15.00*100 = 3000,
gap extend 6.66*100 = 666 (a gap of length k costs 3000 + 666 k).

Output TSV files (tests/.../TestData/Alignment/):
  clustalw_pairwise_oracle.tsv  : kind(acgt|iupac|long) seq1 seq2 clustalwScore identities score
  clustalw_msa_oracle.tsv       : match mismatch gapExtend seqs(;) clustalwScores(i:j:N;...) progressive(;) iterative(;)
Usage: python3 scripts/alignment/generate_clustalw_pairwise_oracle.py OUTDIR      (needs `clustalw` on PATH)
"""
import itertools, os, random, re, subprocess, sys, tempfile
from fractions import Fraction

import numpy as np

CODES = "ABCDEFGHIKLMNOPQRSTUVWXYZ-"
ORDER = "ABCDGHKMNRSTUVWXY"
TRI = [
       10,
       -9, 10,
       -9, 10, 10,
       10, 10, -9, 10,
       -9, 10, -9, 10, 10,
       10, 10, 10, 10, -9, 10,
       -9, 10, -9, 10, 10, 10, 10,
       10, 10, 10, 10, -9, 10, -9, 10,
       10, 10, 10, 10, 10, 10, 10, 10, 10,
       10, 10, -9, 10, 10, 10, 10, 10, 10, 10,
       -9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
       -9, 10, -9, 10, -9, 10, 10, -9, 10, -9, -9, 10,
       -9, 10, -9, 10, -9, 10, 10, -9, 10, -9, -9, 10, 10,
       10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, -9, -9, 10,
       10, 10, -9, 10, -9, 10, 10, 10, 10, 10, -9, 10, 10, 10, 10,
       10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
       -9, 10, 10, 10, -9, 10, 10, 10, 10, -9, 10, 10, 10, 10, 10, 10, 10,
]
assert len(TRI) == 153


def build_matrix():
    n = len(CODES)
    M = [[0] * n for _ in range(n)]
    ix = 0
    for i, ci in enumerate(ORDER):
        for j in range(i + 1):
            k = TRI[ix] * 100
            ix += 1
            a, b = CODES.index(ci), CODES.index(ORDER[j])
            M[a][b] = k
            M[b][a] = k
    t = int(0.5 * M[0][0])
    for a, b in ((0, 4), (2, 11), (2, 12)):   # SubMatrix::getPairwiseMatrix transition-weight lines (residue codes A-E, C-M, C-N)
        M[a][b] = t
        M[b][a] = t
    return M


MAT = build_matrix()


class CW:
    """Port of FullPairwiseAlign (forwardPass, reversePass, diff, tracePath); 1-based residue indices."""

    def __init__(self, a, b, go=3000, ge=666):
        self.a = [0] + [CODES.index(c) for c in a]
        self.b = [0] + [CODES.index(c) for c in b]
        self.n, self.m, self.go, self.ge = len(a), len(b), go, ge
        L = max(self.n, self.m) + 3
        self.HH = [0] * L; self.DD = [0] * L; self.RR = [0] * L; self.SS = [0] * L
        self.displ = [0] * (2 * (self.n + self.m) + 8)
        self.printPtr = 1; self.lastPrint = 0

    def sc(self, i, j):
        return MAT[self.a[i]][self.b[j]]

    def add(self, v):
        if self.lastPrint < 0:
            self.displ[self.printPtr - 1] = v
            self.displ[self.printPtr] = self.lastPrint; self.printPtr += 1
        else:
            self.lastPrint = v; self.displ[self.printPtr] = v; self.printPtr += 1

    def delk(self, k):
        if self.lastPrint < 0:
            self.displ[self.printPtr - 1] -= k; self.lastPrint = self.displ[self.printPtr - 1]
        else:
            self.lastPrint = -k; self.displ[self.printPtr] = -k; self.printPtr += 1

    def tbgap(self, k, tb):
        return 0 if k <= 0 else tb + self.ge * k

    def forward(self):
        go, ge, HH, DD = self.go, self.ge, self.HH, self.DD
        self.maxScore = 0; self.se1 = self.se2 = 0
        for i in range(self.m + 1):
            HH[i] = 0; DD[i] = -go
        for i in range(1, self.n + 1):
            hh = p = 0; f = -go
            for j in range(1, self.m + 1):
                f -= ge
                t = hh - go - ge
                if f < t: f = t
                DD[j] -= ge
                t = HH[j] - go - ge
                if DD[j] < t: DD[j] = t
                hh = p + self.sc(i, j)
                if hh < f: hh = f
                if hh < DD[j]: hh = DD[j]
                if hh < 0: hh = 0
                p = HH[j]; HH[j] = hh
                if hh > self.maxScore:
                    self.maxScore = hh; self.se1 = i; self.se2 = j

    def reverse(self):
        go, ge, HH, DD = self.go, self.ge, self.HH, self.DD
        cost = 0; self.sb1 = self.sb2 = 1
        for i in range(self.se2, 0, -1):
            HH[i] = -1; DD[i] = -1
        for i in range(self.se1, 0, -1):
            hh = f = -1
            p = 0 if i == self.se1 else -1
            for j in range(self.se2, 0, -1):
                f -= ge
                t = hh - go - ge
                if f < t: f = t
                DD[j] -= ge
                t = HH[j] - go - ge
                if DD[j] < t: DD[j] = t
                hh = p + self.sc(i, j)
                if hh < f: hh = f
                if hh < DD[j]: hh = DD[j]
                p = HH[j]; HH[j] = hh
                if hh > cost:
                    cost = hh; self.sb1 = i; self.sb2 = j
                    if cost >= self.maxScore: break
            if cost >= self.maxScore: break

    def diff(self, A, B, M, N, tb, te):
        go, ge, HH, DD, RR, SS = self.go, self.ge, self.HH, self.DD, self.RR, self.SS
        if N <= 0:
            if M > 0: self.delk(M)
            return
        if M <= 1:
            if M <= 0:
                self.add(N); return
            midh = -(tb + ge) - self.tbgap(N, te)
            hh = -(te + ge) - self.tbgap(N, tb)
            if hh > midh: midh = hh
            midj = 0
            for j in range(1, N + 1):
                hh = self.sc(1 + A, j + B) - self.tbgap(N - j, te) - self.tbgap(j - 1, tb)
                if hh > midh: midh = hh; midj = j
            if midj == 0:
                self.delk(1); self.add(N)
            else:
                if midj > 1: self.add(midj - 1)
                self.displ[self.printPtr] = 0; self.printPtr += 1; self.lastPrint = 0
                if midj < N: self.add(N - midj)
            return
        midi = M // 2
        HH[0] = 0; t = -tb
        for j in range(1, N + 1):
            t = t - ge; HH[j] = t; DD[j] = t - go
        t = -tb
        for i in range(1, midi + 1):
            s = HH[0]; t = t - ge; HH[0] = hh = t; f = t - go
            for j in range(1, N + 1):
                hh = hh - go - ge; f = f - ge
                if hh > f: f = hh
                hh = HH[j] - go - ge; e = DD[j] - ge
                if hh > e: e = hh
                hh = s + self.sc(i + A, j + B)
                if f > hh: hh = f
                if e > hh: hh = e
                s = HH[j]; HH[j] = hh; DD[j] = e
        DD[0] = HH[0]
        RR[N] = 0; t = -te
        for j in range(N - 1, -1, -1):
            t = t - ge; RR[j] = t; SS[j] = t - go
        t = -te
        for i in range(M - 1, midi - 1, -1):
            s = RR[N]; t = t - ge; RR[N] = hh = t; f = t - go
            for j in range(N - 1, -1, -1):
                hh = hh - go - ge; f = f - ge
                if hh > f: f = hh
                hh = RR[j] - go - ge; e = SS[j] - ge
                if hh > e: e = hh
                hh = s + self.sc(i + 1 + A, j + 1 + B)
                if f > hh: hh = f
                if e > hh: hh = e
                s = RR[j]; RR[j] = hh; SS[j] = e
        SS[N] = RR[N]
        midh = HH[0] + RR[0]; midj = 0; typ = 1
        for j in range(N + 1):
            hh = HH[j] + RR[j]
            if hh >= midh and (hh > midh or (HH[j] != DD[j] and RR[j] == SS[j])):
                midh = hh; midj = j
        for j in range(N, -1, -1):
            hh = DD[j] + SS[j] + go
            if hh > midh: midh = hh; midj = j; typ = 2
        if typ == 1:
            self.diff(A, B, midi, midj, tb, go)
            self.diff(A + midi, B + midj, M - midi, N - midj, go, te)
        else:
            self.diff(A, B, midi - 1, midj, tb, 0)
            self.delk(2)
            self.diff(A + midi + 1, B + midj, M - midi - 1, N - midj, 0, te)

    def identities(self):
        i1, i2, cnt = self.sb1, self.sb2, 0
        for i in range(1, self.printPtr):
            d = self.displ[i]
            if d == 0:
                if self.a[i1] == self.b[i2]: cnt += 1
                i1 += 1; i2 += 1
            elif d > 0: i2 += d
            else: i1 -= d
        return cnt


def cw_pair(a, b):
    """-> (identities, maxScore, printed score)"""
    if not a or not b:
        return 0, 0, 0
    c = CW(a, b)
    c.forward(); c.reverse()
    c.printPtr = 1; c.lastPrint = 0
    c.diff(c.sb1 - 1, c.sb2 - 1, c.se1 - c.sb1 + 1, c.se2 - c.sb2 + 1, 0, 0)
    cnt = c.identities()
    mm = np.float32(np.float32(100.0) * np.float32(cnt)) / np.float32(min(len(a), len(b)))
    return cnt, c.maxScore, int(np.float32(mm))


def run_clustalw(seqs):
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'x.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs): f.write(f'>s{i}\n{s}\n')
        out = subprocess.run(['clustalw', f'-INFILE={fa}', '-TYPE=DNA', '-ALIGN', f'-OUTFILE={d}/o.aln'],
                             capture_output=True, text=True, cwd=d).stdout
    res = [(int(a) - 1, int(b) - 1, int(c)) for a, b, c in re.findall(r'Sequences \((\d+):(\d+)\) Aligned\. Score:\s+(-?\d+)', out)]
    assert len(res) == len(seqs) * (len(seqs) - 1) // 2, out[-400:]
    return res


rnd = random.Random(20261120)


def gen(kind, L, alpha):
    if kind == 0: return ''.join(rnd.choice(alpha) for _ in range(L))
    if kind == 1: return ''.join(rnd.choice(alpha[:2]) for _ in range(L))          # two-letter, tie heavy
    if kind == 2: return ''.join(rnd.choice([alpha[0]] * 3 + [alpha[1]]) for _ in range(L))
    u = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 4)))               # tandem repeat
    return (u * L)[:L]


def mutate(s, rate, alpha):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: continue
        elif r < rate: out.append(ch); out.append(rnd.choice(alpha))
        else: out.append(ch)
    return ''.join(out) or alpha[0]


def make_set(maxlen, alpha, kmin=3, kmax=7):
    k = rnd.randint(kmin, kmax)
    base = gen(rnd.randint(0, 3), rnd.randint(5, maxlen), alpha)
    seqs = []
    for _ in range(k):
        if rnd.random() < 0.55: seqs.append(mutate(base, rnd.choice([0.05, 0.15, 0.3]), alpha))
        else: seqs.append(gen(rnd.choice([0, 0, 1, 2, 3]), rnd.randint(3, maxlen), alpha))
    return seqs


# ---------------- independent Python MSA (written from the algorithm descriptions) ----------------
def nw(a, b, sc):
    M, X, GO, GE = sc; m, n = len(a), len(b)
    F = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1): F[i][0] = i * GE
    for j in range(n + 1): F[0][j] = j * GE
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            F[i][j] = max(F[i-1][j-1] + (M if a[i-1] == b[j-1] else X), F[i-1][j] + GE, F[i][j-1] + GE)
    i, j = m, n; A = []; B = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i-1][j-1] + (M if a[i-1] == b[j-1] else X):
            A.append(a[i-1]); B.append(b[j-1]); i -= 1; j -= 1
        elif i > 0 and F[i][j] == F[i-1][j] + GE:
            A.append(a[i-1]); B.append('-'); i -= 1
        else:
            A.append('-'); B.append(b[j-1]); j -= 1
    return ''.join(reversed(A)), ''.join(reversed(B))


def cw_dist(a, b):
    if not a or not b: return Fraction(1)
    cnt, _, _ = cw_pair(a, b)
    return 1 - Fraction(cnt, min(len(a), len(b)))


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
                a, b = active[ii], active[jj]
                v = d[(a, b)]
                if best is None or v < best[0]: best = (v, a, b)
        _, a, b = best
        new = (nodes[a], nodes[b]); ns = size[a] + size[b]
        for c in active:
            if c in (a, b): continue
            nd = (d[(a, c)] * size[a] + d[(b, c)] * size[b]) / ns
            d[(a, c)] = nd; d[(c, a)] = nd
        nodes[a] = new; size[a] = ns; active.remove(b)
    return nodes[active[0]]


def colscore(r1, c1, r2, c2, sc):
    M, X, GO, GE = sc; t = 0
    for x in r1:
        for y in r2:
            a, b = x[c1], y[c2]
            if a == '-' and b == '-': continue
            if a == '-' or b == '-': t += GE
            else: t += M if a == b else X
    return Fraction(t, len(r1) * len(r2))


def colgap(r, c, other, sc):
    M, X, GO, GE = sc
    return Fraction(sum(GE * other for x in r if x[c] != '-'), len(r) * other)


def alignprof(p1, p2, sc):
    r1, i1 = p1; r2, i2 = p2
    m = len(r1[0]); n = len(r2[0])
    F = [[Fraction(0)] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1): F[i][0] = F[i-1][0] + colgap(r1, i-1, len(r2), sc)
    for j in range(1, n + 1): F[0][j] = F[0][j-1] + colgap(r2, j-1, len(r1), sc)
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            F[i][j] = max(F[i-1][j-1] + colscore(r1, i-1, r2, j-1, sc), F[i-1][j] + colgap(r1, i-1, len(r2), sc), F[i][j-1] + colgap(r2, j-1, len(r1), sc))
    i, j = m, n; c1 = []; c2 = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i-1][j-1] + colscore(r1, i-1, r2, j-1, sc):
            c1.append(i-1); c2.append(j-1); i -= 1; j -= 1
        elif i > 0 and (j == 0 or F[i][j] == F[i-1][j] + colgap(r1, i-1, len(r2), sc)):
            c1.append(i-1); c2.append(-1); i -= 1
        else:
            c1.append(-1); c2.append(j-1); j -= 1
    c1.reverse(); c2.reverse()
    rows = [''.join('-' if c < 0 else r[c] for c in c1) for r in r1] + [''.join('-' if c < 0 else r[c] for c in c2) for r in r2]
    return rows, i1 + i2


def build(t, seqs, sc):
    if isinstance(t, int): return [seqs[t]], [t]
    return alignprof(build(t[0], seqs, sc), build(t[1], seqs, sc), sc)


def sp(rows, sc):
    M, X, GO, GE = sc; tot = 0
    for c in range(len(rows[0])):
        for i, j in itertools.combinations(range(len(rows)), 2):
            a, b = rows[i][c], rows[j][c]
            if a == '-' and b == '-': continue
            tot += GE if (a == '-' or b == '-') else (M if a == b else X)
    return tot


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


def msa(seqs, sc):
    k = len(seqs)
    D = [[cw_dist(seqs[i], seqs[j]) if i != j else 0 for j in range(k)] for i in range(k)]
    t = upgma(k, D)
    rows, idx = build(t, seqs, sc)
    prog = [None] * k
    for r, i in zip(rows, idx): prog[i] = r
    best = sp(rows, sc)
    for _ in range(16):
        imp = False
        for g in parts(t):
            a = split(rows, idx, lambda i: i in g); b = split(rows, idx, lambda i: i not in g)
            if not a[0] or not b[0]: continue
            cand, cidx = alignprof(a, b, sc)
            s = sp(cand, sc)
            if s > best: rows, idx, best, imp = cand, cidx, s, True
        if not imp: break
    it = [None] * k
    for r, i in zip(rows, idx): it[i] = r
    return prog, it


def main(outdir):
    # ---- pairwise ----
    rows = []
    bad = 0
    def emit(kind, seqs):
        nonlocal bad
        for i, j, sc in run_clustalw(seqs):
            cnt, mx, mine = cw_pair(seqs[i], seqs[j])
            assert mine == sc, (seqs[i], seqs[j], sc, mine)
            rows.append((kind, seqs[i], seqs[j], sc, cnt, mx))
    for _ in range(230): emit('acgt', make_set(45, 'ACGT'))
    for _ in range(40): emit('acgt', make_set(110, 'ACGT'))
    for _ in range(35): emit('iupac', make_set(40, 'ACGTUBDHKMNRSVWXY', 3, 6))
    for _ in range(35): emit('iupac', make_set(40, 'ACGTNRY', 3, 6))
    for _ in range(4): emit('long', make_set(320, 'ACGT', 4, 5))
    with open(os.path.join(outdir, 'clustalw_pairwise_oracle.tsv'), 'w') as f:
        f.write('kind\tseq1\tseq2\tclustalwScore\tidentities\tscore\n')
        for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
    print('pairwise rows', len(rows), 'all equal to clustalw')
    # ---- msa ----
    out = []
    schemes = [(1, -1, -1), (2, -3, -2), (5, -4, -3), (1, -2, -1)]
    for n in range(140):
        seqs = make_set(28, 'ACGT', 3, 6)
        m, x, ge = schemes[n % len(schemes)]
        sc = (m, x, 0, ge)
        scores = run_clustalw(seqs)
        prog, it = msa(seqs, sc)
        out.append((m, x, ge, ';'.join(seqs), ';'.join(f'{i}:{j}:{s}' for i, j, s in scores), ';'.join(prog), ';'.join(it)))
    with open(os.path.join(outdir, 'clustalw_msa_oracle.tsv'), 'w') as f:
        f.write('match\tmismatch\tgapExtend\tseqs\tclustalwScores\tprogressive\titerative\n')
        for r in out: f.write('\t'.join(str(x) for x in r) + '\n')
    print('msa rows', len(out))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
