#!/usr/bin/env python3
"""ASSEMBLY-COMPARE-001 oracle (B10 sweep U4).

Independent transcriptions of the sourced models:
 * Mash (Ondov 2016): exact canonical ACGT-only k-mer sets, J, D = -1/k ln(2J/(1+J)), Mash boundary rules,
   containment |A&B|/|A|.
 * nucmer front end (Delcher 2002 / MUMmer 4): maximal unique matches by DEFINITION (brute force: MEM of
   length >= l whose string occurs once in the reference [, once in the query]), both query strands, then a
   literal port of MUMmer 4 mgaps (src/tigr/mgaps.cc Filter_Matches, include/mummer/mgaps.hh Cluster_each /
   Process_Cluster, Use_Extents=false).
Usage: compare_oracle.py gen N SEED > cases.json   (random + planted cases with expected outputs)
"""
import json, math, random, sys

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '$': '$', '#': '#'}


def rc(s):
    return ''.join(COMP[c] for c in reversed(s))


def norm(s, sym):
    return ''.join(c if c in 'ACGT' else sym for c in s.upper())


# ---------------- Mash ----------------
def canon_kmers(seqs, k):
    out = set()
    for s in seqs:
        s = s.upper()
        for i in range(len(s) - k + 1):
            w = s[i:i + k]
            if set(w) <= set('ACGT'):
                out.add(min(w, rc(w)))
    return out


def compare(a1, a2, k):
    k1, k2 = canon_kmers(a1, k), canon_kmers(a2, k)
    sh = len(k1 & k2)
    un = len(k1 | k2)
    j = sh / un if un else 0.0
    if sh == un:
        d = 0.0
    elif j == 0:
        d = 1.0
    else:
        d = min(1.0, -math.log(2 * j / (1 + j)) / k)
    return dict(k1=len(k1), k2=len(k2), shared=sh, c1=sh / len(k1) if k1 else 0.0,
                c2=sh / len(k2) if k2 else 0.0, jaccard=j, mash=d)


# ---------------- MUMs by definition ----------------
def count_overlapping(text, w):
    c, i = 0, text.find(w)
    while i >= 0:
        c += 1
        i = text.find(w, i + 1)
    return c


def mums(ref, q, l, both):
    out = []
    n, m = len(ref), len(q)
    for i in range(n):
        for j in range(m):
            if ref[i] != q[j]:
                continue
            if i > 0 and j > 0 and ref[i - 1] == q[j - 1]:
                continue  # not left maximal
            L = 0
            while i + L < n and j + L < m and ref[i + L] == q[j + L]:
                L += 1
            if L < l:
                continue
            w = ref[i:i + L]
            if count_overlapping(ref, w) != 1:
                continue
            if both and count_overlapping(q, w) != 1:
                continue
            out.append((i, j, L))
    return out


# ---------------- mgaps port ----------------
class M:
    def __init__(s, a, b, l):
        s.s1, s.s2, s.len = a, b, l
        s.good, s.tent = True, False
        s.score = s.frm = s.adj = 0
        s.cid = 0


def filter_matches(A):
    n = len(A)
    for i in range(n - 1):
        if not A[i].good:
            continue
        idiag = A[i].s2 - A[i].s1
        iend = A[i].s2 + A[i].len
        j = i + 1
        while j < n and A[j].s2 <= iend:
            if not A[j].good:
                j += 1
                continue
            jdiag = A[j].s2 - A[j].s1
            brk = False
            if idiag == jdiag:
                ext = A[j].len + A[j].s2 - A[i].s2
                if ext > A[i].len:
                    A[i].len = ext
                    iend = A[i].s2 + ext
                A[j].good = False
            elif A[i].s1 == A[j].s1 or A[i].s2 == A[j].s2:
                if A[i].s1 == A[j].s1:
                    olap = A[i].s2 + A[i].len - A[j].s2
                else:
                    olap = A[i].s1 + A[i].len - A[j].s1
                if A[i].len < A[j].len:
                    if olap >= A[i].len // 2:
                        A[i].good = False
                        brk = True
                elif A[j].len < A[i].len:
                    if olap >= A[j].len // 2:
                        A[j].good = False
                else:
                    if olap >= A[i].len // 2:
                        A[j].tent = True
                        if A[i].tent:
                            A[i].good = False
                            brk = True
            if brk:
                break
            j += 1
    return [m for m in A if m.good]


def find(uf, a):
    if uf[a] < 0:
        return a
    i = a
    while uf[i] > 0:
        i = uf[i]
    j = a
    while uf[j] != i:
        k = uf[j]
        uf[j] = i
        j = k
    return i


def union(uf, a, b):
    if a == b:
        return
    if uf[a] < uf[b]:
        uf[a] += uf[b]
        uf[b] = a
    else:
        uf[b] += uf[a]
        uf[a] = b


def process_cluster(A, minscore, out):
    while A:
        n = len(A)
        for i in range(n):
            A[i].score, A[i].adj, A[i].frm = A[i].len, 0, -1
            for j in range(i):
                o1 = A[j].s1 + A[j].len - A[i].s1
                o2 = A[j].s2 + A[j].len - A[i].s2
                olap = max(0, o1, o2)
                pen = olap + abs((A[i].s2 - A[i].s1) - (A[j].s2 - A[j].s1))
                if A[j].score + A[i].len - pen > A[i].score:
                    A[i].frm, A[i].score, A[i].adj = j, A[j].score + A[i].len - pen, olap
        best = 0
        for i in range(1, n):
            if A[i].score > A[best].score:
                best = i
        total, i = 0, best
        while i >= 0:
            A[i].good = True
            total += A[i].len
            i = A[i].frm
        if total >= minscore:
            out.append([(m.s1, m.s2, m.len) for m in A if m.good])
        A = [m for m in A if not m.good]


def mgaps(matches, fixed=5, maxsep=90, minscore=65, factor=0.12):
    if not matches:
        return []
    A = [M(a, b, l) for a, b, l in matches]
    A.sort(key=lambda m: (m.s2, m.s1))
    A = filter_matches(A)
    n = len(A)
    uf = [-1] * (n + 1)  # 1-based
    for i in range(1, n):
        iend = A[i - 1].s2 + A[i - 1].len
        idiag = A[i - 1].s2 - A[i - 1].s1
        for j in range(i + 1, n + 1):
            sep = A[j - 1].s2 - iend
            if sep > maxsep:
                break
            dd = abs((A[j - 1].s2 - A[j - 1].s1) - idiag)
            if dd <= max(fixed, int(factor * sep)):
                union(uf, find(uf, i), find(uf, j))
    for i in range(1, n + 1):
        A[i - 1].cid = find(uf, i)
        A[i - 1].good = False
    A.sort(key=lambda m: (m.cid, m.s2, m.s1))
    out = []
    i = 0
    while i < n:
        j = i + 1
        while j < n and A[j].cid == A[i].cid:
            j += 1
        process_cluster(A[i:j], minscore, out)
        i = j
    return out


def blocks(asm1, asm2, min_block, l, max_gap=90, factor=0.12, fixed=5, both=False):
    res = []
    qs = [(i, norm(s, '$')) for i, s in asm2]
    for rid, rs in asm1:
        ref = norm(rs, '#')
        found = []
        for qi, (qid, q) in enumerate(qs):
            n = len(q)
            for strand, qq in ((0, q), (1, rc(q))):
                if len(ref) < l or len(qq) < l:
                    continue
                ms = mums(ref, qq, l, both)
                for ch in mgaps(ms, fixed, max_gap, min_block, factor):
                    s1 = min(c[0] for c in ch); e1 = max(c[0] + c[2] for c in ch) - 1
                    s2 = min(c[1] for c in ch); e2 = max(c[1] + c[2] for c in ch) - 1
                    if strand:
                        s2, e2 = n - 1 - e2, n - 1 - s2
                    found.append((s1, s2, qi, strand, e1, e2))
        for s1, s2, qi, strand, e1, e2 in sorted(found, key=lambda f: (f[0], f[1], f[2], f[3])):
            res.append([rid, s1, e1, qs[qi][0], s2, e2, bool(strand)])
    return res


def rand_seq(r, n):
    return ''.join(r.choice('ACGT') for _ in range(n))


def mutate(r, s, rate):
    return ''.join(r.choice('ACGT') if r.random() < rate else c for c in s)


def gen(n, seed):
    r = random.Random(seed)
    cases = []
    for c in range(n):
        kind = c % 4
        base = [rand_seq(r, r.randint(60, 220)) for _ in range(r.randint(1, 2))]
        if kind == 0:   # mutated copy, same orientation
            other = [mutate(r, s, r.choice([0.0, 0.02, 0.08])) for s in base]
        elif kind == 1:  # one contig reverse complemented + insertion
            other = [rc(s) if i == 0 else s for i, s in enumerate(base)]
            other[0] = other[0][:50] + rand_seq(r, 25) + other[0][50:]
        elif kind == 2:  # rearranged halves
            s = base[0]; h = len(s) // 2
            other = [s[h:] + s[:h]] + base[1:]
        else:            # unrelated + repeats + non-ACGT
            other = [rand_seq(r, 100) + 'N' * r.randint(0, 3) + base[0][:40] * 2]
        if r.random() < 0.3:
            base[0] = base[0][:20] + 'n' + base[0][21:]
        a1 = [("r%d" % i, s) for i, s in enumerate(base)]
        a2 = [("q%d" % i, s) for i, s in enumerate(other)]
        k = r.choice([4, 6, 8, 11])
        l = r.choice([8, 10, 12, 15])
        mb = r.choice([15, 30, 65])
        cases.append(dict(a1=a1, a2=a2, k=k, l=l, minBlock=mb, both=bool(r.getrandbits(1)),
                          compare=compare([s for _, s in a1], [s for _, s in a2], k),
                          blocks=blocks(a1, a2, mb, l, both=False)))
        cases[-1]['both'] = False
    return cases


if __name__ == '__main__':
    if sys.argv[1] == 'gen':
        json.dump(gen(int(sys.argv[2]), int(sys.argv[3])), sys.stdout)
