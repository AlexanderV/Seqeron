#!/usr/bin/env python3
"""GENOMIC-COMMON-001 oracle: independent longest-common-substring / right-maximal-region reference.

For random DNA pairs (alphabets ACGT and AC, planted shared blocks incl. point-mutated blocks, tandem
repeats) this computes, WITHOUT any Seqeron code:
  * the longest common substring by (a) the O(n*m) DP over all (i, j) end cells and (b) a suffix
    automaton of s streamed with t (Gusfield 1997 sect. 7.4 / Chang-Lawler matching statistics) — the script
    asserts the two agree; tie-break = substring whose first occurrence in t ends earliest; positions are the
    first (lowest) 0-based occurrence in each sequence;
  * the right-maximal regions of FindCommonRegions: for every start i in t, the longest t[i:i+L] with
    t[i:i+L] in s (L = matching statistic ms(i), brute-force `in` test), L >= max(1, minLength), distinct
    strings in order of first start, with the first occurrence in s (str.find).
Output TSV embedded as tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/common_region_oracle.tsv:
  id  seq1  seq2  minLength  lcs(seq@pos1@pos2 | -)  regions(length:pos1:pos2,... | -; the region string is seq2[pos2:pos2+length])

Usage: python3 scripts/analysis/generate_common_region_oracle.py > common_region_oracle.tsv
"""
import random

rnd = random.Random(20261118)


def rand(n, alpha):
    return ''.join(rnd.choice(alpha) for _ in range(n))


def lcs_dp(s, t):
    n, m = len(s), len(t)
    prev = [0] * (m + 1)
    best, end_t = 0, -1
    for i in range(1, n + 1):
        cur = [0] * (m + 1)
        for j in range(1, m + 1):
            if s[i - 1] == t[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best or (cur[j] == best and j - 1 < end_t):
                    best, end_t = cur[j], j - 1
        prev = cur
    return best, end_t


class SAM:
    def __init__(self, s):
        self.len = [0]; self.link = [-1]; self.next = [{}]; self.last = 0
        for ch in s:
            cur = len(self.len)
            self.len.append(self.len[self.last] + 1); self.link.append(-1); self.next.append({})
            p = self.last
            while p != -1 and ch not in self.next[p]:
                self.next[p][ch] = cur; p = self.link[p]
            if p == -1:
                self.link[cur] = 0
            else:
                q = self.next[p][ch]
                if self.len[p] + 1 == self.len[q]:
                    self.link[cur] = q
                else:
                    cl = len(self.len)
                    self.len.append(self.len[p] + 1); self.link.append(self.link[q]); self.next.append(dict(self.next[q]))
                    while p != -1 and self.next[p].get(ch) == q:
                        self.next[p][ch] = cl; p = self.link[p]
                    self.link[q] = self.link[cur] = cl
            self.last = cur

    def stream(self, t):
        v, l, out = 0, 0, []
        for ch in t:
            while v != 0 and ch not in self.next[v]:
                v = self.link[v]; l = self.len[v]
            if ch in self.next[v]:
                v = self.next[v][ch]; l += 1
            else:
                v, l = 0, 0
            out.append(l)
        return out  # out[j] = longest common suffix ending at t[j]


def lcs(s, t):
    best, end_t = lcs_dp(s, t)
    ends = SAM(s).stream(t)
    b2 = max(ends, default=0)
    assert best == b2, (s, t, best, b2)
    if best == 0:
        return '-'
    first_end = ends.index(best)
    assert first_end == end_t, (s, t)
    sub = t[first_end - best + 1:first_end + 1]
    return f'{sub}@{s.find(sub)}@{t.find(sub)}'


def regions(s, t, mn):
    mn = max(1, mn)
    ms, l = [], 0
    for i in range(len(t)):
        l = max(l - 1, 0)
        while i + l < len(t) and t[i:i + l + 1] in s:
            l += 1
        ms.append(l)
    seen, out = set(), []
    for i, l in enumerate(ms):
        if l >= mn and t[i:i + l] not in seen:
            seen.add(t[i:i + l]); out.append(f'{l}:{s.find(t[i:i + l])}:{i}')
    return ','.join(out) if out else '-'


def case(n, m, alpha, blocks, mutate, periodic):
    s, t = rand(n, alpha), rand(m, alpha)
    if periodic:
        s = rand(rnd.randint(1, 6), alpha) * rnd.randint(2, max(2, n // 3))
    for _ in range(blocks):
        if not s:
            break
        L = rnd.randint(3, max(3, min(len(s), 250)))
        p = rnd.randint(0, max(0, len(s) - L))
        blk = s[p:p + L]
        if mutate and blk:
            k = rnd.randrange(len(blk)); blk = blk[:k] + rnd.choice(alpha) + blk[k + 1:]
        q = rnd.randint(0, len(t)); t = t[:q] + blk + t[q:]
    return s, t


rows = []
for k in range(260):  # small: exhaustive-ish edge coverage
    alpha = 'AC' if k % 3 == 0 else 'ACGT'
    rows.append(case(rnd.randint(0, 30), rnd.randint(0, 30), alpha, rnd.randint(0, 2), k % 2 == 0, k % 11 == 0))
for k in range(16):   # medium: long matches (exercise the MEM stream), repeats
    alpha = 'AC' if k % 5 == 0 else 'ACGT'
    rows.append(case(rnd.randint(300, 700), rnd.randint(300, 700), alpha, rnd.randint(1, 4), k % 3 == 0, k % 9 == 0))
rows += [('', 'ACGT'), ('ACGT', ''), ('', ''), ('AAAA', 'GGGG'), ('A', 'A'), ('ACGTACGT', 'TTACGTGG'), ('CACAGAG', 'TACATAGAT')]
print('# GENOMIC-COMMON-001 oracle, generated by scripts/analysis/generate_common_region_oracle.py (no Seqeron code)')
for idx, (s, t) in enumerate(rows):
    mn = [0, 1, 2, 3, 4, 6, 9, 15, 30][rnd.randrange(9)]
    print('\t'.join([str(idx), s or '.', t or '.', str(mn), lcs(s, t), regions(s, t, mn)]))
