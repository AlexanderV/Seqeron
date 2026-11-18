#!/usr/bin/env python3
"""GENOMIC-REPEAT-001 oracle: repeated substrings, Gusfield maximal repeats and all longest repeats.

Reference definitions (Gusfield 1997, sec. 7.1 / 7.12.1; CMU 15-451 L10 sec. 2.1), computed TWICE per
(sequence, minLength) by two independent methods that share no code with the C# suffix-tree traversal:
  * brute force: every substring s[i:j], occurrence lists by dictionary (overlaps counted);
  * suffix array + Kasai LCP: repeated substrings = prefixes (length > parent LCP) of lcp-intervals;
    maximal repeats = lcp-intervals whose BWT-style left characters (text start = sentinel) are not all equal.
Emitted order = the order the sorted suffixes would produce: by the lexicographically smallest suffix that
starts with the repeat (a proper prefix sorts first), then increasing length. Positions ascending, 0-based.
Output TSV, embedded as tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/repeat_oracle.tsv:
  sequence  minLength  repeatCount  sha256(first 16 hex) of "seq:p,p\\n"*  maximalRepeats("seq:p,p;..."|-)  allLongest("seq:p,p;..."|-)
Usage: python3 generate_repeat_oracle.py > repeat_oracle.tsv
"""
import hashlib
import random
import sys

rnd = random.Random(20261118)


def occurrences(s, sub):
    out, i = [], s.find(sub)
    while i != -1:
        out.append(i)
        i = s.find(sub, i + 1)
    return out


def brute(s, min_len):
    n = len(s)
    occ = {}
    for i in range(n):
        for j in range(i + 1, n + 1):
            occ.setdefault(s[i:j], []).append(i)
    suffix_rank = {p: r for r, p in enumerate(sorted(range(n), key=lambda p: s[p:]))}
    rep = {k: v for k, v in occ.items() if len(v) >= 2 and len(k) >= max(1, min_len)}
    ordered = sorted(rep, key=lambda k: (min(suffix_rank[p] for p in rep[k]), len(k)))
    repeats = [(k, rep[k]) for k in ordered]
    # Gusfield maximal repeat: >= 2 occurrences, left-diverse and right-diverse (text ends count as distinct).
    maximal = []
    for k, ps in repeats:
        left = {s[p - 1] if p > 0 else '^' for p in ps}
        right = {s[p + len(k)] if p + len(k) < n else '$' for p in ps}
        if len(left) >= 2 and len(right) >= 2:
            maximal.append((k, ps))
    longest = []
    if repeats:
        top = max(len(k) for k, _ in repeats)
        longest = sorted(((k, ps) for k, ps in rep.items() if len(k) == top), key=lambda kv: kv[1][0])
    return repeats, maximal, longest


def sa_lcp(s, min_len):
    n = len(s)
    sa = sorted(range(n), key=lambda p: s[p:])
    rank = [0] * n
    for r, p in enumerate(sa):
        rank[p] = r
    lcp = [0] * n  # lcp[r] = LCP(sa[r-1], sa[r])
    h = 0
    for p in range(n):
        if rank[p] > 0:
            q = sa[rank[p] - 1]
            while p + h < n and q + h < n and s[p + h] == s[q + h]:
                h += 1
            lcp[rank[p]] = h
            if h:
                h -= 1
        else:
            h = 0
    # lcp-intervals by stack: (lcp value, left bound, right bound)
    intervals = []
    stack = [(0, 0)]
    for r in range(1, n + 1):
        cur = lcp[r] if r < n else -1
        left = r - 1
        while stack and stack[-1][0] > cur:
            val, lb = stack.pop()
            parent = max(cur, stack[-1][0]) if stack else cur
            intervals.append((val, parent, lb, r - 1))
            left = lb
        if not stack or stack[-1][0] < cur:
            stack.append((cur, left))
    repeats, maximal = [], []
    for val, parent, lb, rb in intervals:
        if val <= 0:
            continue
        positions = sorted(sa[lb:rb + 1])
        first_suffix = s[sa[lb]:]
        left_chars = {s[p - 1] if p > 0 else '^' for p in positions}
        for length in range(max(parent + 1, max(1, min_len)), val + 1):
            repeats.append((first_suffix[:length], positions, lb, length))
        if val >= max(1, min_len) and len(left_chars) >= 2:
            maximal.append((first_suffix[:val], positions, lb, val))
    repeats.sort(key=lambda t: (t[2], t[3]))
    maximal.sort(key=lambda t: (t[2], t[3]))
    return [(a, b) for a, b, _, _ in repeats], [(a, b) for a, b, _, _ in maximal]


def fmt(items):
    return ';'.join('%s:%s' % (k, ','.join(map(str, ps))) for k, ps in items) or '-'


def digest(items):
    h = hashlib.sha256()
    for k, ps in items:
        h.update(('%s:%s\n' % (k, ','.join(map(str, ps)))).encode())
    return h.hexdigest()[:16]


def rand(n, alpha):
    return ''.join(rnd.choice(alpha) for _ in range(n))


def fibonacci(steps):
    a, b = 'A', 'AC'
    for _ in range(steps):
        a, b = b, b + a
    return b


seqs = ['ATCGATCGA', 'AAAAAAAAAA', 'ATATATA', 'ACGT', 'ACGTACGT', 'ACGTACGTTTTTACGT', 'GATTACA', 'A', 'AA', 'ACA',
        'AAAAAAAAACAAAAAAAAA', 'CAAAAAAAAAAAAAAAAAAAC', 'ACACACACACAC', 'ACGTAGCTAGCTAGCTTACGT', 'GAATTCAAAGAATTCTTTA',
        fibonacci(6), fibonacci(8), ''.join('A' if bin(i).count('1') % 2 == 0 else 'C' for i in range(40)),
        'ACGT' * 6, 'ACGTTGCA' * 3, 'TTTTACGTACGTAC', 'CGCGCGAATTCGCG']
for _ in range(60):
    seqs.append(rand(rnd.randint(3, 40), 'ACGT'))
for _ in range(40):
    seqs.append(rand(rnd.randint(5, 45), 'AC'))
for _ in range(20):
    seqs.append(rand(rnd.randint(5, 40), 'ACG'))
for _ in range(15):
    unit = rand(rnd.randint(8, 14), 'ACGT')
    seqs.append(rand(rnd.randint(3, 12), 'ACGT') + unit + rand(rnd.randint(3, 15), 'ACGT') + unit + rand(rnd.randint(0, 8), 'ACGT') + unit[2:9])

rows = 0
for s in seqs:
    for min_len in (1, 2, 3, 5):
        rep_b, max_b, long_b = brute(s, min_len)
        rep_s, max_s = sa_lcp(s, min_len)
        assert rep_b == rep_s, (s, min_len)
        assert max_b == max_s, (s, min_len, max_b, max_s)
        _, _, lng = brute(s, 1)
        print('\t'.join([s, str(min_len), str(len(rep_b)), digest(rep_b), fmt(max_b), fmt(lng)]))
        rows += 1
print('rows', rows, file=sys.stderr)
