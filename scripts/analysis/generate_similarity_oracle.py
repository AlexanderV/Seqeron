#!/usr/bin/env python3
"""GENOMIC-SIMILARITY-001 oracle: exact k-mer Jaccard index of two DNA sequences (percent), four independent ways.

  (1) Python `set` over literal (non-canonical) k-mers: |A&B| / |A|B|.
  (2) scipy.spatial.distance.jaccard on boolean presence vectors over the union vocabulary: J = 1 - d (scipy 1.18.1).
  (3) Canonical k-mers (min(w, revcomp(w)), Mash `-n`-less default) with Python `set` + Biopython Seq.reverse_complement.
  (4) sourmash 4.9.4 MinHash(n=0, ksize=k, scaled=1).jaccard (exact, canonical DNA hashes) == (3).
Both-empty k-mer sets: Jaccard undefined, convention 0 (same as KmerAnalyzer.JaccardSimilarity).
Output TSV, embedded as tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/similarity_oracle.tsv:
  seq1  seq2  k  jaccard_percent_literal  jaccard_percent_canonical
Usage: python3 generate_similarity_oracle.py > similarity_oracle.tsv
"""
import random
import sys

import numpy as np
from Bio.Seq import Seq
from scipy.spatial.distance import jaccard as scipy_jaccard
import sourmash

rnd = random.Random(20261118)


def kmers(s, k, canonical):
    out = set()
    for i in range(len(s) - k + 1):
        w = s[i:i + k]
        if canonical:
            rc = str(Seq(w).reverse_complement())
            w = min(w, rc)
        out.add(w)
    return out


def pct_set(a, b):
    u = a | b
    return 0.0 if not u else 100.0 * len(a & b) / len(u)


def pct_scipy(a, b):
    vocab = sorted(a | b)
    if not vocab:
        return 0.0
    va = np.array([v in a for v in vocab])
    vb = np.array([v in b for v in vocab])
    return 100.0 * (1.0 - scipy_jaccard(va, vb))


def pct_sourmash(s1, s2, k):
    m1 = sourmash.MinHash(n=0, ksize=k, scaled=1)
    m2 = sourmash.MinHash(n=0, ksize=k, scaled=1)
    m1.add_sequence(s1, force=True)
    m2.add_sequence(s2, force=True)
    if len(m1) == 0 and len(m2) == 0:
        return 0.0
    return 100.0 * m1.jaccard(m2)


def rand(n, alpha='ACGT'):
    return ''.join(rnd.choice(alpha) for _ in range(n))


def mutate(s, rate):
    return ''.join(rnd.choice('ACGT') if rnd.random() < rate else c for c in s)


def rc(s):
    return str(Seq(s).reverse_complement())


pairs = []
for _ in range(60):
    a = rand(rnd.randint(0, 40))
    pairs.append((a, rand(rnd.randint(0, 40))))
    pairs.append((a, mutate(a, 0.1)))
    pairs.append((a, rc(a)))                       # same strand-agnostic content
    pairs.append((a, mutate(rc(a), 0.1)))
for alpha in ('AC', 'AT', 'A'):
    for _ in range(15):
        pairs.append((rand(rnd.randint(1, 30), alpha), rand(rnd.randint(1, 30), alpha)))
pairs += [('ACGTACGT', 'ACGTACGA'), ('ACGT', 'ACGA'), ('AAAAA', 'CCCCC'), ('AAAAAA', 'AAAA'), ('', ''), ('AC', 'GT'),
          ('ACGTAC', ''), ('AAAA', 'TTTT'), ('ATGCATGC', 'ATGCATTT'), ('ACGT', 'ACGT')]
big = rand(2000)
pairs += [(big, mutate(big, 0.02)), (big, rc(big)), (big, rand(2000))]

n = 0
for s1, s2 in pairs:
    for k in (1, 2, 3, 4, 5, 8, 21):
        A, B = kmers(s1, k, False), kmers(s2, k, False)
        Ac, Bc = kmers(s1, k, True), kmers(s2, k, True)
        lit, can = pct_set(A, B), pct_set(Ac, Bc)
        assert abs(lit - pct_scipy(A, B)) < 1e-9, (s1, s2, k)
        assert abs(can - pct_scipy(Ac, Bc)) < 1e-9, (s1, s2, k)
        assert abs(can - pct_sourmash(s1, s2, k)) < 1e-9, (s1, s2, k, can, pct_sourmash(s1, s2, k))
        print(f"{s1 or '-'}\t{s2 or '-'}\t{k}\t{lit!r}\t{can!r}")
        n += 1
print(f"rows={n} scipy={__import__('scipy').__version__} sourmash={__import__('importlib.metadata').metadata.version('sourmash')}", file=sys.stderr)
