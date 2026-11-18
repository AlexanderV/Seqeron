#!/usr/bin/env python3
"""GENOMIC-REPEAT-001 oracle: reverse-complement palindromes (window equals its own reverse complement).

Definition (Rosalind REVP "Locating Restriction Sites"; EMBOSS fuzznuc/Biopython `Seq.reverse_complement`): every
window s[i:i+L] with minLength <= L <= maxLength, L even, equal to its reverse complement, nested and overlapping
windows included, ordered by (0-based position, length). Computed with Biopython 1.88 `Bio.Seq.Seq.reverse_complement`
(independent of the C# centre-expansion code), cross-checked against a pure-str translate/reverse implementation.
Output TSV, embedded as tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/palindrome_oracle.tsv:
  sequence  minLength  maxLength  hits("pos:len,pos:len"|-)
Usage: python3 generate_palindrome_oracle.py > palindrome_oracle.tsv
"""
import random
import sys

from Bio.Seq import Seq

rnd = random.Random(20261118)
TR = str.maketrans('ACGT', 'TGCA')


def rc_str(w):
    return w.translate(TR)[::-1]


def hits(s, lo, hi):
    out = []
    for i in range(len(s)):
        for L in range(lo, min(hi, len(s) - i) + 1):
            if L % 2:
                continue
            w = s[i:i + L]
            bio = str(Seq(w).reverse_complement())
            assert bio == rc_str(w)
            if w == bio:
                out.append((i, L))
    return out


def rand(n, alpha):
    return ''.join(rnd.choice(alpha) for _ in range(n))


seqs = ['TCAATGCATGCGGGTCTATATGCAT', 'GAATTC', 'AAAGAATTCAAA', 'GAATTCAAAAGAATTC', 'AAAGCGGCCGCAAA', 'ACGT', 'AAAA',
        'ACGTACGTACGTACGT', 'GCGCGCGCGCGC', 'ATATATATATATATAT', 'CCCGGGAATTCCCGGG', 'A', 'AC', 'GC', 'ACGTT']
for _ in range(60):
    seqs.append(rand(rnd.randint(4, 60), 'ACGT'))
for _ in range(20):
    seqs.append(rand(rnd.randint(8, 50), 'AT'))
for _ in range(20):
    seqs.append(rand(rnd.randint(8, 50), 'CG'))
for _ in range(15):  # planted palindromes of 6-14 nt in random context
    half = rand(rnd.randint(2, 7), 'ACGT')
    pal = half + rc_str(half)
    seqs.append(rand(rnd.randint(0, 12), 'ACGT') + pal + rand(rnd.randint(0, 12), 'ACGT'))

rows = 0
for s in seqs:
    for lo, hi in ((4, 12), (4, 4), (6, 6), (4, 8), (6, 14), (8, 8), (4, 100)):
        h = hits(s, lo, hi)
        print('\t'.join([s, str(lo), str(hi), ','.join('%d:%d' % x for x in h) or '-']))
        rows += 1
print('rows', rows, file=sys.stderr)
