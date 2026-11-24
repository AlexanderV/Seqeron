#!/usr/bin/env python3
"""Writes inputs/sets.tsv: sets for `mafft --fastapair` / `--fastswpair` (dndfast7: one `fasta34 [-A] -z3 -m10 [-n] -Q -b500000 -E500000 -d500000 query library ktup` per sequence, the library = this sequence and the following ones).
Deterministic (python random.Random, fixed seed). Columns: id, options, seqs ('|' separated). 2 - 40 sequences of 6 - 300 residues, protein and DNA, related families, unknown letters (X / N), low complexity, duplicates, very short sequences, lower case."""
import random
P, N = 'ACDEFGHIKLMNPQRSTVWY', 'ACGT'
rng = random.Random(75072)


def mut(s, alpha, r, ind):
    out = []
    for c in s:
        x = rng.random()
        if x < ind / 2:
            continue
        if x < ind:
            out.append(c)
            out.append(rng.choice(alpha))
            continue
        out.append(rng.choice(alpha) if rng.random() < r else c)
    return ''.join(out) or s[:1] or alpha[0]


def family(alpha, n, k):
    base = ''.join(rng.choice(alpha) for _ in range(n))
    r = rng.choice([0.03, 0.1, 0.25, 0.5])
    ind = rng.choice([0.0, 0.03, 0.1])
    return [mut(base, alpha, r, ind) for _ in range(k)]


OPTS = ['--fastapair', '--fastswpair', '--fastapair --retree 1', '--fastapair --maxiterate 2', '--fastswpair --treeout --distout', '--fastapair --reorder', '--fastapair --anysymbol', '--fastswpair --retree 1',
        '--fastapair --op 2.0 --ep 0.2', '--fastswpair --maxiterate 3', '--fastapair --nofft', '--fastapair --treeout']
rows = []
k = 0
for i in range(180):
    k += 1
    alpha = P if i % 2 == 0 else N
    n = rng.choice([2, 3, 4, 6, 8, 12, 20, 40])
    seqs = []
    while len(seqs) < n:
        seqs += family(alpha, rng.choice([6, 9, 15, 30, 60, 120, 300]), rng.randint(1, max(1, n // 2)))
    seqs = seqs[:n]
    if rng.random() < 0.2:
        seqs[rng.randrange(n)] = seqs[rng.randrange(n)]                                  # a duplicate
    if rng.random() < 0.25:
        xa = 'X' if alpha == P else 'N'
        seqs = [''.join(xa if rng.random() < 0.05 else c for c in s) for s in seqs]
    if rng.random() < 0.1:
        seqs[rng.randrange(n)] = alpha[0] * rng.choice([6, 12, 40])                       # low complexity
    if rng.random() < 0.1:
        j = rng.randrange(n)
        seqs[j] = seqs[j].lower()
    seqs = [s if len(s) >= 2 else s + alpha[0] * 2 for s in seqs]
    rows.append('%d\t%s\t%s' % (k, OPTS[i % len(OPTS)], '|'.join(seqs)))
open('inputs/sets.tsv', 'w').write('\n'.join(rows) + '\n')
print('sets', len(rows))
