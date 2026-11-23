#!/usr/bin/env python3
"""Writes inputs/battery.tsv (the raw `fasta34` battery) and inputs/sets.tsv (sets for `mafft --fastaparttree`, without --quiet).
Deterministic (python random.Random, fixed seeds). battery.tsv columns: id, kind (p | n), query, library ('|' separated); the query and library are written to the files exactly as `splittbfast.c`
(callfasta / getfastascoremtx) writes them (names `>+===========+<i>` + 22 blanks, 60 letters per line) and `fasta34 -z3 -m10 -p|-n -Q -H -b<n> -E<n> -d0 query library 2|6` is run.
sets.tsv columns: id, options, seqs: the 60 sets of the earlier run (sets60.tsv) and 140 more (25 - 260 sequences of 6 - 400 residues)."""
import random, shutil
P, N = 'ACDEFGHIKLMNPQRSTVWY', 'acgt'
rng = random.Random(75071)


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


def fam(alpha, n, k, r=None, ind=None):
    base = ''.join(rng.choice(alpha) for _ in range(n))
    r = r if r is not None else rng.choice([0.03, 0.1, 0.25, 0.5])
    ind = ind if ind is not None else rng.choice([0.0, 0.03, 0.1])
    return [mut(base, alpha, r, ind) for _ in range(k)]


def with_x(s, alpha_x, rate):
    return ''.join(rng.choice(alpha_x) if rng.random() < rate else c for c in s)


rows = []
k = 0
# the 20 homopolymer probes of getfastascoremtx are written by in_container.py itself (they need no randomness)
for kind in 'pn':
    alpha = P if kind == 'p' else N
    for i in range(160):
        k += 1
        nlib = rng.choice([1, 2, 3, 5, 8, 15, 30, 60])
        ln = rng.choice([6, 7, 9, 12, 20, 40, 80, 150, 300, 400])
        lib = fam(alpha, ln, nlib)
        if rng.random() < 0.3:
            lib += [''.join(rng.choice(alpha) for _ in range(rng.choice([6, 15, 60, 200]))) for _ in range(rng.randint(1, 5))]
        if rng.random() < 0.25:
            xa = 'X' if kind == 'p' else 'n'
            lib = [with_x(s, xa, rng.choice([0.02, 0.1, 0.4])) for s in lib]
        if rng.random() < 0.15:
            lib = [s + (('X' if kind == 'p' else 'n') * rng.randint(1, 10)) for s in lib]
        if rng.random() < 0.15:
            lib = [(('X' if kind == 'p' else 'n') * rng.randint(1, 10)) + s for s in lib]
        if rng.random() < 0.15:
            lib = [s[:6] if len(s) > 6 else s for s in lib]
        lib = [s if len(s) >= 6 else s + alpha[0] * (6 - len(s)) for s in lib]
        if rng.random() < 0.5:
            q = lib[rng.randrange(len(lib))]
        else:
            q = mut(lib[0], alpha, 0.2, 0.05) if rng.random() < 0.5 else ''.join(rng.choice(alpha) for _ in range(rng.choice([6, 10, 30, 100, 250])))
        if len(q) < 6:
            q = q + alpha[0] * (6 - len(q))
        if rng.random() < 0.1:
            q = alpha[0] * rng.choice([6, 12, 40])           # low complexity
            lib = [alpha[0] * rng.choice([6, 12, 40])] + lib
        rows.append('%d\t%s\t%s\t%s' % (k, kind, q, '|'.join(lib)))
open('inputs/battery.tsv', 'w').write('\n'.join(rows) + '\n')

sets = []
for line in open('inputs/sets60.tsv'):
    sets.append(line.rstrip('\n'))
k = 60
opts = ['--fastaparttree', '--fastaparttree --partsize 10', '--fastaparttree --groupsize 8', '--fastaparttree --retree 1', '--fastaparttree --reorder', '--fastaparttree --treeout', '--fastaparttree --partsize 5 --reorder']
for i in range(140):
    k += 1
    alpha = P if i % 2 == 0 else 'ACGT'
    n = rng.choice([25, 40, 60, 100, 160, 260]) if i % 4 else rng.choice([6, 10, 15])
    seqs = []
    for _ in range(rng.randint(2, 4)):
        seqs += fam(alpha, rng.choice([8, 30, 60, 120, 300]), max(2, n // 3))
    rng.shuffle(seqs)
    seqs = [s for s in seqs if len(s) >= 6] or [alpha[0] * 6, alpha[1] * 6]
    if i % 9 == 0:
        xa = 'XBZ' if alpha == P else 'N'
        seqs = [with_x(s, xa, 0.05) for s in seqs]
    sets.append('%d\t%s\t%s' % (k, opts[i % len(opts)], '|'.join(seqs)))
open('inputs/sets.tsv', 'w').write('\n'.join(sets) + '\n')
print('battery', len(rows), 'sets', len(sets))
