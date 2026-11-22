#!/usr/bin/env python3
"""Writes inputs/sets.tsv: deterministic sets (python random.Random, fixed seed) for `mafft --fastaparttree` and `--parttree` (the same sets for both: the port is checked against the fasta
variant here and against the installed program locally). Columns: id, options, seqs ('|' separated). 60 sets of 25 - 60 protein / DNA sequences of 40 - 150 residues, 4 options groups."""
import random
rng = random.Random(75052)
P, N = 'ACDEFGHIKLMNPQRSTVWY', 'ACGT'


def fam(alpha, n, k):
    base = ''.join(rng.choice(alpha) for _ in range(n))
    out = []
    for _ in range(k):
        s = []
        for c in base:
            x = rng.random()
            if x < 0.03:
                continue
            s.append(rng.choice(alpha) if x < 0.20 else c)
            if rng.random() < 0.03:
                s.append(rng.choice(alpha))
        out.append(''.join(s) or base[:1])
    return out


rows = []
k = 0
for opt in ['--fastaparttree', '--fastaparttree --partsize 10', '--fastaparttree --groupsize 8', '--fastaparttree --retree 1']:
    for i in range(15):
        k += 1
        alpha = P if i % 2 == 0 else N
        rows.append((k, opt, fam(alpha, rng.choice([40, 80, 150]), rng.randint(25, 60))))
with open('inputs/sets.tsv', 'w') as f:
    for k, opt, seqs in rows:
        f.write('%d\t%s\t%s\n' % (k, opt, '|'.join(seqs)))
print('wrote', len(rows), 'sets')
