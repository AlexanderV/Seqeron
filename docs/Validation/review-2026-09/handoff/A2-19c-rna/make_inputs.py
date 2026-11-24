#!/usr/bin/env python3
"""Writes inputs/sets.tsv: one set per line: id <TAB> option line <TAB> sequences ('|' separated). Deterministic (fixed seed). RNA families with stems (so that the base-pairing probabilities carry a signal) and plain nucleotide / protein families."""
import os, random
rng = random.Random(780019)
PAIR = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


def stem_loop(stem, loop):
    return stem + loop + ''.join(PAIR[c] for c in reversed(stem))


def family(n, length_stems):
    base = ''
    for _ in range(length_stems):
        stem = ''.join(rng.choice('ACGU') for _ in range(rng.randint(4, 7)))
        loop = ''.join(rng.choice('ACGU') for _ in range(rng.randint(3, 6)))
        base += ''.join(rng.choice('ACGU') for _ in range(rng.randint(2, 5))) + stem_loop(stem, loop)
    out = []
    for _ in range(n):
        s = []
        for c in base:
            x = rng.random()
            if x < 0.03:
                continue
            if x < 0.06:
                s.append(c + rng.choice('ACGU'))
                continue
            s.append(rng.choice('ACGU') if rng.random() < 0.08 else c)
        out.append(''.join(s))
    return out


OPTS = [
    '--mccaskill --globalpair', '--mccaskill --localpair', '--mccaskill --genafpair', '--mccaskill --globalpair --maxiterate 2', '--mccaskill --localpair --maxiterate 1000', '--mccaskill --scarnapair',
    '--scarnapair', '--scarnapair --maxiterate 2', '--rnaalifold --globalpair', '--rnaalifold --localpair', '--rnaalifold --globalpair --maxiterate 2',
    '--dafs --localpair', '--dafspair', '--contrafold --globalpair', '--contrafold --localpair', '--mccaskill --globalpair --clustalout', '--mccaskill --globalpair --retree 1', '--mccaskill --localpair --thread 1',
]
os.makedirs('inputs', exist_ok=True)
with open('inputs/sets.tsv', 'w') as f:
    k = 0
    for opts in OPTS:
        for rep in range(8):
            n = rng.choice([2, 3, 4, 5, 6, 8])
            seqs = family(n, rng.choice([1, 2, 3]))
            f.write('rna%03d\t%s\t%s\n' % (k, opts, '|'.join(seqs)))
            k += 1
    # the same options on DNA letters (T instead of U) and on a protein set: the script's own refusals
    for opts in ('--mccaskill --globalpair', '--rnaalifold --localpair', '--scarnapair'):
        seqs = [s.replace('U', 'T') for s in family(4, 2)]
        f.write('dna%03d\t%s\t%s\n' % (k, opts, '|'.join(seqs)))
        k += 1
        prot = [''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(40)) for _ in range(3)]
        f.write('prot%03d\t%s\t%s\n' % (k, opts, '|'.join(prot)))
        k += 1
print('sets', k)
