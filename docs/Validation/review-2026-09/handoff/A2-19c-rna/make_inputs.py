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

# run 4 (B09 WP86): sequences for the direct probes of the two programs whose output the port has to equal exactly: CONTRAfold 2.02 (`contrafold predict IN --posteriors 0.01 OUT`) and MXSCARNA (`mxscarnamod -m -writebpp IN`):
# edge lengths 1 - 12, random RNA of 15 - 330 letters with planted stems, a few with N, DNA letters. A separate random generator (seed 86) so that inputs/sets.tsv above stays the one of runs 2 and 3.
prng = random.Random(86)
PAIRS = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
def probe_seq(n, nstem):
    s = [prng.choice('ACGU') for _ in range(n)]
    for _ in range(nstem):
        k = prng.randint(3, 8)
        if n < 2 * k + 6:
            continue
        a = prng.randint(0, n - 2 * k - 5)
        b = prng.randint(a + k + 3, n - k)
        for t in range(k):
            s[a + t] = prng.choice('ACGU')
            s[b + k - 1 - t] = PAIRS[s[a + t]]
    return ''.join(s)
probes = [probe_seq(n, 0) for n in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
for n in [15, 20, 25, 30, 40, 50, 60, 75, 90, 110, 130, 160, 200, 250, 330]:
    for rep in range(6):
        probes.append(probe_seq(n + prng.randint(0, 9), prng.randint(1, 4)))
for _ in range(8):
    s = list(probe_seq(prng.randint(30, 120), 2))
    for _ in range(prng.randint(1, 4)):
        s[prng.randrange(len(s))] = 'N'
    probes.append(''.join(s))
probes.append(probe_seq(60, 2).replace('U', 'T'))
with open('inputs/probe_seqs.tsv', 'w') as f:
    for i, s in enumerate(probes):
        f.write('p%03d\t%s\n' % (i, s))
print('probes', len(probes))
