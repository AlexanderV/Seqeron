#!/usr/bin/env python3
"""Writes inputs/sets.tsv for the second run of `mafft --fastapair` / `--fastswpair` (B09 WP75): the first run (handoff/A2-19c-fastapair, 180 sets) had 63 DNA sets that end in a segmentation fault of dndfast7 and 27 that stop with "count (fasta34 -z 3)":
only ONE DNA set ran to the end, so the success path of the nucleotide reader (ReadFasta34m10_nuc: blocks of the reverse strand, several blocks per pair, cutal on both strands, -A on DNA) has one example. This set gives it more: DNA families of 150 - 1500 nt
(every pair shares 6-mers on both strands: reverse-complemented members, so that every library sequence has a block in both frames), 3 - 15 sequences, and more option combinations for the progressive stage after dndfast7.
Deterministic (python random.Random, fixed seed). Columns: id, options, seqs ('|' separated). Lower case and N for some sets."""
import random
rng = random.Random(75073)
N, P = 'ACGT', 'ACDEFGHIKLMNPQRSTVWY'
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def rc(s):
    return ''.join(COMP[c] for c in reversed(s))


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
    return ''.join(out) or alpha[0] * 2


OPTS_DNA = ['--fastapair', '--fastswpair', '--fastapair --retree 1', '--fastswpair --retree 1', '--fastapair --maxiterate 2', '--fastswpair --maxiterate 3', '--fastapair --treeout --distout', '--fastswpair --treeout',
            '--fastapair --reorder', '--fastapair --anysymbol', '--fastapair --op 1.0 --ep 0.1', '--fastapair --nuc', '--fastswpair --nuc --retree 2', '--fastapair --retree 2', '--fastapair --large', '--fastswpair --stepadd',
            '--fastapair --nofft', '--fastapair --leavegappyregion', '--fastapair --kappa 4', '--fastapair --weighti 0']
OPTS_PROT = ['--fastapair', '--fastswpair', '--fastapair --retree 2', '--fastapair --maxiterate 5', '--fastapair --bl 80', '--fastswpair --jtt 100', '--fastapair --large', '--fastapair --weighti 0', '--fastswpair --retree 3', '--fastapair --op 2.0']
rows = []
for i in range(240):
    prot = i % 6 == 5
    alpha = P if prot else N
    n = rng.choice([3, 4, 5, 6, 8, 10, 12, 15])
    length = rng.choice([40, 80, 150, 300, 600, 1000, 1500]) if not prot else rng.choice([40, 100, 300, 500])
    base = ''.join(rng.choice(alpha) for _ in range(length))
    seqs = []
    for k in range(n):
        s = mut(base, alpha, rng.choice([0.02, 0.08, 0.2, 0.4]), rng.choice([0.0, 0.02, 0.08]))
        if not prot and rng.random() < 0.35:
            s = rc(s)                                    # the reverse strand: blocks of the second frame
        if not prot and rng.random() < 0.25:
            a = rng.randrange(len(s) // 2 + 1)
            s = s[a:a + max(20, rng.randrange(len(s) // 2 + 1) + len(s) // 4)]    # a fragment
        seqs.append(s)
    if rng.random() < 0.2:
        seqs[rng.randrange(n)] = seqs[rng.randrange(n)]
    if rng.random() < 0.2:
        xa = 'X' if prot else 'N'
        seqs = [''.join(xa if rng.random() < 0.03 else c for c in s) for s in seqs]
    if rng.random() < 0.1:
        j = rng.randrange(n)
        seqs[j] = seqs[j].lower()
    seqs = [s if len(s) >= 2 else s + alpha[0] * 2 for s in seqs]
    opts = (OPTS_PROT if prot else OPTS_DNA)
    rows.append('%d\t%s\t%s' % (i + 1, opts[(i // 6 * 5 + i % 6) % len(opts)] if not prot else opts[(i // 6) % len(opts)], '|'.join(seqs)))
open('inputs/sets.tsv', 'w').write('\n'.join(rows) + '\n')
print('sets', len(rows))
