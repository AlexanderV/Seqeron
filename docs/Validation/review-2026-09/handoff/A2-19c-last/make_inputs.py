#!/usr/bin/env python3
"""Writes inputs/sets.tsv for the run of `mafft --lastpair` / `--lastmultipair` (B09 WP77, package A2-19c-last): 206 sets, deterministic (python random.Random, fixed seed).
Columns (tab separated): id, mode, options, old sequences ('|' separated; the file in.fa; for the modes with `add` an alignment of them is made by the installed mafft inside the container), new sequences ('|' separated, the file new.fa; empty without --add*), names flag.
Modes: plain (in.fa = the unaligned sequences), add (in.fa = `mafft --quiet` alignment of the old sequences, new.fa = the sequences to add; the options hold the --addfragments / --addfull / --addlong / --addtotop option and `new.fa`).
The mafft script reads the nucleotide sequences only for --lastpair / --lastmultipair (pairlocalalign stops with "Not yet supported" for amino acids): some protein sets record that message, too."""
import random, sys
rng = random.Random(77019)
N, P = 'ACGT', 'ACDEFGHIKLMNPQRSTVWY'
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}


def rc(s):
    return ''.join(COMP.get(c, c) for c in reversed(s))


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


def family(alpha, n, length, r=None, ind=None):
    base = ''.join(rng.choice(alpha) for _ in range(length))
    r = rng.choice([0.03, 0.08, 0.18]) if r is None else r
    ind = rng.choice([0.0, 0.02, 0.06]) if ind is None else ind
    return base, [mut(base, alpha, r, ind) for _ in range(n)], r, ind


PAIR = ['--lastpair', '--lastmultipair']
PLAIN_EXTRAS = [[], [], [], ['--retree', '1'], ['--retree', '2'], ['--maxiterate', '2'], ['--maxiterate', '1000'], ['--treeout', '--distout'], ['--reorder'], ['--op', '1.5', '--ep', '0.2'], ['--nofft'], ['--anysymbol'], ['--thread', '2'],
                ['--clustalout'], ['--leavegappyregion'], ['--kappa', '4'], ['--weighti', '0'], ['--bl', '80'], ['--nuc'], ['--quiet']]
ADD_KINDS = ['--addfragments', '--addfull', '--addlong', '--addtotop']
ADD_EXTRAS = [[], [], ['--keeplength'], ['--mapout'], ['--reorder'], ['--treeout', '--distout'], ['--op', '1.5', '--ep', '0.2'], ['--retree', '1'], ['--nofft'], ['--thread', '2'], ['--anysymbol'], ['--compactmapout', '--keeplength']]
rows = []


def add(mode, opts, old, new=()):
    rows.append((('set%03d' % len(rows)), mode, ' '.join(opts), '|'.join(old), '|'.join(new)))


# 1. plain DNA families (70)
for i in range(70):
    n = rng.choice([2, 3, 4, 5, 6, 8, 10])
    length = rng.choice([60, 120, 250, 400, 600])
    base, fam, r, ind = family(N, n, length)
    if rng.random() < 0.3:
        fam = [rc(s) if rng.random() < 0.4 else s for s in fam]      # reverse complemented members (lastal -s 1: forward strand only)
    if rng.random() < 0.1:
        fam = [s.lower() for s in fam]
    if rng.random() < 0.1:
        fam[0] = fam[0][:5] + 'N' * rng.randint(1, 8) + fam[0][5:]
    add('plain', ['--lastmultipair' if i % 9 == 8 else '--lastpair'] + rng.choice(PLAIN_EXTRAS), fam)       # --lastmultipair needs --addfragments: the script stops with "Specify --addfragments too"
# 2. add modes (60): an alignment plus new sequences, --lastpair / --lastmultipair with the four fragment options
for i in range(60):
    base, fam, r, ind = family(N, rng.randint(2, 6), rng.choice([80, 160, 300]))
    kind = ADD_KINDS[i % 4]
    new = []
    for _ in range(rng.randint(1, 5)):
        if kind in ('--addfragments', '--addtotop') and rng.random() < 0.8:
            a = rng.randrange(max(1, len(base) - 30))
            piece = base[a:a + rng.randint(25, max(26, len(base) // 2))]
        else:
            piece = base
        new.append(mut(piece, N, r, ind))
    add('add', [kind, 'new.fa', PAIR[(i // 4) % 2]] + rng.choice(ADD_EXTRAS), fam, new)
# 3. protein sets (20): the programs stop with "Not yet supported"
for i in range(20):
    base, fam, r, ind = family(P, rng.randint(2, 6), rng.choice([40, 100, 250]))
    add('plain', ['--lastpair'] + rng.choice(PLAIN_EXTRAS[:6]), fam)
for i in range(6):
    base, fam, r, ind = family(P, rng.randint(2, 4), 80)
    add('add', [ADD_KINDS[i % 4], 'new.fa', PAIR[i % 2]], fam, [mut(base, P, 0.1, 0.02)])
# 4. odd sets (30): one sequence, two sequences, identical, short, N rich, long, many
for i in range(8):
    add('plain', ['--lastpair'], [''.join(rng.choice(N) for _ in range(rng.randint(1, 25))) for _ in range(rng.randint(2, 6))])
for i in range(6):
    s = ''.join(rng.choice(N) for _ in range(rng.choice([50, 150, 400])))
    add('plain', ['--lastpair'], [s] * rng.randint(2, 5))
for i in range(6):
    base, fam, r, ind = family(N, 4, 200)
    fam = [''.join('N' if rng.random() < 0.3 else c for c in s) for s in fam]
    add('plain', ['--lastpair'], fam)
for i in range(4):
    base, fam, r, ind = family(N, rng.choice([15, 25]), rng.choice([300, 800]), 0.05, 0.02)
    add('plain', ['--lastpair'] + rng.choice([[], ['--retree', '1']]), fam)
for i in range(4):
    base, fam, r, ind = family(N, 3, 1500, 0.08, 0.03)
    add('plain', ['--lastpair'], fam)
# 5. two related sequences (20): `pairlocalalign -R` puts the pair (0, 1) into row 0 of its table, which has room for it, so the program does not overflow its heap (valgrind, README) and the run reaches the end: the success path
for i in range(20):
    base, fam, r, ind = family(N, 2, rng.choice([60, 150, 300, 600, 1200]))
    add('plain', ['--lastpair'] + rng.choice(PLAIN_EXTRAS), fam)
with open('inputs/sets.tsv', 'w') as f:
    for r in rows:
        f.write('\t'.join(r) + '\n')
print('sets', len(rows), file=sys.stderr)
