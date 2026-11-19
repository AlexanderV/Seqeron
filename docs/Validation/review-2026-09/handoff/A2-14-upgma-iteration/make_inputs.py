#!/usr/bin/env python3
"""Writes inputs/sets.tsv (id, mode, type, seqs, options) for the ClustalW 2.1 `-CLUSTERING=UPGMA` run of the paths A2-14-clustalw-upgma did not reach (deterministic: seed 20261120; the committed file is
reproduced exactly). 800 sets, alternating DNA / protein; modes: `align-alignment` (-ITERATION=alignment -NUMITER=n, 250 sets), `align-tree` (-ITERATION=tree, 200 sets), `profile` (-PROFILE1 -PROFILE2 -PROFILE, 175 sets),
`sequences` (-PROFILE1 -PROFILE2 -SEQUENCES, 175 sets). `seqs` = the unaligned sequences (';'); for `profile` / `sequences` the first `n1` (the `n1=` option) form profile 1, the rest profile 2 (aligned by the program
itself, NJ, before the command: run_sets.py). Options: -NOWEIGHTS -MAXDIV -GAPOPEN -GAPEXT -MATRIX -NUMITER -QUICKTREE. Usage: python3 make_inputs.py
"""
import os, random

rnd = random.Random(20261120)
FREQ = dict(A=8.25, R=5.53, N=4.06, D=5.45, C=1.37, Q=3.93, E=6.75, G=7.07, H=2.27, I=5.96, L=9.66, K=5.84, M=2.42, F=3.86, P=4.70, S=6.56, T=5.34, W=1.08, Y=2.92, V=6.87)
AA = list(FREQ); W = list(FREQ.values())


def res(t): return rnd.choices(AA, W)[0] if t == 'protein' else rnd.choice('ACGT')


def mutate(t, s, rate):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.55: out.append(res(t))
        elif r < rate * 0.75: continue
        elif r < rate:
            out.append(ch); out += [res(t) for _ in range(rnd.choice([1, 1, 2, 4]))]
        else: out.append(ch)
    return ''.join(out) or res(t)


def make_set(t, kmin, kmax):
    k = rnd.randint(kmin, kmax)
    base = ''.join(res(t) for _ in range(rnd.randint(20, 100)))
    div = rnd.choice([0.03, 0.1, 0.2, 0.35, 0.55])
    return [mutate(t, base, rnd.choice([div, div / 2, div * 1.3])) if rnd.random() < 0.9 else ''.join(res(t) for _ in range(rnd.randint(15, 60))) for _ in range(k)]


def options(t, mode):
    o = {}
    if mode.startswith('align') and rnd.random() < 0.4: o['numiter'] = rnd.choice([1, 2, 3, 5])
    if rnd.random() < 0.15: o['noweights'] = None
    if rnd.random() < 0.15: o['maxdiv'] = rnd.choice([0, 10, 30, 60, 100])
    if rnd.random() < 0.2: o['gapopen'] = rnd.choice([5, 10, 15, 25] if t == 'protein' else [8, 15, 25])
    if rnd.random() < 0.2: o['gapext'] = rnd.choice([0.1, 0.5, 1, 3.3, 6.66])
    if t == 'protein' and rnd.random() < 0.25: o['matrix'] = rnd.choice(['BLOSUM', 'PAM', 'ID', 'GONNET'])
    if rnd.random() < 0.1: o['quicktree'] = None
    return o


def main():
    modes = ['align-alignment'] * 250 + ['align-tree'] * 200 + ['profile'] * 175 + ['sequences'] * 175
    rnd.shuffle(modes)
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(os.path.join(here, 'inputs'), exist_ok=True)
    with open(os.path.join(here, 'inputs', 'sets.tsv'), 'w') as f:
        f.write('id\tmode\ttype\tseqs\toptions\n')
        for i, mode in enumerate(modes):
            t = 'dna' if i % 2 == 0 else 'protein'
            seqs = make_set(t, 4, 10) if mode.startswith('align') else make_set(t, 3, 9)
            o = options(t, mode)
            if not mode.startswith('align'):
                o['n1'] = rnd.randint(1, len(seqs) - 1)
            f.write(f'{i}\t{mode}\t{t}\t{";".join(seqs)}\t{",".join(k if v is None else f"{k}={v}" for k, v in o.items())}\n')
    print(len(modes), 'sets written to inputs/sets.tsv')


if __name__ == '__main__':
    main()
