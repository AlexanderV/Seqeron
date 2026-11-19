#!/usr/bin/env python3
"""Writes inputs/sets.tsv (id, type, seqs, options) for the ClustalW 2.1 `-CLUSTERING=UPGMA -ALIGN` oracle run (deterministic: seed 20261119; the committed file is reproduced exactly).
600 sets, alternating DNA / protein, 2..9 sequences of 8..90 residues from a common ancestor with substitutions, indels and (10 %) unrelated sequences, random options
(-NOWEIGHTS -MAXDIV -GAPOPEN -GAPEXT -MATRIX -PWMATRIX -QUICKTREE). Usage: python3 make_inputs.py
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261119)
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


def make_set(t):
    k = rnd.choice([2, 3, 3, 4, 4, 5, 6, 7, 9])
    base = ''.join(res(t) for _ in range(rnd.randint(8, 90)))
    div = rnd.choice([0.03, 0.1, 0.2, 0.35, 0.55, 0.8])
    return [mutate(t, base, rnd.choice([div, div / 2, div * 1.3])) if rnd.random() < 0.9 else ''.join(res(t) for _ in range(rnd.randint(8, 60))) for _ in range(k)]


def options(t):
    o = {}
    if rnd.random() < 0.15: o['noweights'] = None
    if rnd.random() < 0.25: o['maxdiv'] = rnd.choice([0, 10, 30, 60, 100])
    if rnd.random() < 0.2: o['gapopen'] = rnd.choice([5, 10, 15, 25] if t == 'protein' else [8, 15, 25])
    if rnd.random() < 0.2: o['gapext'] = rnd.choice([0.1, 0.5, 1, 3.3, 6.66])
    if t == 'protein' and rnd.random() < 0.25: o['matrix'] = rnd.choice(['BLOSUM', 'PAM', 'ID', 'GONNET'])
    if t == 'protein' and rnd.random() < 0.1: o['pwmatrix'] = rnd.choice(['BLOSUM', 'PAM', 'ID', 'GONNET'])
    if rnd.random() < 0.1: o['quicktree'] = None
    return o


def main():
    n = 600
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(os.path.join(here, 'inputs'), exist_ok=True)
    with open(os.path.join(here, 'inputs', 'sets.tsv'), 'w') as f:
        f.write('id\ttype\tseqs\toptions\n')
        for i in range(n):
            t = 'dna' if i % 2 == 0 else 'protein'
            seqs = make_set(t)
            o = options(t)
            f.write(f'{i}\t{t}\t{";".join(seqs)}\t{",".join(k if v is None else f"{k}={v}" for k, v in o.items())}\n')
    print(n, 'sets written to inputs/sets.tsv')


if __name__ == '__main__':
    main()
