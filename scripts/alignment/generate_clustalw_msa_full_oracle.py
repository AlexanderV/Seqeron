#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-10b) oracle: ClustalW 2.1 multiple-alignment stage for DNA.

Seqeron API: SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWOptions {...}).

Reference (run here): the installed `clustalw` 2.1 (Debian 2.1+lgpl-7) `-TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` with
the option variants -GAPOPEN -GAPEXT -TRANSWEIGHT -MAXDIV -NOWEIGHTS -NEGATIVE -QUICKTREE (-KTUPLE -WINDOW -TOPDIAGS -PAIRGAP)
-PWGAPOPEN -PWGAPEXT. The binary computes its own pairwise distances, neighbour-joining guide tree, sequence weights and profile
alignments; the C# port must reproduce the FASTA rows exactly (input order).

Output TSV (tests/.../TestData/Alignment/clustalw_msa_full_oracle.tsv):  id  seqs(;)  rows(;)  options(k=v,k=v)
Usage: python3 scripts/alignment/generate_clustalw_msa_full_oracle.py OUTDIR [nSets]      (needs `clustalw` on PATH)
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261123)
ALPHA = "ACGT"


def gen(kind, L, alpha):
    if kind == 0: return ''.join(rnd.choice(alpha) for _ in range(L))
    if kind == 1: return ''.join(rnd.choice(alpha[:2]) for _ in range(L))
    if kind == 2: return ''.join(rnd.choice([alpha[0]] * 3 + [alpha[1]]) for _ in range(L))
    u = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 4)))
    return (u * L)[:L]


def mutate(s, rate, alpha):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: continue
        elif r < rate:
            out.append(ch)
            for _ in range(rnd.choice([1, 1, 1, 2, 3])): out.append(rnd.choice(alpha))
        else: out.append(ch)
    return ''.join(out) or alpha[0]


def make_set(maxlen, alpha, kmin=2, kmax=9):
    k = rnd.randint(kmin, kmax)
    base = gen(rnd.choice([0, 0, 0, 1, 2, 3]), rnd.randint(8, maxlen), alpha)
    seqs = []
    div = rnd.choice([0.03, 0.1, 0.2, 0.35, 0.5])
    for _ in range(k):
        if rnd.random() < 0.8: seqs.append(mutate(base, rnd.choice([div, div * 0.5, div * 1.5]), alpha))
        else: seqs.append(gen(rnd.choice([0, 0, 1, 2, 3]), rnd.randint(5, maxlen), alpha))
    return seqs


def options():
    r = rnd.random()
    o = {}
    if r < 0.45: return o
    if rnd.random() < 0.3: o['gapopen'] = rnd.choice([5, 10, 15, 20, 25, 7.5])
    if rnd.random() < 0.3: o['gapext'] = rnd.choice([0.5, 1, 3, 6.66, 8, 10.25])
    if rnd.random() < 0.25: o['transweight'] = rnd.choice([0, 0.25, 0.5, 0.75, 1])
    if rnd.random() < 0.25: o['maxdiv'] = rnd.choice([0, 10, 30, 50, 70, 100])
    if rnd.random() < 0.15: o['noweights'] = None
    if rnd.random() < 0.2: o['negative'] = None
    if rnd.random() < 0.25:
        o['quick'] = None
        if rnd.random() < 0.5:
            k = rnd.choice([1, 2, 3, 4]); o['ktuple'] = k
            if rnd.random() < 0.5: o['window'] = rnd.randint(k + 1, k + 8)
            if rnd.random() < 0.5: o['topdiags'] = rnd.randint(k + 1, k + 8)
            if rnd.random() < 0.5: o['pairgap'] = rnd.randint(k + 1, k + 10)
    if rnd.random() < 0.15: o['pwgapopen'] = rnd.choice([5, 10, 15, 20])
    if rnd.random() < 0.15: o['pwgapext'] = rnd.choice([1, 3.33, 6.66, 9])
    return o


def run(seqs, o):
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'x.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs): f.write(f'>s{i}\n{s}\n')
        out = os.path.join(d, 'o.fa')
        args = ['clustalw', f'-INFILE={fa}', '-TYPE=DNA', '-ALIGN', '-OUTPUT=FASTA', '-OUTORDER=INPUT', f'-OUTFILE={out}']
        for k, v in o.items():
            flag = {'gapopen': 'GAPOPEN', 'gapext': 'GAPEXT', 'transweight': 'TRANSWEIGHT', 'maxdiv': 'MAXDIV', 'noweights': 'NOWEIGHTS',
                    'negative': 'NEGATIVE', 'quick': 'QUICKTREE', 'ktuple': 'KTUPLE', 'window': 'WINDOW', 'topdiags': 'TOPDIAGS',
                    'pairgap': 'PAIRGAP', 'pwgapopen': 'PWGAPOPEN', 'pwgapext': 'PWGAPEXT'}[k]
            args.append(f'-{flag}' + ('' if v is None else f'={v}'))
        r = subprocess.run(args, capture_output=True, text=True, cwd=d)
        rows = {}
        name = None
        for line in open(out):
            line = line.rstrip('\n')
            if line.startswith('>'): name = line[1:].strip(); rows[name] = ''
            elif name: rows[name] += line.strip()
    return [rows[f's{i}'] for i in range(len(seqs))]


def ambiguous(seqs):
    return [''.join(c if rnd.random() > 0.06 else rnd.choice("NRYKMSWBDHVU") for c in s) for s in seqs]


def main(outdir, nsets):
    os.makedirs(outdir, exist_ok=True)
    kinds = [(0.55, 60), (0.3, 120), (0.1, 250), (0.05, 600)]
    with open(os.path.join(outdir, 'clustalw_msa_full_oracle.tsv'), 'w') as f:
        f.write('id\tseqs\trows\toptions\n')
        for i in range(nsets):
            r = rnd.random(); acc = 0
            for p, ml in kinds:
                acc += p
                if r < acc: break
            seqs = make_set(ml, ALPHA, 2, rnd.choice([6, 9, 9, 14, 22]))
            if rnd.random() < 0.1: seqs = ambiguous(seqs)
            if rnd.random() < 0.05: seqs = seqs + [seqs[rnd.randrange(len(seqs))] for _ in range(rnd.randint(1, 3))]   # duplicates
            if rnd.random() < 0.04: seqs = [gen(rnd.randint(0, 3), rnd.randint(1, 6), ALPHA) for _ in range(rnd.randint(2, 6))]  # tiny
            o = options()
            rows = run(seqs, o)
            assert all(len(x) == len(rows[0]) for x in rows)
            ostr = ','.join(k if v is None else f'{k}={v}' for k, v in o.items())
            f.write(f"{i}\t{';'.join(seqs)}\t{';'.join(rows)}\t{ostr}\n")
    print(nsets, 'sets')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.', int(sys.argv[2]) if len(sys.argv) > 2 else 1000)
