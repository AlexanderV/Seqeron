#!/usr/bin/env python3
"""A1-10c-muscle hand-off inputs: 420 deterministic sets for MUSCLE 3.8.1551 (protein 270, DNA 110, RNA 40) with the option variants the port supports.

Writes inputs/sets.tsv:  id  type(protein|dna|rna)  seqs(;)  options(k=v,...)
options: maxiters (1|2), model (le|sp|sv|spn), gapopen, gapextend, center, distance1 (kmer6_6|kmer4_6), cluster1, cluster2 (upgmb|upgma|upgmamin|upgmamax), sueff, maxtrees (0).
The sets cover: families of 2-85 % divergence, 2 sequences (stops after the progressive stage), 3+ (tree refinement), identical sequences (tied distances), sequences shorter than
5 residues (k-mer distance fallback), sequences starting with M (the methionine hack), wildcards B / X / Z / J / O (protein) and N / R / Y / X / K (nucleotide: the AX_R / AX_Y quirk),
low-complexity runs, lower-case input letters, U in DNA.
Usage: python3 make_inputs.py            (deterministic: the committed inputs/sets.tsv is reproduced exactly)
"""
import os, random

rnd = random.Random(31337)
HERE = os.path.dirname(os.path.abspath(__file__))
FREQ = dict(A=8.25, R=5.53, N=4.06, D=5.45, C=1.37, Q=3.93, E=6.75, G=7.07, H=2.27, I=5.96, L=9.66, K=5.84, M=2.42, F=3.86, P=4.70, S=6.56, T=5.34, W=1.08, Y=2.92, V=6.87)
AA = list(FREQ)
W = [FREQ[a] for a in AA]


def rprot(): return rnd.choices(AA, W)[0]


def mutate_p(s, rate):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.55: out.append(rprot())
        elif r < rate * 0.75: continue
        elif r < rate:
            out.append(ch)
            for _ in range(rnd.choice([1, 1, 1, 2, 3, 5])): out.append(rprot())
        else: out.append(ch)
    return ''.join(out) or 'A'


def mutate_n(s, rate, alpha='ACGT'):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.55: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: continue
        elif r < rate:
            out.append(ch)
            for _ in range(rnd.choice([1, 1, 2, 3, 6])): out.append(rnd.choice(alpha))
        else: out.append(ch)
    return ''.join(out) or 'A'


def family(typ, maxlen, kmin, kmax):
    k = rnd.randint(kmin, kmax)
    div = rnd.choice([0.02, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85])
    alpha = 'ACGU' if typ == 'rna' else 'ACGT'
    if typ == 'protein':
        kind = rnd.choice([0, 0, 0, 0, 1, 2])
        L = rnd.randint(12, maxlen)
        if kind == 0: base = ''.join(rprot() for _ in range(L))
        elif kind == 1: base = ''.join(rnd.choice('AAAL') for _ in range(L))
        else: base = (''.join(rprot() for _ in range(rnd.randint(2, 6))) * L)[:L]
        mut = mutate_p
    else:
        L = rnd.randint(20, maxlen * 2)
        base = ''.join(rnd.choice(alpha) for _ in range(L))
        mut = lambda s, r: mutate_n(s, r, alpha)
    seqs = []
    for _ in range(k):
        if rnd.random() < 0.88: seqs.append(mut(base, rnd.choice([div, div * 0.5, div * 1.3])))
        else: seqs.append(''.join(rprot() for _ in range(rnd.randint(5, maxlen))) if typ == 'protein' else ''.join(rnd.choice(alpha) for _ in range(rnd.randint(8, maxlen))))
    return seqs


def decorate(typ, seqs):
    r = rnd.random()
    if typ == 'protein':
        if r < 0.25:   # leading M (methionine hack)
            seqs = [('M' + s if rnd.random() < 0.7 else s) for s in seqs]
        elif r < 0.35:  # wildcards and invalid letters
            seqs = [''.join(c if rnd.random() > 0.04 else rnd.choice('BXZJOU') for c in s) for s in seqs]
        elif r < 0.40:  # lower case
            seqs = [s.lower() if rnd.random() < 0.5 else s for s in seqs]
    else:
        if r < 0.15:
            seqs = [''.join(c if rnd.random() > 0.04 else rnd.choice('NRYXKMSW') for c in s) for s in seqs]
        elif r < 0.22:
            seqs = [s.lower() if rnd.random() < 0.5 else s for s in seqs]
    if rnd.random() < 0.06 and len(seqs) >= 2:   # duplicates (tied distances)
        seqs = seqs + [seqs[rnd.randrange(len(seqs))] for _ in range(rnd.randint(1, 3))]
    if rnd.random() < 0.05:                      # very short sequences (< 5 residues)
        seqs = seqs + [''.join(rnd.choice(AA if typ == 'protein' else 'ACGT') for _ in range(rnd.randint(1, 4))) for _ in range(rnd.randint(1, 2))]
    return seqs


def options(typ):
    o = {}
    o['maxiters'] = rnd.choice([1, 2, 2, 2])
    if rnd.random() < 0.15: o['maxtrees'] = 0
    if typ == 'protein':
        if rnd.random() < 0.2: o['model'] = rnd.choice(['sp', 'sv', 'le'])
        if rnd.random() < 0.2: o['gapopen'] = rnd.choice(['-1.5', '-2.9', '-4', '-6'])
        if rnd.random() < 0.1: o['gapextend'] = rnd.choice(['-0.1', '-0.5', '-1'])
        if rnd.random() < 0.1: o['center'] = rnd.choice(['-0.3', '-0.52', '-1'])
    else:
        if rnd.random() < 0.2: o['gapopen'] = rnd.choice(['-200', '-400', '-700'])
        if rnd.random() < 0.1: o['gapextend'] = rnd.choice(['-10', '-30', '-60'])
    if rnd.random() < 0.15: o['cluster1'] = rnd.choice(['upgma', 'upgmamin', 'upgmamax', 'upgmb'])
    if rnd.random() < 0.15: o['cluster2'] = rnd.choice(['upgma', 'upgmamin', 'upgmamax', 'upgmb'])
    if rnd.random() < 0.1: o['sueff'] = rnd.choice(['0.0', '0.3', '0.5', '1.0'])
    if rnd.random() < 0.08: o['distance1'] = 'kmer6_6'
    return o


def main():
    rows = []
    for i in range(420):
        typ = 'protein' if i < 270 else 'dna' if i < 380 else 'rna'
        maxlen = rnd.choice([40, 70, 70, 130, 130, 220])
        seqs = family(typ, maxlen, 2, rnd.choice([3, 4, 5, 6, 8, 8, 12, 20]))
        if rnd.random() < 0.1: seqs = seqs[:2]
        seqs = decorate(typ, seqs)
        rows.append((f'M{i}', typ, ';'.join(seqs), ','.join(f'{k}={v}' for k, v in options(typ).items())))
    os.makedirs(os.path.join(HERE, 'inputs'), exist_ok=True)
    with open(os.path.join(HERE, 'inputs', 'sets.tsv'), 'w') as f:
        f.write('id\ttype\tseqs\toptions\n')
        for r in rows:
            f.write('\t'.join(r) + '\n')
    print(len(rows), 'sets')


if __name__ == '__main__':
    main()
