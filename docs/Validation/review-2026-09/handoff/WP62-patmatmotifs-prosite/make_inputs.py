#!/usr/bin/env python3
"""Writes the inputs of the run (inputs/*.fa) deterministically: protein sequences (seed 20261122) with planted instances of frequent PROSITE patterns (zinc finger C2H2, P-loop, N-glycosylation, protein kinase C / casein kinase II /
cAMP phosphorylation, N-myristoylation, leucine zipper, EF-hand calcium-binding, ATP/GTP-binding site A, tyrosine phosphorylation) and random residues, one file with many sequences, one file per option set, damaged input.
Usage: python3 make_inputs.py DIR"""
import os, random, sys

AA = 'ACDEFGHIKLMNPQRSTVWY'
PLANTS = [
    'CAACAAAFLAAAAAAAALAAAAAAAAHAAAH',                 # C2H2 zinc finger (PS00028)
    'GAAAAGKS', 'AGAAAAGKT',                           # P-loop / ATP-GTP binding site A (PS00017)
    'NASA', 'NTSA', 'NVTL',                            # N-glycosylation (PS00001)
    'SAR', 'TAK', 'SAAD', 'TAAE', 'RRAS', 'KRAT',      # PKC / CK2 / cAMP phospho sites (PS00005 PS00006 PS00004)
    'GGAAAS', 'GNAAAT',                                # N-myristoylation (PS00008)
    'LAAAAAALAAAAAAALAAAAAAALAAAAAAAL',                # leucine zipper-like (PS00029)
    'DKDGDGTITTKELGTVMRSLGQNPTEAELQDMING',            # EF-hand (PS00018)
    'KAAAAAAEY', 'RAAAAAEEEY',                         # tyrosine kinase phosphorylation (PS00007)
]


def seq(rnd, n):
    return ''.join(rnd.choice(AA) for _ in range(n))


def main(d):
    rnd = random.Random(20261122)
    os.makedirs(d, exist_ok=True)
    recs = []
    for i in range(60):
        parts = [seq(rnd, rnd.randint(5, 60))]
        for _ in range(rnd.choice([0, 1, 2, 3])):
            parts += [rnd.choice(PLANTS), seq(rnd, rnd.randint(3, 50))]
        s = ''.join(parts)
        if rnd.random() < 0.15: s = s.lower()
        recs.append(('seq%03d' % i, 'planted test sequence %d' % i, s))
    recs += [('short1', 'one residue', 'M'), ('short2', 'two', 'MK'), ('ambig', 'with B Z X', 'MKBZXAAAGKSAAA' + seq(rnd, 30)), ('star', 'stop', 'MKV*AAAGKS')]
    def fa(rs, w=60):
        return ''.join('>%s %s\n%s\n' % (n, desc, '\n'.join(s[j:j + w] for j in range(0, len(s), w)) or '') for n, desc, s in rs)
    open(os.path.join(d, 'many.fa'), 'w').write(fa(recs))
    for i, name in enumerate(['one.fa', 'two.fa', 'five.fa']):
        open(os.path.join(d, name), 'w').write(fa(recs[i * 7: i * 7 + (1, 2, 5)[i]]))
    open(os.path.join(d, 'damaged.fa'), 'w').write('>x\nMKV\x00AAA\n>y\n\n>z desc only\n')
    open(os.path.join(d, 'empty.fa'), 'w').write('')


if __name__ == '__main__':
    main(sys.argv[1])
