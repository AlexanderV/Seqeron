#!/usr/bin/env python3
"""A1-10c-muscle ROUND 3 hand-off inputs (a small follow-up of round 2, 630 of 638 sets of which are equal to the port).

Writes inputs/sets.tsv with the columns  id  kind  type  data  options  (same format as round 2, see ../A1-10c-muscle-r2/make_inputs.py):
  1. the 8 sets of round 2 whose console output the round-2 harness could not decode (-matrix blosum62nostar: the program prints 0xff bytes in its "Matrix is not symmetrical" warning); same ids and
     data (read from ../A1-10c-muscle-r2/inputs/sets.tsv) - the round-3 harness records such bytes as \\xNN;
  2. -matrix with nucleotide sequences (the program reads the matrix with the nucleotide alphabet, then switches g_Alpha to amino acid and the score to SP: the behaviour is only known from the source);
  3. one- and two-sequence alignments and profiles in every command line (-root1 / -root2 non-pseudo on one-sequence profiles crashed in round 2: all other combinations are recorded);
  4. align with two sequences (and one) with the tree files, -usetree, roots, clusterings, -computeweights;
  5. the score flags -le / -sp / -sv / -spn on nucleotides and proteins.
Usage: python3 make_inputs.py   (deterministic; needs ../A1-10c-muscle-r2/make_inputs.py and its inputs/sets.tsv)
"""
import importlib.util, os, random

HERE = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.join(HERE, '..', 'A1-10c-muscle-r2')
spec = importlib.util.spec_from_file_location('m2', os.path.join(R2, 'make_inputs.py'))
m2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m2)
m2.rnd = random.Random(20261119)
rows_out = []


def add(kind, typ, data, d):
    rows_out.append((f'R{638 + len(rows_out) - 8}', kind, typ, data, m2.opt(d)))


def main():
    r2 = {}
    for line in list(open(os.path.join(R2, 'inputs', 'sets.tsv')))[1:]:
        sid, kind, typ, data, opts = line.rstrip('\n').split('\t')
        r2[sid] = (sid, kind, typ, data, opts)
    for sid in ['R330', 'R335', 'R336', 'R342', 'R344', 'R347', 'R348', 'R350']:
        assert 'blosum62nostar' in r2[sid][4], sid
        rows_out.append(r2[sid])
    U, F = m2.unaligned, m2.family
    # 2. -matrix with nucleotides
    for i, typ in enumerate(['dna', 'dna', 'dna', 'dna', 'rna', 'rna', 'dna', 'rna']):
        add('align', typ, ';'.join(U(typ, 3, 6, 60)), {'maxiters': [1, 2, 4, 8][i % 4], 'matrix': ['@blosum62', '@pam250', '@blosum62nostar'][i % 3]})
    # 3. one- and two-sequence profiles, alignments, scores, trees
    def msa(typ, k, n=40): return m2.msa_family(typ, n, k, 0.2)
    for typ, k1, k2, root in [('protein', 1, 3, 'minavgleafdist'), ('dna', 3, 1, 'minavgleafdist'), ('protein', 1, 1, 'midlongestspan'), ('dna', 1, 2, 'midlongestspan'), ('protein', 2, 2, 'minavgleafdist')]:
        both = msa(typ, k1 + k2)
        add('profile', typ, ';'.join(both[:k1]) + '|' + ';'.join(both[k1:]), {'root1': root})
    for typ, k, d in [('protein', 1, {}), ('protein', 1, {'root2': 'midlongestspan'}), ('dna', 2, {}), ('protein', 2, {'root2': 'minavgleafdist'}), ('dna', 1, {'cluster2': 'neighborjoining'})]:
        add('maketree', typ, ';'.join(msa(typ, k)), d)
    for typ, k, d in [('protein', 1, {}), ('dna', 2, {}), ('protein', 2, {'weight1': 'gsc'}), ('dna', 1, {'root2': 'midlongestspan'})]:
        add('spscore', typ, ';'.join(msa(typ, k)), d)
    for typ, k, d in [('protein', 1, {'maxiters': 3}), ('protein', 2, {'maxiters': 3}), ('dna', 2, {'maxiters': 8}), ('protein', 2, {'maxiters': 4, 'noanchors': None}), ('dna', 1, {})]:
        add('refine', typ, ';'.join(msa(typ, k)), d)
    for typ, k1, k2, d in [('dna', 1, 1, {'root1': 'midlongestspan'}), ('dna', 2, 1, {}), ('rna', 1, 2, {'root1': 'minavgleafdist'}), ('dna', 3, 1, {'cluster2': 'neighborjoining'})]:
        both = msa(typ, k1 + k2)
        add('ppscore', typ, ';'.join(both[:k1]) + '|' + ';'.join(both[k1:]), d)
    # 4. align with one and two sequences
    for typ, k, d in [('protein', 2, {}), ('protein', 2, {'maxiters': 8}), ('dna', 2, {'maxiters': 3, 'cluster1': 'neighborjoining'}), ('protein', 2, {'root1': 'midlongestspan'}),
                      ('protein', 2, {'distance1': 'kmer20_3'}), ('dna', 2, {'tree1': 1, 'tree2': 1, 'distmx1': 1, 'distmx2': 1, 'maxiters': 2}), ('protein', 2, {'usetree': '@own'}),
                      ('protein', 2, {'weight1': 'gsc', 'maxiters': 3}), ('dna', 2, {'computeweights': 1}), ('protein', 1, {'tree1': 1}), ('dna', 1, {'maxiters': 4, 'matrix': None}),
                      ('protein', 3, {'root1': 'minavgleafdist', 'maxiters': 1}), ('protein', 3, {'usetree': '@own', 'maxiters': 8})]:
        d = {k_: v for k_, v in d.items() if not (k_ == 'matrix')}
        add('align', typ, ';'.join(U(typ, k, k, 50)), d)
    # 5. score flags
    for typ, flag in [('dna', 'le'), ('dna', 'sp'), ('dna', 'sv'), ('rna', 'sp'), ('protein', 'spn'), ('protein', 'sv'), ('protein', 'sp'), ('protein', 'le')]:
        add('align', typ, ';'.join(U(typ, 3, 6, 60)), {'maxiters': 2, flag: None})
    os.makedirs(os.path.join(HERE, 'inputs'), exist_ok=True)
    with open(os.path.join(HERE, 'inputs', 'sets.tsv'), 'w') as f:
        f.write('id\tkind\ttype\tdata\toptions\n')
        for r in rows_out:
            f.write('\t'.join(r) + '\n')
    from collections import Counter
    print(len(rows_out), 'sets:', dict(Counter(r[1] for r in rows_out)))


if __name__ == '__main__':
    main()
