#!/usr/bin/env python3
"""A1-10d hand-off inputs: deterministic selection of committed ClustalW oracle sets (rows of the installed clustalw 2.1 are the control).

Writes inputs/sets.tsv:  id  type(DNA|PROTEIN)  seqs(;)  options(k=v,...)  controlRows(;)
Usage: python3 make_inputs.py            (run from this directory; reads the committed oracle TSVs of tests/.../TestData/Alignment)
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', '..', '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment')


def rows(name):
    with open(os.path.join(DATA, name)) as f:
        f.readline()
        for line in f:
            line = line.rstrip('\n')
            if line:
                p = line.split('\t')
                yield p[0], p[1], p[2], (p[3] if len(p) > 3 else '')


def main():
    os.makedirs(os.path.join(HERE, 'inputs'), exist_ok=True)
    out = []
    # protein: skip -NOVGAP (uninitialised member), -QUICKTREE variants are kept
    n = 0
    for i, seqs, r, o in rows('clustalw_protein_msa_oracle.tsv'):
        if 'novgap' in o.split(','):
            continue
        out.append((f'P{i}', 'PROTEIN', seqs, o, r))
        n += 1
        if n == 220:
            break
    n = 0
    for i, seqs, r, o in rows('clustalw_msa_full_oracle.tsv'):
        out.append((f'D{i}', 'DNA', seqs, o, r))
        n += 1
        if n == 120:
            break
    with open(os.path.join(HERE, 'inputs', 'sets.tsv'), 'w') as f:
        f.write('id\ttype\tseqs\toptions\tcontrol_rows\n')
        for t in out:
            f.write('\t'.join(t) + '\n')
    print(len(out), 'sets')


if __name__ == '__main__':
    main()
