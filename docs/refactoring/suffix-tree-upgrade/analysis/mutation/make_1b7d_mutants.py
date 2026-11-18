#!/usr/bin/env python3
"""Writes the hand-made single-site mutants of the new production code of work package 1b.7d (src/SuffixTree/Algorithms/SuffixTree.Esa/External: Dcs/*, Sort/*, Sa/SaBuilder, SeededRandom (with WorkPool), and the parts
   of MappedTextReader and Dist that the package added) as one JSON list per source file (1b7d_<file>_listed.json). Each pattern must occur exactly once in its file (checked here).
   usage: make_1b7d_mutants.py <repo root> [round prefix, default ""] """
import json, os, sys

ROOT = sys.argv[1]
ROUND = sys.argv[2] if len(sys.argv) > 2 else ''
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/External/'


def m(i, desc, old, new):
    return {'id': i, 'desc': desc, 'old': old, 'new': new}


FILES = {}
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mutants_1b7d_%sdata.py' % ROUND)).read())

bad = 0
for f, lst in FILES.items():
    path = os.path.join(ROOT, SRC, f)
    text = open(path).read()
    out = []
    for mm in lst:
        c = text.count(mm['old'])
        if c != 1:
            print('PATTERN', f, mm['id'], 'occurs', c)
            bad += 1
        mm = dict(mm, file=SRC + f)
        out.append(mm)
    name = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1b7d_%s%s_listed.json' % (ROUND, f.replace('/', '_').replace('.cs', '').lower()))
    json.dump(out, open(name, 'w'), indent=1)
    print(name, len(out))
sys.exit(1 if bad else 0)
