#!/usr/bin/env python3
"""Mutation audit of work package 1b.9h: the result files results_1b9h_<file>[_<stage>].tsv of run_1b9h_mutants.py against the lists 1b9h_<file>_listed.json: one line per mutant (its verdict, the closure by a test or an argued equivalence
   from 1b9h_closed.py) and the counts. usage: summarize_1b9h.py   (reads every results_1b9h_*.tsv of this directory except the G1 files)"""
import glob, json, os, re, collections
HERE = os.path.dirname(os.path.abspath(__file__))
CLOSED = {}
p = os.path.join(HERE, '1b9h_closed.py')
if os.path.exists(p):
    exec(open(p).read())  # defines CLOSED = {id: ('test' | 'equivalent', text)}
listed = {}
for f in sorted(glob.glob(os.path.join(HERE, '1b9h_*_listed.json'))):
    if 'g1_' in f:
        continue
    for m in json.load(open(f)):
        listed[m['id']] = (os.path.basename(f), m['desc'])
verdict = {}
for f in sorted(glob.glob(os.path.join(HERE, 'results_1b9h_*.tsv'))):
    if 'g1_' in f:
        continue
    for line in open(f):
        q = line.rstrip('\n').split('\t')
        if len(q) >= 3 and q[0] in listed:
            old = verdict.get(q[0], '')
            new = q[2] + (' | stage 2: ' + q[3] if len(q) > 3 and q[3] else '')
            if not old.startswith('killed') or new.startswith('killed'):
                verdict[q[0]] = new
cnt = collections.Counter()
print('id\tfile\tdesc\tverdict\tclosure')
for k, (f, d) in sorted(listed.items()):
    v = verdict.get(k, 'not run')
    status = 'killed' if v.startswith('killed') or 'killed(' in v else ('closed' if k in CLOSED else 'OPEN')
    cnt[(f, status)] += 1
    print('\t'.join([k, f, d, v, CLOSED.get(k, ('', ''))[0] + ': ' + CLOSED.get(k, ('', ''))[1][:120] if k in CLOSED else '']))
print()
for (f, s), n in sorted(cnt.items()):
    print(f, s, n)
