#!/usr/bin/env python3
"""Merges the result files of the E6 audit (round 1: results_e6_<A|B><file>.tsv with the lists e6_<file>_mutants.json; round 2, the replacements of the mutants that did not compile: results_e6_R2*_*.tsv with e6_round2_*_mutants.json)
into results_e6.tsv: one row per mutant that was run, verdict and note, and prints the counts. usage: merge_e6.py   (run in this directory)"""
import glob, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
REPLACED = {'P13': 'P13b', 'P15': 'P15b', 'R01': 'R01b', 'R05': 'R05b'}
EQUIVALENT = {}
descs, listed = {}, []
for jf in sorted(glob.glob(os.path.join(HERE, 'e6_*_mutants.json'))):
    for m in json.load(open(jf)):
        descs.setdefault(m['id'], m['desc'])
        if m['id'] not in listed:
            listed.append(m['id'])
r1, r2 = {}, {}
for rf in sorted(glob.glob(os.path.join(HERE, 'results_e6_*.tsv'))):
    base = os.path.basename(rf)
    if base == 'results_e6.tsv':
        continue
    target = r2 if base.startswith('results_e6_R2') else r1
    for line in open(rf):
        c = line.rstrip('\n').split('\t')
        if len(c) >= 3:
            target[c[0]] = (c[2], c[3] if len(c) > 3 else '')
rows = ['id\tdescription\tround 1 (stage 1 / stage 2)\tround 2\tverdict\tnote']
counts = {'listed': 0, 'compiled': 0, 'killed': 0, 'equivalent': 0, 'survivors': 0, 'not compiled': 0}
for i in listed:
    a, b = r1.get(i, ('', '')), r2.get(i, ('', ''))
    counts['listed'] += 1
    notcompiled = a[0].startswith('BUILD-') or b[0].startswith('BUILD-')
    killed = any(x.startswith('killed') for x in (a[0], a[1], b[0], b[1]))
    note = ''
    if a[0].startswith('BUILD-ERROR') and i in REPLACED:
        verdict, note = 'not compiled', 'replaced by ' + REPLACED[i] + ' (analyzer error ' + a[0].split()[-1] + ')'
        counts['not compiled'] += 1
    elif killed:
        verdict = 'killed'
        counts['compiled'] += 1
        counts['killed'] += 1
    elif i in EQUIVALENT:
        verdict, note = 'equivalent', EQUIVALENT[i]
        counts['compiled'] += 1
        counts['equivalent'] += 1
    elif a[0] == '' and b[0] == '':
        verdict = 'not run'
    else:
        verdict = 'SURVIVED'
        counts['compiled'] += 1
        counts['survivors'] += 1
    rows.append('\t'.join([i, descs[i], ' / '.join(x for x in a if x), ' / '.join(x for x in b if x), verdict, note]))
open(os.path.join(HERE, 'results_e6.tsv'), 'w').write('\n'.join(rows) + '\n')
print(counts)
