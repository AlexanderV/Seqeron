#!/usr/bin/env python3
"""Merges the result files of the E5 audit (round 1: results_e5_<file>.tsv, results_e5_p2_*.tsv, results_e5_q_*.tsv; round 2: results_e5_r2*_*.tsv, round 3: results_e5_r3_*.tsv) with the mutant lists (e5_*_mutants.json, e5_round2_*_mutants.json) into
results_e5.tsv: one row per mutant of the lists, verdict and note, and prints the counts. usage: merge_e5.py   (run in this directory)"""
import glob, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
EQUIVALENT = {
    'S06': 'equivalent: Advance is called only from Push, which sets the dirty flag of the top block right after it (`_curDirty = true`), so the value given to the new top block by Advance is never read',
}
REPLACED = {'S03': 'S03b', 'S04': 'S04b', 'S16': 'S16b', 'S17': 'S17b', 'S19': 'S19b', 'W07': 'W07b', 'W09': 'W09b', 'W13': 'W13b', 'W14': 'W14b', 'A07': 'A07b', 'A17': 'A17b', 'X04': 'X04b', 'E07': 'E07b', 'E12': 'E12b', 'D02': 'D02b', 'D04': 'D04b'}
descs, src = {}, {}
listed = []
for jf in sorted(glob.glob(os.path.join(HERE, 'e5_*_mutants.json'))):
    name = os.path.basename(jf)[3:-13]
    for m in json.load(open(jf)):
        descs.setdefault(m['id'], m['desc'])
        src.setdefault(m['id'], name)
        if not name.startswith('round2_') or m['id'] not in listed:
            if m['id'] not in listed:
                listed.append(m['id'])
r1, r2, r3 = {}, {}, {}
for rf in sorted(glob.glob(os.path.join(HERE, 'results_e5_*.tsv'))):
    base = os.path.basename(rf)
    if base == 'results_e5.tsv':
        continue
    target = r3 if base.startswith('results_e5_r3') else r2 if base.startswith('results_e5_r2') else r1
    for line in open(rf):
        c = line.rstrip('\n').split('\t')
        if len(c) >= 3:
            target[c[0]] = (c[2], c[3] if len(c) > 3 else '')
rows = ['id\tfile\tdescription\tround 1 (stage 1 / stage 2)\tround 2\tround 3\tverdict\tnote']
counts = {'listed': 0, 'compiled': 0, 'killed': 0, 'equivalent': 0, 'survivors': 0, 'not compiled': 0}
for i in listed:
    a = r1.get(i, ('', ''))
    b = r2.get(i, ('', ''))
    c3 = r3.get(i, ('', ''))
    counts['listed'] += 1
    killed = any(x.startswith('killed') for x in (a[0], a[1], b[0], b[1], c3[0], c3[1]))
    note = ''
    if killed:
        verdict = 'killed'
    elif i in EQUIVALENT:
        verdict, note = 'equivalent', EQUIVALENT[i]
    elif a[0].startswith('BUILD-ERROR') or a[0].startswith('PATTERN'):
        verdict = 'not compiled'
        note = ('replaced by ' + REPLACED[i]) if i in REPLACED else ''
    elif i.endswith('b') and b[0].startswith('BUILD-ERROR'):
        verdict = 'not compiled'
    elif a[0] == '' and b[0] == '':
        verdict = 'NOT RUN'
    else:
        verdict = 'SURVIVOR'
    if i in r2 and not (a[0].startswith('BUILD-ERROR') and i in REPLACED) and a[0] and not killed and verdict == 'SURVIVOR':
        note = 'survived again in round 2'
    if verdict == 'killed':
        counts['killed'] += 1
        counts['compiled'] += 1
    elif verdict == 'equivalent':
        counts['equivalent'] += 1
        counts['compiled'] += 1
    elif verdict == 'not compiled':
        counts['not compiled'] += 1
    else:
        counts['survivors'] += 1
        counts['compiled'] += 1
    rows.append('\t'.join([i, src.get(i, ''), descs.get(i, ''), (a[0] + ((' / ' + a[1]) if a[1] else '')), (b[0] + ((' / ' + b[1]) if b[1] else '')), (c3[0] + ((' / ' + c3[1]) if c3[1] else '')), verdict, note]))
open(os.path.join(HERE, 'results_e5.tsv'), 'w').write('\n'.join(rows) + '\n')
print(counts)
