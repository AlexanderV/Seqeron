#!/usr/bin/env python3
"""Merges the result files of the E4 audit (round 1: results_e4_<file>.tsv, round 2: results_e4_round2_<file>.tsv) with the mutant lists into results_e4.tsv: one row per mutant of the lists, verdict and note.
usage: merge_e4.py   (run in this directory)"""
import glob, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
# arguments for the survivors that are equivalent, and where a mutant that did not compile was replaced
EQUIVALENT = {
    'P06': 'equivalent: when the predecessor of j is the empty suffix, PLCP[j] = 0, and the Kasai invariant PLCP[j] >= PLCP[j-1] - 1 then gives PLCP[j-1] <= 1, so the carried length (PLCP[j-1] - 1, at least 0) is already 0 at j; at the first position of a chunk l starts at 0 anyway',
    'B29': 'equivalent (time only): a worker that does not look at the error flag goes on with the remaining chunks; the first error is still stored and rethrown after all tasks have ended, scratch is removed by the spool disposal; no observable result changes',
}
REPLACED = {'B12': 'B12b', 'B26': 'B26b', 'E07': 'E07b', 'E08': 'E08b', 'E09': 'E09b', 'P09': 'P09b', 'E05': 'H02'}
rows = []
r1, r2, files, descs = {}, {}, {}, {}
for jf in sorted(glob.glob(os.path.join(HERE, 'e4_*_mutants.json'))):
    base = os.path.basename(jf)
    round2 = base.startswith('e4_round2_')
    for m in json.load(open(jf)):
        (files if not round2 else files)[m['id']] = files.get(m['id'], base)
        descs[m['id']] = m['desc']
for rf in sorted(glob.glob(os.path.join(HERE, 'results_e4_*.tsv'))):
    if os.path.basename(rf) == 'results_e4.tsv':
        continue
    target = r2 if os.path.basename(rf).startswith('results_e4_round2_') else r1
    for line in open(rf):
        c = line.rstrip('\n').split('\t')
        if len(c) >= 3:
            target[c[0]] = (c[2], c[3] if len(c) > 3 else '')
src = {}
for jf in sorted(glob.glob(os.path.join(HERE, 'e4_*_mutants.json'))):
    name = os.path.basename(jf)[3:-13]
    for m in json.load(open(jf)):
        src.setdefault(m['id'], name)
out = ['id\tfile\tdescription\tround 1 (stage 1)\tround 1 (stage 2)\tround 2\tverdict\tnote']
ids = [i for i in r1] + [i for i in r2 if i not in r1]
for i in ids:
    a = r1.get(i, ('', ''))
    b = r2.get(i, ('', ''))
    killed = a[0].startswith('killed') or a[1].startswith('killed') or b[0].startswith('killed') or b[1].startswith('killed')
    note = ''
    if killed:
        verdict = 'killed'
    elif i in EQUIVALENT:
        verdict, note = 'equivalent', EQUIVALENT[i]
    elif a[0].startswith('BUILD-ERROR') and i in REPLACED:
        verdict, note = 'not compiled', 'replaced by ' + REPLACED[i]
    elif a[0].startswith('BUILD-ERROR'):
        verdict = 'not compiled'
    elif b[0].startswith('BUILD-ERROR'):
        verdict = 'not compiled'
    else:
        verdict = 'SURVIVOR'
    out.append('\t'.join([i, src.get(i, ''), descs.get(i, ''), a[0], a[1], b[0] + ((' ' + b[1]) if b[1] else ''), verdict, note]))
open(os.path.join(HERE, 'results_e4.tsv'), 'w').write('\n'.join(out) + '\n')
v = [l.split('\t')[6] for l in out[1:]]
print({k: v.count(k) for k in set(v)}, len(v), 'mutants')
