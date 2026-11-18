#!/usr/bin/env python3
"""Merges the four rounds of the verifier audit into one table: the last result of every mutant id, with the argued equivalents marked. usage: summarize_e8_mutants.py > results_e8_final.tsv"""
import glob, json, re
EQUIVALENT = {
    'B06': 'a count that is too high is always accompanied by one that is too low (both sums are n), so "<" and "!=" report the same inputs',
    'C13': 'the symbols at the mismatch differ, so "<=" and "<" give the same order',
    'E06': 'a counts list that is too long also fails the sum check of the metadata (the sum of the counts is n)',
    'E10': 'sigma 257 needs 9 bits, which the bits clause (at most 8) and the bits-for-sigma check refuse',
    'F06': 'the escape of the 4-byte base is read only when the base has an escape (it has none)',
    'D14': 'vacuous: string concatenation with the empty string (written by mistake, not a mutant)',
}
SURVIVOR = {'C19d': 'defence in depth: no input was found in which the spread probes are the only check that fails (a corruption in the middle of a long match is also seen by the exact checks of the short pairs around it)'}
desc = {}
for f in glob.glob('e8_*_listed.json'):
    for m in json.load(open(f)):
        desc.setdefault(m['id'], m['desc'])
last = {}
for rnd in (1, 2, 3, 4):
    for f in sorted(glob.glob('results_round%d_*.tsv' % rnd)):
        for line in open(f):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3:
                last[p[0]] = (rnd, p[2], p[3] if len(p) > 3 else '', p[1])
rows = []
for i in sorted(last):
    rnd, r1, r2, d = last[i]
    if i in EQUIVALENT: status = 'equivalent (argued): ' + EQUIVALENT[i]
    elif i in SURVIVOR: status = 'survivor: ' + SURVIVOR[i]
    elif r1.startswith('killed') or r2.startswith('killed'): status = 'killed'
    elif r1.startswith('BUILD-ERROR') or r1.startswith('PATTERN'): status = r1
    else: status = 'SURVIVED (not argued)'
    rows.append((i, rnd, d, status))
print('id\tlast round\tmutant\tstatus')
for r in rows: print('\t'.join(map(str, r)))
import collections
c = collections.Counter(s.split(':')[0].split(' (')[0] for _, _, _, s in rows)
import sys
print('# %d mutant ids; %s' % (len(rows), ', '.join('%s %d' % kv for kv in sorted(c.items()))), file=sys.stderr)
