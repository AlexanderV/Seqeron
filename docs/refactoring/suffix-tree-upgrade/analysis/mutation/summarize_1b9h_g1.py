#!/usr/bin/env python3
"""Gate G1 of work package 1b.9h: the 28 mutants of 1b.9c whose kill was the death or the stop of the test host (14 crashes, 14 hangs), run again by rerun_1b9h_g1.py against the tests of this package.
   usage: summarize_1b9h_g1.py <result.tsv> [...]: one line per mutant (id, what it does, the verdict of 1b.9c, the verdict now, the number of failed tests, the first failing test, seconds) and the counts.
   The verdict 'assertion' means: the run ended (the runner printed its summary), at least one test failed and is named, the host neither crashed nor hung nor was stopped by the limit of the run."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
was = {}
for line in open(os.path.join(HERE, 'results_1b9c_hang.tsv')):
    p = line.rstrip('\n').split('\t')
    was[p[0]] = p[2]
kinds = {}
for name in ('1b9h_g1_crash_listed.json', '1b9h_g1_hang_listed.json'):
    for m in json.load(open(os.path.join(HERE, name))):
        kinds[m['package'] + ':' + m['id']] = ('crash' if 'crash' in name else 'hang', m['desc'])
rows = {}
for f in sys.argv[1:]:
    for line in open(f):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 7:
            rows[p[0]] = p
count = {}
print('mutant\tkind in 1b.9c\twhat it does\tnow\tfailed tests\tfirst failing test\tseconds')
for k in sorted(kinds):
    p = rows.get(k)
    kind, desc = kinds[k]
    now = p[2] if p else 'not run'
    n = p[7] if p and len(p) > 7 else '?'
    first = p[3].split(',')[0] if p and p[3] else ''
    print('\t'.join([k, kind, desc, now, n, first, p[6] if p else '']))
    count[(kind, now)] = count.get((kind, now), 0) + 1
print()
for (kind, now), n in sorted(count.items()):
    print('%s mutants, verdict %s: %d' % (kind, now, n))
