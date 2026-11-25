#!/usr/bin/env python3
"""Summary of the audit of work package S2.4: python3 summarize_s24.py results_s24_audit.tsv [results_s24_r2_audit.tsv]. The last row of an id decides (a later round replaces an earlier one); a survivor is reported as killed(r2) / survived (argued) from s24_equivalent.py."""
import collections, os, re, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from s24_equivalent import EQUIVALENT

final, first, source = {}, {}, {}
for rnd, path in enumerate(sys.argv[1:], 1):
    for line in open(path):
        if not line.strip():
            continue
        mid, desc, verdict = line.rstrip('\n').split('\t')[:3]
        first.setdefault(mid, verdict)
        final[mid] = (desc, verdict)
        source[mid] = rnd
c = collections.Counter()
rows = []
for mid, (desc, verdict) in sorted(final.items()):
    if verdict.startswith('killed'):
        st = re.search(r'stage (\d)', verdict)
        k = 'killed' + (' in stage ' + st.group(1) if st else '')
        round_first = 'killed' if first[mid].startswith('killed') else 'survived round 1'
        rows.append((mid, desc, round_first, 'killed (round %d)' % source[mid]))
    elif mid in EQUIVALENT:
        k = 'equivalent (argued)'
        rows.append((mid, desc, 'survived round 1', 'equivalent'))
    else:
        k = 'SURVIVED (a gap)'
        rows.append((mid, desc, 'survived round 1' if not first[mid].startswith('killed') else 'killed', 'SURVIVED'))
    c[k] += 1
r1 = collections.Counter('killed' if v.startswith('killed') else 'survived' for v in first.values())
print('mutants listed: %d; round 1 (the first verdict of every id): %s' % (len(final), dict(r1)))
for k in sorted(c):
    print('%-26s %d' % (k, c[k]))
print('\nid\tdescription\tround 1\tfinal')
for r in rows:
    if r[2] != 'killed' or not r[3].startswith('killed'):
        print('\t'.join(r))
