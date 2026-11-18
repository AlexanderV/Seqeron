#!/usr/bin/env python3
"""Final status of every mutant of the audit of 1b.7a over the four rounds (results_1b7c_*.tsv): killed in some round, argued equivalent (1b7c_equivalent.py), or open; also the counts per round. Writes summary_1b7c.tsv."""
import glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, '1b7c_equivalent.py')).read())
desc, rounds = {}, {}
for f in sorted(glob.glob(os.path.join(HERE, 'results_1b7c_*.tsv'))):
    name = os.path.basename(f)[len('results_1b7c_'):-4]
    m = re.match(r'(r\d)_(.*)', name)
    rnd = m.group(1) if m else 'r1'
    base = m.group(2) if m else name
    for line in open(f):
        p = line.rstrip('\n').split('\t')
        if len(p) < 3:
            continue
        i = p[0]
        st = p[2]
        st2 = p[3] if len(p) > 3 else ''
        desc[(base, i)] = p[1]
        rounds.setdefault((base, i), []).append((rnd, st, st2))
final = {}
per_round = {}
for key, lst in rounds.items():
    outcome = None
    for rnd, st, st2 in lst:
        kind = 'killed' if (st.startswith('killed') or st2.startswith('killed')) else ('not compiled' if st.startswith('BUILD') else 'survived')
        per_round.setdefault(rnd, {}).setdefault(kind, 0)
        per_round[rnd][kind] += 1
        if kind == 'killed':
            outcome = 'killed'
    base, i = key
    if outcome is None:
        outcome = 'equivalent' if i in EQUIVALENT else ('not compiled, replaced or argued' if all(x[1].startswith('BUILD') for x in lst) else 'OPEN')
    final[key] = outcome
ids = {}
for (base, i), o in final.items():
    ids.setdefault(i, []).append(o)
tally = {}
for i, os_ in ids.items():
    o = 'killed' if 'killed' in os_ else os_[0]
    if o == 'not compiled, replaced or argued' and i in EQUIVALENT:
        o = 'equivalent'
    tally[o] = tally.get(o, 0) + 1
with open(os.path.join(HERE, 'summary_1b7c.tsv'), 'w') as out:
    for (base, i), o in sorted(final.items(), key=lambda x: (x[0][0], x[0][1])):
        out.write('%s\t%s\t%s\t%s\n' % (base, i, o, desc[(base, i)]))
print('per round (runs):', json.dumps(per_round, sort_keys=True))
print('distinct mutant ids with a result:', len(ids), tally)
print('equivalent ids argued:', len(EQUIVALENT), sorted(EQUIVALENT))
