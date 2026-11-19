#!/usr/bin/env python3
"""Final status of every mutant of the audit of 1b.7e over all rounds (results_1b7e_*.tsv; the round is the prefix r1a_, r1b_, r1s_, r2_, ... of the file name, no prefix = the aborted first attempt that is not part of the count): killed in some run, argued equivalent
(1b7e_equivalent.py: EQUIVALENT), kept as a second line of defence that no injectable fault reaches (DEFENCE), or open; also the counts per round. Writes summary_1b7e.tsv."""
import glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, '1b7e_equivalent.py')).read())
desc, rounds = {}, {}
for f in sorted(glob.glob(os.path.join(HERE, 'results_1b7e_*.tsv'))):
    name = os.path.basename(f)[len('results_1b7e_'):-4]
    m = re.match(r'(r\d+[a-z]?)_(.*)', name)
    if not m:
        continue
    rnd, base = m.group(1), m.group(2)
    for line in open(f):
        p = line.rstrip('\n').split('\t')
        if len(p) < 3:
            continue
        desc[(base, p[0])] = p[1]
        rounds.setdefault((base, p[0]), []).append((rnd, p[2], p[3] if len(p) > 3 else ''))
final, per_round = {}, {}
for key, lst in rounds.items():
    outcome = None
    for rnd, st, st2 in lst:
        kind = 'killed' if (st.startswith('killed') or st2.startswith('killed')) else ('not compiled' if st.startswith('BUILD') else 'survived')
        per_round.setdefault(rnd, {}).setdefault(kind, 0)
        per_round[rnd][kind] += 1
        if kind == 'killed':
            outcome = 'killed'
    if outcome is None:
        outcome = 'equivalent' if key[1] in EQUIVALENT else ('defence' if key[1] in DEFENCE else ('not compiled' if all(x[1].startswith('BUILD') for x in lst) else 'OPEN'))
    final[key] = outcome
tally = {}
for o in final.values():
    tally[o] = tally.get(o, 0) + 1
with open(os.path.join(HERE, 'summary_1b7e.tsv'), 'w') as out:
    for (base, i), o in sorted(final.items()):
        out.write('%s\t%s\t%s\t%s\n' % (base, i, o, desc[(base, i)]))
print('per round (runs):', json.dumps(per_round, sort_keys=True))
print('distinct mutants with a result:', len(final), tally)
print('open:', sorted(i for (b, i), o in final.items() if o == 'OPEN'))
print('equivalent ids argued:', len(EQUIVALENT), 'defence:', len(DEFENCE))
