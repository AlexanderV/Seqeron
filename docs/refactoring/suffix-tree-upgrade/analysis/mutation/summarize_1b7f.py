#!/usr/bin/env python3
"""Final status of every mutant of the audit of 1b.7f over all rounds (results_1b7f_<file>_<round>.tsv; the round is the suffix of the file name: '' for the first run, r2, r3, ...). A mutant is KILLED when a run ended with a
failing test (killed(N), N >= 1), KILLED BY THE WATCHDOG when only runs that hit the time limit or crashed the host say so (an infinite loop of the mutant), BUILD-ERROR when no run compiled it (not counted as a mutant that ran),
EQUIVALENT / DEFENCE when no run killed it and 1b7f_equivalent.py argues it, OPEN otherwise; NOT RUN when it is in the lists and no round ran it. Writes summary_1b7f.tsv and prints the counts."""
import glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, '1b7f_equivalent.py')).read())
listed = {}
for f in glob.glob(os.path.join(HERE, '1b7f_*_listed.json')):
    for m in json.load(open(f)):
        listed[m['id']] = m['desc']
runs, desc = {}, {}
for f in sorted(glob.glob(os.path.join(HERE, 'results_1b7f_*.tsv'))):
    name = os.path.basename(f)[len('results_1b7f_'):-4]
    rnd = name.rsplit('_', 1)[1] if re.search(r'_(r\d+)?$', name) else ''
    for line in open(f):
        p = line.rstrip('\n').split('\t')
        if len(p) < 3:
            continue
        desc[p[0]] = p[1]
        runs.setdefault(p[0], []).append((name, p[2]))
final = {}
for i in listed:
    lst = runs.get(i, [])
    states = [s for _, s in lst]
    if not lst:
        final[i] = 'NOT RUN'
    elif any(re.fullmatch(r'killed\(\d+\)', s) for s in states):
        final[i] = 'killed'
    elif any(s.startswith('killed(hang') or s.startswith('killed(host') for s in states):
        final[i] = 'killed by the watchdog'
    elif all(s.startswith('BUILD') or s.startswith('PATTERN') for s in states):
        final[i] = 'build error'
    elif i in EQUIVALENT:
        final[i] = 'equivalent'
    elif i in DEFENCE:
        final[i] = 'defence'
    else:
        final[i] = 'OPEN'
with open(os.path.join(HERE, 'summary_1b7f.tsv'), 'w') as out:
    for i in sorted(listed, key=lambda x: (re.sub(r'\d+', '', x), int(re.sub(r'\D', '', x) or 0))):
        out.write('%s\t%s\t%s\t%s\n' % (i, final[i], listed[i], '; '.join('%s=%s' % (n, s) for n, s in runs.get(i, []))))
from collections import Counter
c = Counter(final.values())
print(len(listed), 'listed:', dict(c))
for i, st in final.items():
    if st == 'OPEN':
        print('OPEN', i, listed[i])
