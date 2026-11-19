#!/usr/bin/env python3
"""Final status of every mutant of the audit of 1b.7e over all rounds (results_1b7e_<round>_*.tsv; the round is the prefix r1_ .. r8_ of the file name): r1 and r2 are the two halves of the list (even and odd positions) with the narrow
filters of each file, r3 .. r8 are the reruns (survivors after the tests that closed them, and the runs that ended by the watchdog, which are run again with a longer limit and without StopOnError, because a test host that stops at the first
failure can hang). A mutant is KILLED when a run ended with a failing test (killed(N), N >= 1), KILLED BY THE WATCHDOG when only runs that hit the time limit or crashed the host say so (an infinite loop of the mutant, a test host that
does not end), EQUIVALENT / DEFENCE when no run killed it and 1b7e_equivalent.py argues it, OPEN otherwise; NOT RUN when it is in the lists and no round ran it. Writes summary_1b7e.tsv."""
import glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, '1b7e_equivalent.py')).read())
listed = {}
for f in glob.glob(os.path.join(HERE, '1b7e_r[12]_*_listed.json')):
    for m in json.load(open(f)):
        listed[m['id']] = m['desc']
desc, runs = {}, {}
for f in sorted(glob.glob(os.path.join(HERE, 'results_1b7e_r*_*.tsv'))):
    name = os.path.basename(f)[len('results_1b7e_'):-4]
    m = re.match(r'(r\d+)_(.*)', name)
    if not m:
        continue
    rnd = m.group(1)
    for line in open(f):
        p = line.rstrip('\n').split('\t')
        if len(p) < 3:
            continue
        desc[p[0]] = p[1]
        runs.setdefault(p[0], []).append((rnd, p[2], p[3] if len(p) > 3 else ''))
final = {}
for i, lst in runs.items():
    states = [s for _, a, b in lst for s in (a, b)]
    clean = any(re.fullmatch(r'killed\(\d+\)', s) for s in states)
    # the runs of r1 .. r6 stop at the first failing test (StopOnError), and a test host that has to stop can hang: a time limit there is not evidence; from r7 on the runs do not stop at the first failure
    late = [(a, b) for r, a, b in lst if int(r[1:]) >= 7]
    late_dog = any(a.startswith('killed(hang') or a.startswith('killed(host') for a, _ in late)
    early_dog = any(a.startswith('killed(hang') or a.startswith('killed(host') or b.startswith('killed(hang') for r, a, b in lst if int(r[1:]) < 7)
    built = all(a.startswith('BUILD') for _, a, _ in lst)
    if clean:
        final[i] = 'killed'
    elif late_dog:
        final[i] = 'killed by the watchdog'
    elif late:
        final[i] = 'equivalent' if i in EQUIVALENT else ('defence' if i in DEFENCE else 'OPEN')
    elif early_dog:
        final[i] = 'killed by the watchdog (unconfirmed)'
    elif built:
        final[i] = 'not compiled'
    elif i in EQUIVALENT:
        final[i] = 'equivalent'
    elif i in DEFENCE:
        final[i] = 'defence'
    else:
        final[i] = 'OPEN'
tally = {}
for o in final.values():
    tally[o] = tally.get(o, 0) + 1
notrun = sorted(set(listed) - set(final))
with open(os.path.join(HERE, 'summary_1b7e.tsv'), 'w') as out:
    for i, o in sorted(final.items()):
        out.write('%s\t%s\t%s\n' % (i, o, desc[i]))
    for i in notrun:
        out.write('%s\tnot run\t%s\n' % (i, listed[i]))
second = [i for i, lst in runs.items() if final[i] == 'killed' and not any(re.fullmatch(r'killed\(\d+\)', a) for _, a, _ in lst) and any(re.fullmatch(r'killed\(\d+\)', b) for _, _, b in lst)]
print('listed', len(listed), 'ran', len(final), 'not run', len(notrun))
print('status', json.dumps(tally, sort_keys=True))
print('killed in stage 2 only (second line):', sorted(second))
print('open:', sorted(i for i, o in final.items() if o == 'OPEN'))
print('equivalent argued:', len(EQUIVALENT), 'defence argued:', len(DEFENCE))
print('argued but killed or not run:', sorted(i for i in list(EQUIVALENT) + list(DEFENCE) if final.get(i) not in ('equivalent', 'defence')))
