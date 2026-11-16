#!/usr/bin/env python3
# Final status of the E9 mutation audit from the three rounds of results (results_e9_<file>.tsv round 1, results_e9_round2_<file>.tsv, results_e9_round3_<file>.tsv): one status per mutant id and file.
# A mutant that survived is "killed after round N" when a later round (run after the tests were strengthened) killed it; the equivalent ones carry the justification below; the rest are reported as survivors.
# usage: python3 -I tables_e9_mutation.py
import glob, os

H = os.path.dirname(os.path.abspath(__file__))
EQUIV = {
    ('PsiStream', 'P20'): 'equivalent: the group start pass * SymbolsPerPass and pass * group differ only when sigma < SymbolsPerPass, and then there is one pass (pass = 0)',
    ('BwtPass', 'B08'): 'equivalent: row 0 must hold n (checked), so the first row that holds 0 is row >= 1 (for n = 0 there is one row): "primary >= 1" and "primary >= 0" agree whenever a second row holds 0',
    ('IsaByDist', 'D07'): 'equivalent: rows entries in [0, n] without a missing position are a permutation; a duplicate always leaves a position missing, which the missing-position check refuses (same exception type; only the message differs)',
    ('IsaByDist', 'D08'): 'equivalent: the converse: with the duplicate check in place, no duplicate and rows entries in [0, n] leave nothing missing',
}
NOTE = {
    ('LinkSmoke', 'S04'): 'SURVIVOR, not killed: an exception leaves okB = false, which already differs from a successful reference run; the term matters only when the reference run gave up and the other run threw, which the tests could not construct',
    ('LinkSmoke', 'S02'): 'SURVIVOR, not killed: the step counts are compared with <= instead of !=; the tables the tests corrupt never give a link that differs only by more steps in the table under test',
}
FILES = {'psistream': 'PsiStream', 'bwtpass': 'BwtPass', 'isabydist': 'IsaByDist', 'psifiles': 'PsiFiles', 'linksmoke': 'LinkSmoke'}

def load(pattern):
    out = {}
    for path in sorted(glob.glob(os.path.join(H, pattern))):
        base = os.path.basename(path)[len('results_e9_'):-4]
        for prefix in ('round2_', 'round3_'):
            if base.startswith(prefix):
                base = base[len(prefix):]
        for line in open(path):
            f = line.rstrip('\n').split('\t')
            if len(f) >= 3 and f[0]:
                out[(FILES[base], f[0])] = (f[1], f[2], f[3] if len(f) > 3 else '')
    return out

r1 = {}
for p in glob.glob(os.path.join(H, 'results_e9_*.tsv')):
    if '_round' in p:
        continue
    r1.update(load(os.path.basename(p)))
r2 = {}
for p in glob.glob(os.path.join(H, 'results_e9_round2_*.tsv')):
    r2.update(load(os.path.basename(p)))
r3 = {}
for p in glob.glob(os.path.join(H, 'results_e9_round3_*.tsv')):
    r3.update(load(os.path.basename(p)))

# round 2 replaces the mutants that did not compile in round 1 (L01, L02, S02, S06) and re-runs survivors: its entry for the same id is the later run
final = {}
for key in sorted(set(r1) | set(r2)):
    a, b, c = r1.get(key), r2.get(key), r3.get(key)
    steps = [(n, x) for n, x in (('1', a), ('2', b), ('3', c)) if x]
    status = steps[-1][1][1]
    first = steps[0][1][1]
    final[key] = (steps, status)

def kind(s):
    return 'killed' if s.startswith('killed') else s.split(' ')[0] if s.startswith('BUILD') else s

killed = equiv = surv = nocompile_only = stale = 0
rows = []
for key, (steps, status) in sorted(final.items()):
    desc = (r2.get(key) or r1.get(key))[0]
    st = status
    if key in EQUIV:
        verdict = 'equivalent (round 1: survived): ' + EQUIV[key].split(': ', 1)[1]
        equiv += 1
    elif status.startswith('killed'):
        verdict = 'killed' + ('' if len(steps) == 1 else ' (round %s, after round %s: %s)' % (steps[-1][0], steps[-2][0], steps[-2][1][1].split(' ')[0]))
        killed += 1
    elif key in NOTE:
        verdict = NOTE[key]
        surv += 1
    elif status == 'PATTERN-NOT-UNIQUE' or status.startswith('BUILD'):
        verdict = 'not run: ' + status
        stale += 1
    else:
        verdict = 'SURVIVOR: ' + status
        surv += 1
    rows.append((key[0], key[1], desc, verdict))

print('| file | id | mutant | final status |')
print('|---|---|---|---|')
for r in rows:
    print('| %s | %s | %s | %s |' % r)
print()
print('distinct mutants with a final status: %d; killed %d; equivalent (justified) %d; survivors not killed %d; never run (no compile / stale pattern, not replaced) %d' % (len(rows), killed, equiv, surv, stale))
