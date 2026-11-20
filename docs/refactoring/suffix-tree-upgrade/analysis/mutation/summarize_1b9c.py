#!/usr/bin/env python3
"""Final status of every listed mutant of 1b.7d (561) and 1b.7e (600) after work package 1b.9c, over the earlier audits (summary_1b7d.tsv, summary_1b7e.tsv) and the runs of this package:
     results_1b9c_<d|e>_<file>.tsv      stage 1: the 409 mutants that the earlier audits did not run (make_1b9c_lists.py), the tests of the list and the classes that name the file
     results_1b9c_s2_<d|e>_<file>.tsv   stage 2: the survivors of stage 1 against the whole default set of SuffixTree.Esa.Tests, at the head of the audit (before the new tests)
     results_1b9c_hang.tsv              the one-by-one re-run of every kill that only a watchdog or a crash of the host counted (the 17 of the earlier audits and those of this package)
     results_1b9c_c_<d|e>_<file>.tsv    the mutants that a new test closes, run against the final code and the new test classes only: the kill must reproduce
   Statuses (the order is the order of the rules):
     killed (earlier audit)          killed by a failing assertion in a run of 1b.7d / 1b.7e
     killed (1b.9c, stage 1 / 2)     killed by a failing assertion in a run of this package
     killed, assertion (re-run)      a kill by a watchdog or a crash that the one-by-one re-run found to be a failing assertion on its own
     killed, hang / crash (re-run)   a kill that the re-run confirmed as a hang or a crash of the host (legitimate: the mutant does not terminate or kills the process); see ADR-004 section 1b.9c for the table
     closed by a test                a survivor of stage 2 that a new test closes (1b9c_closed.py) and whose re-run against the final tests was killed by a test of this package
     equivalent                      argued in 1b7d_equivalent.py / 1b7e_equivalent.py (EQUIVALENT earlier, EQUIVALENT_1B9C in this package)
     OPEN                            anything else; the gate is zero
   Writes summary_1b9c.tsv (package, id, status, description) and prints the arithmetic per file and in total.
   usage: summarize_1b9c.py"""
import collections, glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, '1b9c_closed.py')).read())  # CLOSED
ns = {'d': {}, 'e': {}}
exec(open(os.path.join(HERE, '1b7d_equivalent.py')).read(), ns['d'])
exec(open(os.path.join(HERE, '1b7e_equivalent.py')).read(), ns['e'])
EQ = {k: dict(ns[k]['EQUIVALENT'], **ns[k]['EQUIVALENT_1B9C']) for k in 'de'}
DEFENCE = {k: ns[k]['DEFENCE'] for k in 'de'}

earlier = {}
for l in open(os.path.join(HERE, 'summary_1b7d.tsv')):
    b, i, st, d = l.rstrip('\n').split('\t')
    earlier[('d', i)] = st
for l in open(os.path.join(HERE, 'summary_1b7e.tsv')):
    i, st, d = l.rstrip('\n').split('\t')
    earlier[('e', i)] = st
listed = collections.OrderedDict()
for f in sorted(glob.glob(os.path.join(HERE, '1b7d_*_listed.json'))):
    b = os.path.basename(f)
    if not re.match(r'1b7d_r\d', b):
        for m in json.load(open(f)):
            listed.setdefault('1b7d/' + b[len('1b7d_'):-len('_listed.json')], []).append(('d', m['id'], m['desc']))
tmp = collections.OrderedDict()
for f in sorted(glob.glob(os.path.join(HERE, '1b7e_r[12]_*_listed.json'))):
    for m in json.load(open(f)):
        tmp.setdefault('1b7e/' + os.path.basename(m['file']), []).append(('e', m['id'], m['desc']))
listed.update(tmp)


def read(pattern, rx):
    out = collections.defaultdict(list)
    for f in sorted(glob.glob(os.path.join(HERE, pattern))):
        m = re.match(rx, os.path.basename(f))
        if not m:
            continue
        for line in open(f):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3:
                out[(m.group(1), p[0])].append(p[2:])
    return out


s1 = read('results_1b9c_[de]_*.tsv', r'results_1b9c_([de])_.*\.tsv')
s2 = read('results_1b9c_s2_[de]_*.tsv', r'results_1b9c_s2_([de])_.*\.tsv')
cl = read('results_1b9c_c_[de]_*.tsv', r'results_1b9c_c_([de])_.*\.tsv')
hang = {}
if os.path.exists(os.path.join(HERE, 'results_1b9c_hang.tsv')):
    for line in open(os.path.join(HERE, 'results_1b9c_hang.tsv')):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3:
            k, i = p[0].split(':')
            hang[(k, i)] = p[2:]


def assertion(rs):
    return any(re.search(r'killed\(\d+ failed', ' '.join(r)) for r in rs)


def dog(rs):
    return any(('killed(hang' in ' '.join(r) or 'killed(host crash' in ' '.join(r)) for r in rs)


def survived(rs):
    return bool(rs) and all(r[0] == 'survived' or r[0].startswith('NO-RESULT') for r in rs) and not dog(rs)


out = open(os.path.join(HERE, 'summary_1b9c.tsv'), 'w')
tally = collections.OrderedDict()
open_ids = []
for fname, ms in listed.items():
    c = collections.Counter()
    for k, i, d in ms:
        key = (k, i)
        st0 = earlier.get(key, 'not run')
        status = None
        if st0 == 'killed':
            status = 'killed (earlier audit)'
        elif st0 == 'not run':
            if assertion(s1.get(key, [])):
                status = 'killed (1b.9c, stage 1)'
            elif key in hang and hang[key][0] == 'assertion':
                status = 'killed, assertion (re-run)'
            elif dog(s1.get(key, [])) or st0 == 'killed by the watchdog':
                status = 'killed, hang / crash (re-run)' if key in hang and hang[key][0] in ('hang', 'crash') else 'OPEN (hang or crash not re-run, or survived the re-run)'
            elif s1.get(key) and assertion(s2.get(key, [])):
                status = 'killed (1b.9c, stage 2)'
        elif st0 == 'killed by the watchdog':
            status = 'killed, hang / crash (re-run)' if key in hang and hang[key][0] in ('hang', 'crash') else ('killed, assertion (re-run)' if key in hang and hang[key][0] == 'assertion' else 'OPEN (hang or crash not re-run)')
        if status is None and key in CLOSED and assertion(cl.get(key, [])):
            status = 'closed by a test'
        if status is None and i in EQ[k] and (st0 in ('not run', 'equivalent', 'defence') or key in CLOSED or True):
            if st0 != 'killed':
                status = 'equivalent'
        if status is None and st0 == 'equivalent':
            status = 'equivalent'
        if status is None:
            status = 'OPEN (%s)' % ('closed by a test, not yet reproduced' if key in CLOSED else 'survivor without a closure' if s1.get(key) else 'not run')
        if hang.get((k, i)) and hang[(k, i)][0] == 'assertion' and st0 == 'killed by the watchdog':
            status = 'killed, assertion (re-run)'
        c[status.split(' (')[0] if status.startswith('OPEN') else status] += 1
        out.write('%s\t%s\t%s\t%s\n' % (k, i, status, d))
        if status.startswith('OPEN'):
            open_ids.append((k, i, status))
    tally[fname] = c
out.close()
total = collections.Counter()
for fname, c in tally.items():
    total.update(c)
    print('%-34s %4d  %s' % (fname, sum(c.values()), dict(c)))
print('TOTAL', sum(total.values()), dict(total))
print('open:', len(open_ids), open_ids[:20])
