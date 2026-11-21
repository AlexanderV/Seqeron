#!/usr/bin/env python3
"""Summary of the mutation audit of work package 1b.9e.1 (run_1b9e1_audit.py): one line per mutant with its verdict in stage 1 (the quick set for the format, the lock and the reopening, the whole set of the package for
the rest), stage 1b (the survivors of the quick set against the whole set of the package), the re-run at the head that holds the tests written for the survivors (closed), stage 2 (the whole default set of
SuffixTree.Esa.Tests for a survivor), and the closure (a test, or the argument in 1b9e1_equivalent.py). Writes summary_1b9e1.tsv and prints the counts.
usage: summarize_1b9e1.py"""
import glob, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
try:
    from importlib import import_module
    EQ = import_module('1b9e1_equivalent').EQUIVALENT
except ImportError:
    EQ = {}
SKIP = ('1b9e1_s2_', '1b9e1_s2b_', '1b9e1_b_', '1b9e1_c_', '1b9e1_h_')
listed = {}
for f in sorted(glob.glob(os.path.join(HERE, '1b9e1_*_listed.json'))):
    if os.path.basename(f).startswith(SKIP):
        continue
    for m in json.load(open(f)):
        listed[m['id']] = m


def results(prefix):
    out = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'results_1b9e1_%s*.tsv' % prefix))):
        base = os.path.basename(f)
        if prefix == '' and base.startswith(tuple('results_' + p for p in SKIP)):
            continue
        for line in open(f):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3:
                out[p[0]] = p
    return out


def verdict(p):
    return p[2] if p else ''


def is_kill(v):
    return v.startswith('killed')


s1, s1b, s2, s2b, cl, hg = results(''), results('b_'), results('s2_'), results('s2b_'), results('c_'), results('h_')
rows = []
counts = {}
for i, m in sorted(listed.items(), key=lambda kv: (kv[1]['file'], kv[0])):
    v1 = verdict(s1.get(i)) or 'NOT RUN'
    vb = verdict(s1b.get(i))
    vc = verdict(cl.get(i))
    v2 = s2.get(i, ['', '', '', ''])[3] if i in s2 else ''
    if i in s2b:
        v2 = s2b[i][2]  # the stage-2 run without the one allocation test that fails with the load of the machine
    if v1 == 'NOT RUN':
        status = 'not run'
    elif v1.startswith('BUILD'):
        status = 'does not compile'
    elif is_kill(v1) and ('host crash' in v1 or 'hang' in v1) and is_kill(verdict(hg.get(i))):
        status = 'killed (stage 1, by a hang) and again by an assertion (re-run at the head that holds the guarded test)'
    elif is_kill(v1):
        status = 'killed (stage 1)' + (' by the death or the stop of the host' if ('host crash' in v1 or 'hang' in v1) else '')
    elif is_kill(vb):
        status = 'killed (stage 1b, the whole set of the package)'
    elif is_kill(vc):
        status = 'closed by a test (red against the mutant at the head that holds the test)'
    elif v2.startswith('killed'):
        status = 'killed (stage 2)'
    elif i in EQ:
        status = 'equivalent (argued)'
    elif v2 == '':
        status = 'survived, stage 2 not run'
    else:
        status = 'OPEN'
    key = status.split(' by the')[0]
    counts[key] = counts.get(key, 0) + 1
    rows.append((i, os.path.basename(m['file']), m['desc'], v1[:90], vb[:60] or '-', vc[:60] or '-', v2[:60] or '-', status))
with open(os.path.join(HERE, 'summary_1b9e1.tsv'), 'w') as f:
    f.write('id\tfile\tdescription\tstage 1\tstage 1b\tclosed run\tstage 2\tstatus\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')
print(len(listed), 'listed;', ', '.join('%s: %d' % kv for kv in sorted(counts.items())))
