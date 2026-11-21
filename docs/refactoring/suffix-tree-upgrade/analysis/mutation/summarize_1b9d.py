#!/usr/bin/env python3
"""Summary of the mutation audit of work package 1b.9d (run_1b9d_audit.py): one line per mutant with its verdict in stage 1 (the tests of the package and the older tests of the same code), stage 2 (the whole default set, for
a survivor of stage 1), the re-run of a kill by the death of the host with the tests that do not crash it, and the closure (a test, or the argument in 1b9d_equivalent.py). Writes summary_1b9d.tsv and prints the counts.
usage: summarize_1b9d.py"""
import glob, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
try:
    from importlib import import_module
    EQ = import_module('1b9d_equivalent').EQUIVALENT
except ImportError:
    EQ = {}
listed = {}
for f in sorted(glob.glob(os.path.join(HERE, '1b9d_*_listed.json'))):
    name = os.path.basename(f)
    if name.startswith(('1b9d_s2_', '1b9d_h_', '1b9d_c_')):
        continue
    for m in json.load(open(f)):
        listed[m['id']] = m

def results(prefix):
    out = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'results_1b9d_%s*.tsv' % prefix))):
        base = os.path.basename(f)
        if prefix == '' and base.startswith(('results_1b9d_s2_', 'results_1b9d_h_', 'results_1b9d_c_')):
            continue
        for line in open(f):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3:
                out[p[0]] = p
    return out

s1, s2, h, cl = results(''), results('s2_'), results('h_'), results('c_')
rows = []
counts = {}
for i, m in sorted(listed.items(), key=lambda kv: (kv[1]['file'], kv[0])):
    v1 = s1.get(i, ['', '', 'NOT RUN'])[2]
    v2 = s2.get(i, ['', '', '', ''])[3] if i in s2 else ''
    vh = h.get(i, ['', '', ''])[2] if i in h else ''
    if v1 == 'NOT RUN':
        status = 'not run'
    elif v1.startswith('killed') and 'host crash' not in v1 and 'hang' not in v1:
        status = 'killed (stage 1, assertion)'
    elif v1.startswith('killed'):
        status = 'killed by the death or stop of the host in stage 1; re-run: ' + (vh.split('[')[0] if vh else 'not re-run')
        if vh.startswith('killed') and 'host crash' not in vh and 'hang' not in vh:
            status = 'killed (stage 1, host death; also by an assertion alone)'
        elif vh.startswith('killed'):
            status = 'killed (stage 1, host death; the tests that do not crash it also see it only as a host death)'
    elif v1.startswith('BUILD'):
        status = 'does not compile'
    else:  # survived stage 1
        vc = cl.get(i, ['', '', ''])[2] if i in cl else ''
        if vc.startswith('killed'):
            status = 'closed by a test (red against the mutant at the head that holds the test)'
        elif v2.startswith('killed'):
            status = 'killed (stage 2)'
        elif i in EQ:
            status = 'equivalent (argued)'
        elif v2 == '':
            status = 'survived stage 1, stage 2 not run'
        else:
            status = 'OPEN'
    counts[status.split(';')[0].split(' (')[0] if status.startswith('killed') else status] = counts.get(status.split(';')[0].split(' (')[0] if status.startswith('killed') else status, 0) + 1
    rows.append((i, os.path.basename(m['file']), m['desc'], v1[:90], (cl.get(i, ['', '', ''])[2][:60] if i in cl else '-'), v2[:60], vh[:60], status))
with open(os.path.join(HERE, 'summary_1b9d.tsv'), 'w') as f:
    f.write('id\tfile\tdescription\tstage 1\tclosed run\tstage 2\trerun without the crashing tests\tstatus\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')
print(len(listed), 'listed;', ', '.join('%s: %d' % kv for kv in sorted(counts.items())))
