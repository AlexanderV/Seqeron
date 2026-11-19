#!/usr/bin/env python3
"""Summary of the audit of work package 1b.8: reads the results of every round (results_1b8_*.tsv), the arguments (1b8_equivalent.py) and the sidecar files of the runs (results_*.tsv.<id>.rep.log: the progress log of the harness
   with the first failure and its minimized text), writes summary_1b8.tsv (one row per listed mutant: file, description, round-1 outcome, final status with its evidence) and fails_1b8.tsv (the first failure line of every kill by
   the harness, in full), and with --clean removes the sidecars (they are folded into fails_1b8.tsv; nothing else is kept of a run).
   usage: summarize_1b8.py [--clean]"""
import glob, importlib.util, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('eq', os.path.join(HERE, '1b8_equivalent.py'))
eq = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eq)

listed = {}
for f in sorted(glob.glob(os.path.join(HERE, '1b8_*_listed.json'))):
    name = os.path.basename(f)
    if re.match(r'1b8_r\d+u?_', name):
        continue
    for m in json.load(open(f)):
        listed[m['id']] = (os.path.basename(m['file']), m['desc'])


def read(pattern):
    out = {}
    for f in sorted(glob.glob(os.path.join(HERE, pattern))):
        for line in open(f, errors='replace'):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3:
                out[p[0]] = p[2]
    return out


round1 = {}
for f in sorted(glob.glob(os.path.join(HERE, 'results_1b8_*.tsv'))):
    if re.match(r'results_1b8_r\d+u?_', os.path.basename(f)):
        continue
    for line in open(f, errors='replace'):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3:
            round1[p[0]] = p[2]
later = {}
for r in ('r2', 'r3', 'r4'):
    for k, v in read('results_1b8_%s_*.tsv' % r).items():
        later[k] = (r, v)
unit = read('results_1b8_r3u_*.tsv')

fails = {}
for f in sorted(glob.glob(os.path.join(HERE, 'results_*.tsv.*.rep.log'))):
    mid = re.search(r'\.tsv\.(1b8_\w+)\.rep\.log$', f).group(1)
    for line in open(f, errors='replace'):
        if line.startswith('FAIL '):
            fails.setdefault(mid, []).append(line.strip())

counts = {}
rows = []
for mid, (fname, desc) in listed.items():
    r1 = round1.get(mid, 'NOT RUN')
    if r1.startswith('killed'):
        status = 'killed by the harness in round 1'
    elif mid in later and later[mid][1].startswith('killed'):
        status = 'killed by the harness after an extension (round %s)' % later[mid][0][1:]
    elif mid in unit and unit[mid].startswith('killed'):
        status = 'killed by the existing unit tests, not by the harness: ' + eq.KILLED_BY_UNIT_TESTS.get(mid, '')
    elif mid in eq.EQUIVALENT:
        status = 'argued equivalent: ' + eq.EQUIVALENT[mid]
    elif mid in eq.DEFENCE:
        status = 'second-line defence: ' + eq.DEFENCE[mid]
    else:
        status = 'OPEN (' + r1 + ')'
    key = status.split(':')[0].split(' (')[0]
    counts[key] = counts.get(key, 0) + 1
    rows.append((mid, fname, desc, r1[:60], status))

with open(os.path.join(HERE, 'summary_1b8.tsv'), 'w') as o:
    o.write('id\tfile\tdescription\tround 1\tfinal status\n')
    for r in rows:
        o.write('\t'.join(x.replace('\t', ' ').replace('\n', ' ') for x in r) + '\n')
with open(os.path.join(HERE, 'fails_1b8.tsv'), 'w') as o:
    o.write('id\tfirst failure of the harness (the minimized text as hex code units at the end of the line)\n')
    for mid in sorted(fails):
        o.write('%s\t%s\n' % (mid, fails[mid][0]))
print('listed', len(listed), 'rows', len(rows))
for k, v in sorted(counts.items()):
    print(v, k)
if '--clean' in sys.argv:
    n = 0
    for f in glob.glob(os.path.join(HERE, 'results_*.tsv.*.rep*')) + glob.glob(os.path.join(HERE, 'results_1b8_*.log')) + glob.glob(os.path.join(HERE, '1b8_*audit.done')):
        os.remove(f)
        n += 1
    print('removed', n, 'sidecar files')
