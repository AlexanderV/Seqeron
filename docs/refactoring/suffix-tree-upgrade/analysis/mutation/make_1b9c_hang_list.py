#!/usr/bin/env python3
"""Work package 1b.9c: the list of the mutants that an earlier audit counted as killed only by a watchdog or a crash of the test host (the lines of results_1b7d_*.tsv and results_1b7e_*.tsv that say killed(hang) or killed(host crash),
   of mutants that no run of those audits killed by a failing assertion), with the edit of each (from the lists of 1b.7d and 1b.7e), for rerun_1b9c_hang.py. Writes 1b9c_h_all_listed.json."""
import glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
pool = {}
for f in sorted(glob.glob(os.path.join(HERE, '1b7d_*_listed.json'))):
    if not re.match(r'1b7d_r\d', os.path.basename(f)):
        for m in json.load(open(f)):
            pool[('d', m['id'], os.path.basename(m['file']))] = m
for f in sorted(glob.glob(os.path.join(HERE, '1b7e_r[12]_*_listed.json'))):
    for m in json.load(open(f)):
        pool[('e', m['id'], os.path.basename(m['file']))] = m
found = {}
for pre, kind in (('results_1b7d_', 'd'), ('results_1b7e_', 'e')):
    runs = {}
    for f in sorted(glob.glob(os.path.join(HERE, pre + 'r*.tsv'))):
        base = re.sub(r'^r\d+[a-z]?_', '', os.path.basename(f)[len(pre):-4])
        for line in open(f):
            p = line.rstrip('\n').split('\t')
            if len(p) < 3:
                continue
            runs.setdefault((p[0], 'all'), []).append((base, ' '.join(p[2:5])))
    for (i, _), lst in runs.items():
        dog = any('killed(hang' in s or 'killed(host crash' in s for _, s in lst)
        clean = any(re.search(r'killed\(\d+\)', s) for _, s in lst)
        if dog and not clean:
            found[(kind, i)] = next(b for b, s in lst if 'killed(hang' in s or 'killed(host crash' in s)
out = []
for (kind, i), base in sorted(found.items()):
    cands = [m for (k, mid, fn), m in pool.items() if k == kind and mid == i and base.replace('-1', '').replace('-2', '').split('_')[-1] in fn.lower()]
    assert len(cands) == 1, (kind, i, base, [c['file'] for c in cands])
    m = dict(cands[0])
    m['package'] = kind
    out.append(m)
n_earlier = len(out)
# the kills of this package's own runs (results_1b9c_<d|e>_<file>.tsv, stage 1) that only a watchdog or a crash counted
for f in sorted(glob.glob(os.path.join(HERE, '1b9c_[de]_*_listed.json'))):
    name = os.path.basename(f)[len('1b9c_'):-len('_listed.json')]
    res = os.path.join(HERE, 'results_1b9c_%s.tsv' % name)
    if not os.path.exists(res):
        continue
    hung = set()
    for line in open(res):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3 and ('killed(hang' in ' '.join(p[2:5]) or 'killed(host crash' in ' '.join(p[2:5])) and not re.search(r'killed\(\d+ failed', ' '.join(p[2:5])):
            hung.add(p[0])
    for m in json.load(open(f)):
        if m['id'] in hung:
            m = dict(m)
            m['package'] = name[0]
            out.append(m)
json.dump(out, open(os.path.join(HERE, '1b9c_h_all_listed.json'), 'w'), indent=1)
print(len(out), 'mutants:', n_earlier, 'of the earlier audits and', len(out) - n_earlier, 'of this package:', [(m['package'], m['id']) for m in out])
