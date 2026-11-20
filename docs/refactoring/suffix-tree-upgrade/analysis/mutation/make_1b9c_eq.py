#!/usr/bin/env python3
"""Work package 1b.9c: the lists for the check of the argued equivalents (EQUIVALENT_1B9C of 1b7d_equivalent.py and 1b7e_equivalent.py): an equivalent mutant must survive the whole default set of SuffixTree.Esa.Tests as well, so each of them is run
   with the stage 2 of the audit (run_1b9c_audit.sh with MUTANT_STAGES=2) at the final head. One list per source file: 1b9c_q_<d|e>_<file>_listed.json, results_1b9c_q_<d|e>_<file>.tsv. A pattern that the seams of this package moved is read from
   1b9c_fixups.py. usage: make_1b9c_eq.py"""
import glob, json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '../../../../..'))
exec(open(os.path.join(HERE, '1b9c_fixups.py')).read())
ns = {'d': {}, 'e': {}}
exec(open(os.path.join(HERE, '1b7d_equivalent.py')).read(), ns['d'])
exec(open(os.path.join(HERE, '1b7e_equivalent.py')).read(), ns['e'])
pool = {}
for f in sorted(glob.glob(os.path.join(HERE, '1b7d_*_listed.json'))):
    if not re.match(r'1b7d_r\d', os.path.basename(f)):
        for m in json.load(open(f)):
            pool[('d', m['id'])] = m
for f in sorted(glob.glob(os.path.join(HERE, '1b7e_r[12]_*_listed.json'))):
    for m in json.load(open(f)):
        pool[('e', m['id'])] = m
by_file = {}
for k in 'de':
    for i in sorted(ns[k]['EQUIVALENT_1B9C']):
        m = dict(pool[(k, i)])
        if (k, i) in FIXUPS:
            m['old'], m['new'] = FIXUPS[(k, i)]['old'], FIXUPS[(k, i)]['new']
        if open(os.path.join(ROOT, m['file'])).read().count(m['old']) != 1:
            print('STALE', k, i, m['desc'])
        by_file.setdefault((k, m['file']), []).append({'id': i, 'desc': m['desc'], 'file': m['file'], 'old': m['old'], 'new': m['new'], 'filter': ''})
for (k, file), ms in by_file.items():
    name = '1b9c_q_%s_%s_listed.json' % (k, os.path.splitext(os.path.basename(file))[0].lower())
    json.dump(ms, open(os.path.join(HERE, name), 'w'), indent=1)
    print(name, [m['id'] for m in ms])
