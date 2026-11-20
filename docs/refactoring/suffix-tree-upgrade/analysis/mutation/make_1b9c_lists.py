#!/usr/bin/env python3
"""Work package 1b.9c: the lists of the mutants that the audits of 1b.7d and 1b.7e did not run, one list per source file (the runner mutates one file at a time).
   1b.7d: the full lists 1b7d_<file>_listed.json (561 mutants) minus the mutants that summary_1b7d.tsv holds (285, the even positions): 276.
   1b.7e: the lists 1b7e_r1_* and 1b7e_r2_* (600 mutants) minus the mutants that summary_1b7e.tsv holds with a verdict: the lines 'not run' (133).
   A mutant whose pattern no longer matches the production code at the head exactly once (the code was edited by 1b.7e, 1b.9a or 1b.9b after the list was written) is written again in 1b9c_fixups.py (id -> old, new).
   Writes 1b9c_arith.json (the arithmetic per file, read by summarize_1b9c.py) and the lists 1b9c_<d|e>_<file>_listed.json.
   usage: make_1b9c_lists.py [--arith-only]"""
import collections, glob, json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '../../../../..'))
exec(open(os.path.join(HERE, '1b9c_fixups.py')).read())  # FIXUPS = {('d'|'e', id): {'old':..., 'new':..., 'desc'?:...}}


def lists(pre):
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, pre + '_*_listed.json'))):
        b = os.path.basename(f)
        if pre == '1b7d' and re.match(r'1b7d_r\d', b):
            continue
        if pre == '1b7e' and not re.match(r'1b7e_r[12]_', b):
            continue
        out += json.load(open(f))
    return out


def src(m):
    return open(os.path.join(ROOT, m['file'])).read()


arith = collections.OrderedDict()
by_file = collections.OrderedDict()
# 1b.7d
status_d = {}
for l in open(os.path.join(HERE, 'summary_1b7d.tsv')):
    p = l.rstrip('\n').split('\t')
    status_d[(p[0], p[1])] = p[2]
listed_d = collections.OrderedDict()
for f in sorted(glob.glob(os.path.join(HERE, '1b7d_*_listed.json'))):
    b = os.path.basename(f)
    if re.match(r'1b7d_r\d', b):
        continue
    base = b[len('1b7d_'):-len('_listed.json')]
    listed_d[base] = json.load(open(f))
for base, l in listed_d.items():
    c = collections.Counter(status_d.get((base, m['id']), 'unrun') for m in l)
    arith['1b7d/' + base] = {'listed': len(l), **c}
    for m in l:
        if (base, m['id']) not in status_d:
            by_file.setdefault(('d', m['file']), []).append(dict(m, id=m['id'], base=base))
# 1b.7e
status_e = {}
for l in open(os.path.join(HERE, 'summary_1b7e.tsv')):
    p = l.rstrip('\n').split('\t')
    status_e[p[0]] = p[1]
per_file_e = collections.OrderedDict()
for m in lists('1b7e'):
    per_file_e.setdefault(os.path.basename(m['file']), []).append(m)
for fn, l in per_file_e.items():
    c = collections.Counter('unrun' if status_e[m['id']] == 'not run' else status_e[m['id']] for m in l)
    arith['1b7e/' + fn] = {'listed': len(l), **c}
    for m in l:
        if status_e[m['id']] == 'not run':
            by_file.setdefault(('e', m['file']), []).append(dict(m, base=fn))
json.dump(arith, open(os.path.join(HERE, '1b9c_arith.json'), 'w'), indent=1)
if '--arith-only' in sys.argv:
    sys.exit(0)
n = 0
for (k, file), ms in by_file.items():
    out = []
    for m in ms:
        fx = FIXUPS.get((k, m['id']))
        if fx:
            m = dict(m, old=fx['old'], new=fx['new'])
        if src(m).count(m['old']) != 1:
            print('STALE', k, m['id'], src(m).count(m['old']), m['desc'])
        out.append({'id': m['id'], 'desc': m['desc'], 'file': m['file'], 'old': m['old'], 'new': m['new'], 'filter': m.get('filter', ''), 'orig': ('d_' if k == 'd' else 'e_') + m['id']})
    name = '1b9c_%s_%s_listed.json' % (k, os.path.splitext(os.path.basename(file))[0].lower())
    json.dump(out, open(os.path.join(HERE, name), 'w'), indent=1)
    n += len(out)
    print(name, len(out))
print('total', n)
