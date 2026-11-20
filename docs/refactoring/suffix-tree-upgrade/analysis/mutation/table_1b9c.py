#!/usr/bin/env python3
"""Work package 1b.9c: the tables of the report in ADR-004 section 5a (1b.9c), from the files of the audit: the arithmetic per file (summary_1b9c.tsv and 1b9c_arith.json), the kills by a watchdog or a crash and their re-run
   (results_1b9c_hang.tsv), the survivors closed by a test (1b9c_closed.py with the result of the re-run against the final tests), the argued equivalents (1b7d_equivalent.py, 1b7e_equivalent.py), the fresh mutants
   (results_1b9c_f_*.tsv). Prints markdown; nothing is written. usage: table_1b9c.py [arith|hang|closed|equivalent|fresh]"""
import collections, glob, json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
what = sys.argv[1] if len(sys.argv) > 1 else 'arith'
exec(open(os.path.join(HERE, '1b9c_closed.py')).read())
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


def short(text, n=170):
    text = text.replace('|', '/').replace('\n', ' ')
    return text if len(text) <= n else text[:n - 3].rsplit(' ', 1)[0] + ' ...'


if what == 'hang':
    print('| mutant | what it does | the earlier audit | re-run alone | the test that failed or the last one running | wall s |')
    print('|---|---|---|---|---|---|')
    for line in open(os.path.join(HERE, 'results_1b9c_hang.tsv')):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 7:
            print('| `%s` | %s | %s | **%s** | %s | %s |' % (p[0], short(p[1], 80), '', p[2], short(p[3] or p[5], 110), p[6]))
elif what == 'closed':
    print('| list | mutant | what it does | the test (SaGap1b9cTests / SesaGap1b9cTests) |')
    print('|---|---|---|---|')
    for (k, i), t in sorted(CLOSED.items()):
        m = pool[(k, i)]
        print('| %s | `%s` | %s | `%s` |' % ('1b.7d' if k == 'd' else '1b.7e', i, short(m['desc'], 90), t.split('.')[-1]))
elif what == 'equivalent':
    print('| list | mutant | what it does | argument |')
    print('|---|---|---|---|')
    for k in 'de':
        for i, a in sorted(ns[k]['EQUIVALENT_1B9C'].items()):
            print('| %s | `%s` | %s | %s |' % ('1b.7d' if k == 'd' else '1b.7e', i, short(pool[(k, i)]['desc'], 80), short(a, 400)))
elif what == 'fresh':
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, 'results_1b9c_f_*.tsv'))):
        for line in open(f):
            p = line.rstrip('\n').split('\t')
            rows.append(p)
    for p in rows:
        print('| `%s` | %s | %s | %s |' % (p[0], short(p[1], 80), short(' '.join(x for x in p[2:5] if x), 110), p[5]))
else:
    arith = json.load(open(os.path.join(HERE, '1b9c_arith.json')))
    print(json.dumps(arith, indent=1)[:200])
