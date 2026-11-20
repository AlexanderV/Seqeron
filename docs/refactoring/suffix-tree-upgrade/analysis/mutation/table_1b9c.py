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
    lst = json.load(open(os.path.join(HERE, '1b9c_h_all_listed.json')))
    earlier = {m['package'] + ':' + m['id'] for m in lst[:17]}
    checked = {}
    cf = os.path.join(HERE, 'results_1b9c_hang_checked.tsv')
    if os.path.exists(cf):
        for line in open(cf):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 7:
                checked[p[0]] = p
    def label(p):
        """crash = the host died within seconds; hang = it was killed by the time limit of a test (the runtime prints "Aborted." when vstest kills it after the 60 s of --blame-hang-timeout)"""
        kind = p[2]
        secs = int(p[6]) if p[6].isdigit() else 0
        return 'hang' if kind == 'crash' and secs >= 40 and 'Aborted' in p[4] else kind

    print('| mutant | what it does | counted as killed by | re-run alone, production build | the test that failed, or the last one running | wall s | checked build (`EsaBounds`) |')
    print('|---|---|---|---|---|---|---|')
    for line in open(os.path.join(HERE, 'results_1b9c_hang.tsv')):
        p = line.rstrip('\n').split('\t')
        if len(p) >= 7:
            c = checked.get(p[0])
            ctext = '' if c is None else '**%s**: %s' % (label(c), short(c[3] or c[5].split('.')[-1], 60))
            why = p[4].replace('Fatal error. System.', '').replace('Fatal error.', '')
            print('| `%s` | %s | %s | **%s** %s | %s | %s | %s |' % (p[0], short(p[1], 80), 'a watchdog or a crash (1b.7%s)' % p[0][0] if p[0] in earlier else 'a watchdog or a crash (this package)', label(p), short(why, 50), short((p[3] or p[5]).split('.')[-1], 90), p[6], ctext))
elif what == 'closed':
    print('| list | mutant | what it does | the test (SaGap1b9cTests / SesaGap1b9cTests) |')
    print('|---|---|---|---|')
    for (k, i), t in sorted(CLOSED.items()):
        m = pool[(k, i)]
        print('| %s | `%s` | %s | `%s` |' % ('1b.7d' if k == 'd' else '1b.7e', i, short(m['desc'], 90), t.split('.')[-1]))
elif what == 'closed-by-test':
    by_test = collections.OrderedDict()
    for (k, i), t in sorted(CLOSED.items(), key=lambda kv: (kv[1], kv[0])):
        by_test.setdefault(t, []).append(('' if k == 'd' else 'e:') + i)
    print('| test | mutants it kills (1b.7d ids; `e:` = 1b.7e) |')
    print('|---|---|')
    for t, ids in by_test.items():
        print('| `%s` | %s |' % (t, ', '.join('`%s`' % x for x in ids)))
elif what == 'equivalent':
    print('| list | mutant | what it does | argument |')
    print('|---|---|---|---|')
    for k in 'de':
        for i, a in sorted(ns[k]['EQUIVALENT_1B9C'].items()):
            print('| %s | `%s` | %s | %s |' % ('1b.7d' if k == 'd' else '1b.7e', i, short(pool[(k, i)]['desc'], 80), short(a, 400)))
elif what == 'fresh':
    exec(open(os.path.join(HERE, '1b9c_fresh_fate.py')).read())
    lst = {}
    for f in sorted(glob.glob(os.path.join(HERE, '1b9c_f_*_listed.json'))):
        for m in json.load(open(f)):
            lst[m['id']] = m
    p1, p2 = {}, {}
    for f in sorted(glob.glob(os.path.join(HERE, 'results_1b9c_f_*.tsv'))):
        for line in open(f):
            p = line.rstrip('\n').split('\t')
            p1[p[0]] = p[2] if len(p) > 2 else ''
    for f in sorted(glob.glob(os.path.join(HERE, 'results_1b9c_f2_*.tsv'))):
        for line in open(f):
            p = line.rstrip('\n').split('\t')
            p2[p[0]] = p[2] if len(p) > 2 else ''
    kills = [i for i in sorted(lst) if i not in FRESH_FATE and p1.get(i, '').startswith('killed')]
    print('| mutant | what it does | first run: the new tests alone | closed by | the same run after the closing test (final tests, new classes alone) |')
    print('|---|---|---|---|---|')
    for i in sorted(FRESH_FATE):
        now = p1.get(i, '')
        print('| `%s` | %s | survived%s | %s | %s |' % (i, short(lst[i]['desc'], 80), ', then ' + p2[i].split('(')[0] + ' in the whole default set' if p2.get(i) and not p2[i].startswith('killed') else ('' if not p2.get(i) else ', then killed by the default set'), FRESH_FATE[i], short(now.replace('[only]', ''), 70) if now.startswith('killed') else now))
    print()
    print('killed in pass 1 (%d): %s' % (len(kills), ', '.join('`%s`' % i for i in kills)))
else:
    # the arithmetic per file: summary_1b9c_files.tsv (summarize_1b9c.py)
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(HERE, 'summary_1b9c_files.tsv'))]
    head = rows[0]
    print('| file | listed | run by the earlier audit | not run | killed by an assertion (earlier audit) | killed (stage 1) | killed (stage 2) | watchdog / crash kill: failing assertion on its own | watchdog / crash kill: hang or crash | closed by a test | equivalent (argued) | open |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|')
    sub = {'1b7d': [0] * 11, '1b7e': [0] * 11}
    for r in rows[1:]:
        if r[0] == 'TOTAL':
            continue
        for i in range(11):
            sub[r[0].split('/')[0]][i] += int(r[1 + i])
        print('| %s | %s |' % (r[0], ' | '.join(r[1:])))
    for k in ('1b7d', '1b7e'):
        print('| **%s total** | %s |' % (k, ' | '.join('**%d**' % v for v in sub[k])))
    t = rows[-1]
    print('| **both** | %s |' % ' | '.join('**%s**' % v for v in t[1:]))
