#!/usr/bin/env python3
"""Tables of gate G3 of work package 1b.7b from rows_1b7b_g3.tsv (written by run_1b7b_g3.sh): the reader (production against prototype, medians of the per-run medians), the writer (production sparse and fallocate against
   the prototype, device writes per output byte and per symbol, time, RssAnon) and the 64 KB check of the harness.   usage: tables_1b7b_g3.py <rows file>"""
import collections, statistics, sys

rows = [l.rstrip('\n').split('\t') for l in open(sys.argv[1]) if l.strip()]
reader = collections.defaultdict(list)   # (input, part, what, variant) -> medians of the runs
for r in rows:
    if r[0].startswith('#') or len(r) < 9 or r[2] not in ('proto', 'prod', 'proto-vector', 'prod-vector', 'proto-scalar', 'prod-scalar'):
        continue
    label = r[0]   # <input>.<part>.run<k>
    inp, part = label.split('.')[0], label.split('.')[1]
    reader[(inp, part, r[1], r[2])].append(float(r[4]))
print('READER (ns per operation; median over runs of the median of 5 passes; ratio = production / prototype)')
done = set()
for (inp, part, what, variant), v in sorted(reader.items()):
    if variant.startswith('prod'):
        continue
    other = variant.replace('proto', 'prod')
    a, b = statistics.median(v), statistics.median(reader[(inp, part, what, other)])
    print('%-8s %-9s %-10s %-13s proto %9.3f  prod %9.3f  ratio %+6.1f %%   runs %s / %s' % (inp, part, what, variant.replace('proto', '').strip('-') or 'scalar', a, b, (b / a - 1) * 100,
          ','.join('%.2f' % x for x in v), ','.join('%.2f' % x for x in reader[(inp, part, what, other)])))
w = collections.defaultdict(dict)
for r in rows:
    if r[0].startswith('w') and len(r) == 10 and not r[0].startswith('#'):
        w[(r[0].split('.')[0], r[1], r[2])].setdefault(r[0], {}).update({'created': float(r[4]), 'seconds': float(r[5]), 'rss_before': int(r[6]), 'rss_max': int(r[7]), 'bytes': int(r[8]), 'windows': int(r[9])})
    if r[0] == '#dev':
        key = (r[1].split('.')[0], r[2], r[3])
        w[key].setdefault(r[1], {}).update({'dev_w': int(r[6]), 'dev_r': int(r[8]), 'wall': float(r[10]), 'wall_sync': float(r[12])})
print()
print('WRITER (400 MB output for 10^8 symbols; medians over the runs; device = block layer after sync, whole process; time = creation of the file to the end of the writes)')
print('%-6s %-10s %-5s %5s %9s %9s %10s %10s %9s %10s' % ('label', 'variant', 'pat', 'runs', 'create s', 'total s', 'dev W B/ch', 'dev R B/ch', 'wall+sync', 'RssAnon MB'))
for (lab, var, pat), runs in sorted(w.items()):
    def med(f):
        v = [x[f] for x in runs.values() if f in x]
        return statistics.median(v) if v else float('nan')
    n = 1e8
    print('%-6s %-10s %-5s %5d %9.4f %9.3f %10.3f %10.3f %9.3f %10.1f' % (lab, var, pat, len(runs), med('created'), med('seconds'), med('dev_w') / n, med('dev_r') / n, med('wall_sync'), (med('rss_max') - med('rss_before')) / 1000.0))
