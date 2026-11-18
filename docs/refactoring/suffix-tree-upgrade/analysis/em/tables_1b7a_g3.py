#!/usr/bin/env python3
"""Reads rows_1b7a_g3.tsv (written by run_1b7a_g3.sh) and prints the G3 table of work package 1b.7a: median of the Build call per cell and variant, the ratio production / prototype, RssAnon, arena peak, I/O, the
   counter overhead and the identity checks. usage: tables_1b7a_g3.py <rows file>"""
import statistics, sys
from collections import defaultdict

rows = []
cmp_rows = []
ext = {}
for line in open(sys.argv[1]):
    line = line.rstrip('\n')
    if line.startswith('#cmp'):
        cmp_rows.append(line.split('\t'))
    elif line.startswith('#ext'):
        p = line.split('\t')
        ext[(p[1], p[2])] = int(p[3])
    elif line.startswith('#'):
        continue
    else:
        rows.append(line.split('\t'))

cells = defaultdict(lambda: defaultdict(list))
for r in rows:
    label, variant = r[0], r[1]
    cell = label.rsplit('.', 1)[0]
    cells[cell][variant].append(r)

def med(xs):
    return statistics.median(xs)

print('%-16s %-8s %-7s %s' % ('cell', 'variant', 'median s', 'runs (s)'))
for cell in sorted(cells):
    for v in sorted(cells[cell]):
        xs = [float(r[4]) for r in cells[cell][v]]
        print('%-16s %-8s %-7.3f %s' % (cell, v, med(xs), ' '.join('%.3f' % x for x in xs)))
print()
print('ratio production / prototype (median of the Build call)')
for cell in sorted(cells):
    if 'proto' in cells[cell] and 'prod' in cells[cell]:
        a = med([float(r[4]) for r in cells[cell]['proto']])
        b = med([float(r[4]) for r in cells[cell]['prod']])
        print('  %-14s %.3f  (%+.1f %%)' % (cell, b / a, (b / a - 1) * 100))
print()
print('counter overhead (production with counters / without), warm, parameters A')
c = cells.get('Aoverhead.warm')
if c and 'prod' in c and 'prodnc' in c:
    a = med([float(r[4]) for r in c['prodnc']])
    b = med([float(r[4]) for r in c['prod']])
    print('  %.3f s / %.3f s = %.3f  (%+.1f %%)' % (b, a, b / a, (b / a - 1) * 100))
print()
print('RssAnon (kB): in-process 10 ms production sampler, maximum over runs; before = at the start of the test, external = python sampler of /proc/<testhost>/status every 20 ms')
for cell in sorted(cells):
    for v in sorted(cells[cell]):
        mx = max(int(r[6]) for r in cells[cell][v])
        bf = max(int(r[5]) for r in cells[cell][v])
        ap = max(int(r[8]) for r in cells[cell][v])
        e = max([x for (lab, vv), x in ext.items() if vv == v and lab.rsplit('.', 1)[0] == cell] or [0])
        print('  %-16s %-8s max %7d  before %6d  delta %7d  arena peak %10d B  external max %7d' % (cell, v, mx, bf, mx - bf, ap, e))
print()
print('I/O (bytes) of the first run of each cell: read / write / mapped / scratch peak')
for cell in sorted(cells):
    for v in sorted(cells[cell]):
        r = cells[cell][v][0]
        print('  %-16s %-8s %s' % (cell, v, ' / '.join(r[9:13])))
print()
eq = sum(1 for p in cmp_rows if p[3] == 'equal-to-reference')
print('cmp of every output with the in-RAM ISA of the E9 setup: %d equal, %d different' % (eq, len(cmp_rows) - eq))
