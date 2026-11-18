#!/usr/bin/env python3
"""Reads rows_1b7c_g3.tsv (written by run_1b7c_g3.sh) and prints the G3 table of work package 1b.7c: median of the build call per cell and variant (proto = the prototype's fused PH4 + PH5, prod = LcpChildBuilder, prodnr = the
   builder without its RSS sampler), the ratio production / prototype, the RssAnon of the test host (own sampler and the independent one), the counted I/O, the spill and window figures, the cost of the builder's RSS sampler
   and the identity checks. usage: tables_1b7c_g3.py <rows file>"""
import statistics, sys
from collections import defaultdict

rows, cmp_rows, ext = [], [], {}
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
    cell = r[0].rsplit('.', 1)[0]
    cells[cell][r[1]].append(r)


def med(xs):
    return statistics.median(xs)


print('%-26s %-7s %-9s %s' % ('cell', 'variant', 'median s', 'runs (s)'))
for cell in sorted(cells):
    for v in sorted(cells[cell]):
        xs = [float(r[2]) for r in cells[cell][v]]
        print('%-26s %-7s %-9.3f %s' % (cell, v, med(xs), ' '.join('%.3f' % x for x in xs)))
print()
print('ratio production / prototype (median of the build call)')
for cell in sorted(cells):
    if 'proto' in cells[cell] and 'prod' in cells[cell]:
        a = med([float(r[2]) for r in cells[cell]['proto']])
        b = med([float(r[2]) for r in cells[cell]['prod']])
        print('  %-26s %.3f  (%+.1f %%)  proto %.3f s, prod %.3f s' % (cell, b / a, (b / a - 1) * 100, a, b))
print()
print('cost of the builder RSS sampler (prod / prodnr, warm, parameters A)')
for cell in sorted(cells):
    if cell.endswith('Aoverhead.warm') and 'prod' in cells[cell] and 'prodnr' in cells[cell]:
        a = med([float(r[2]) for r in cells[cell]['prodnr']])
        b = med([float(r[2]) for r in cells[cell]['prod']])
        print('  %-26s %.3f s / %.3f s = %.3f  (%+.1f %%)' % (cell, b, a, b / a, (b / a - 1) * 100))
print()
print('RssAnon of the test host in kB: own sampler (max over the runs of a cell) and the independent one; before = before the build')
for cell in sorted(cells):
    for v in sorted(cells[cell]):
        mx = max(int(r[4]) for r in cells[cell][v])
        before = med([int(r[3]) for r in cells[cell][v]])
        extmax = max([x for (lab, var), x in ext.items() if lab.rsplit('.', 1)[0] == cell and var == v] or [0])
        print('  %-26s %-7s own max %7d  independent max %7d  before %7d' % (cell, v, mx, extmax, before))
print()
print('counted I/O and the rest (first run of the cell; all runs of a cell are equal for these columns unless marked)')
names = ['arena_peak', 'child_arena_peak', 'read', 'write', 'mapped_lcp', 'child_write', 'spill_written', 'spill_read', 'child_flushed', 'child_reopened', 'stack_peak', 'workers', 'R', 'scratch_peak', 'b1_s', 'b2_s', 'b3_s', 'compare_steps']
for cell in sorted(cells):
    for v in sorted(cells[cell]):
        first = cells[cell][v][0]
        same = all(r[6:20] == first[6:20] for r in cells[cell][v])
        print('  %-26s %-7s %s%s' % (cell, v, ' '.join('%s=%s' % (n, x) for n, x in zip(names, first[6:])), '' if same else '  (runs differ)'))
print()
print('equality of the counted figures, production against prototype (read, write, mapped_lcp, child_write, spill_written, spill_read, child_flushed, child_reopened, stack_peak, workers, R)')
for cell in sorted(cells):
    if 'proto' in cells[cell] and 'prod' in cells[cell]:
        a, b = cells[cell]['proto'][0], cells[cell]['prod'][0]
        idx = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        print('  %-26s %s' % (cell, 'EQUAL' if [a[i] for i in idx] == [b[i] for i in idx] else 'DIFFERENT: proto %s prod %s' % ([a[i] for i in idx], [b[i] for i in idx])))
print()
bad = [c for c in cmp_rows if c[3] != 'equal-to-reference']
print('outputs cmp\'ed with the reference files: %d runs, %d different' % (len(cmp_rows), len(bad)))
for c in bad:
    print('  ', c)
