#!/usr/bin/env python3
"""The navigator rows of S2.5 (run_s25_nav_ab.sh): one line per cell file. The figure of a process is the minimum of its warm calls (and, second, the median of its last five calls); the pooled ratio is the median over the processes of the arm
   NEW divided by the median over the processes of the arm BASE (a time ratio: the speed is its inverse, the gate is a speed of at least 0.97, i.e. a time ratio of at most 1.031); the paired ratio is the median (with quartiles) of NEW / BASE of the two
   processes of one round. The null is the second copy of BASE (arm base2) treated the same way; its paired range (q25 - q75 and min - max) is the spread of the instrument for the cell. Checks: every row of every arm has the same count of results and checksum.
   usage: summarize_s25_nav.py <rows file>... > summary.tsv"""
import re, sys, statistics as st


def q(xs, f):
    xs = sorted(xs)
    return xs[int(f * (len(xs) - 1))]


print('cell\tn_proc\tpooled_new/base\tspeed\tpooled_null/base\tpaired_new/base\tq25\tq75\tpaired_null\tnull_q25\tnull_q75\tnull_min\tnull_max\tsteady_pooled\tsteady_null\tresults_checksum_equal\tfinished')
for path in sys.argv[1:]:
    mins, steady, res = {}, {}, set()
    finished = False
    for line in open(path):
        if line.startswith('done'):
            finished = True
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 3:
            continue
        t = '\t'.join(p[2:])
        m = re.search(r'warm_median_s\t[0-9.]+\ts\tmedian of \d+ warm calls \([^)]*\) \[([0-9.]+) \.\. ([0-9.]+)\]', t)
        m2 = re.search(r'steady_median_s\t([0-9.]+)\ts', t)
        r = re.search(r'_results\t(\d+)\tcount', t)
        c = re.search(r'_checksum\t(\d+)\tcount', t)
        if m:
            mins.setdefault(p[1], {})[int(p[0])] = float(m.group(1))
        if m2:
            steady.setdefault(p[1], {})[int(p[0])] = float(m2.group(1))
        if r and c:
            res.add((r.group(1), c.group(1)))
    b, b2, n = mins['base'], mins['base2'], mins['new']
    pooled = st.median(n.values()) / st.median(b.values())
    pooled0 = st.median(b2.values()) / st.median(b.values())
    rr = [n[i] / b[i] for i in b if i in n]
    r0 = [b2[i] / b[i] for i in b if i in b2]
    sp = st.median(steady['new'].values()) / st.median(steady['base'].values())
    sp0 = st.median(steady['base2'].values()) / st.median(steady['base'].values())
    cell = path.split('/')[-1].replace('rows_', '').replace('.tsv', '')
    print('%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s' % (
        cell, len(rr), pooled, 1 / pooled, pooled0, st.median(rr), q(rr, .25), q(rr, .75), st.median(r0), q(r0, .25), q(r0, .75), min(r0), max(r0), sp, sp0, 'yes' if len(res) == 1 else 'NO %s' % sorted(res), finished))
