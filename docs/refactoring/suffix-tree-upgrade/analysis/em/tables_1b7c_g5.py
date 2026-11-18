#!/usr/bin/env python3
"""Reads rows_1b7c_g5.tsv (written by run_1b7c_g5.sh) and prints the G5 table of work package 1b.7c: per phase of a process (cold = the first pass, warm = the later passes) and per measured part the median over the
   processes (cold) or over all warm passes of all processes (warm) for the Packed and the Tagged codec and the ratio Packed / Tagged. usage: tables_1b7c_g5.py <rows file>"""
import statistics, sys
from collections import defaultdict

cols = ['walk_s', 'count_present_s', 'count_absent_s', 'locate12_s']
data = defaultdict(list)
steps = 0
for line in open(sys.argv[1]):
    if line.startswith('#'):
        continue
    p = line.rstrip('\n').split('\t')
    codec, phase = p[1], p[2]
    steps = int(p[5])
    vals = [float(p[4]), float(p[6]), float(p[7]), float(p[8])]
    for c, v in zip(cols, vals):
        data[(phase, c, codec)].append(v)
print('child steps per walk: %d' % steps)
print('%-6s %-16s %9s %9s %8s   runs (packed | tagged)' % ('phase', 'part', 'packed s', 'tagged s', 'ratio'))
for phase in ['cold', 'warm']:
    for c in cols:
        a, b = data[(phase, c, 'packed')], data[(phase, c, 'tagged')]
        ma, mb = statistics.median(a), statistics.median(b)
        print('%-6s %-16s %9.4f %9.4f %8.3f   n=%d|%d  min-max packed %.3f-%.3f tagged %.3f-%.3f' % (phase, c, ma, mb, ma / mb, len(a), len(b), min(a), max(a), min(b), max(b)))
walk_a, walk_b = statistics.median(data[('warm', 'walk_s', 'packed')]), statistics.median(data[('warm', 'walk_s', 'tagged')])
print('\nper child step, warm: packed %.1f ns, tagged %.1f ns (+%.1f %%)' % (walk_a / steps * 1e9, walk_b / steps * 1e9, (walk_a / walk_b - 1) * 100))
