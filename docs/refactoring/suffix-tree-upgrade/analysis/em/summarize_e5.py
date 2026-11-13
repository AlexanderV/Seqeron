#!/usr/bin/env python3
"""Medians of the raw rows of E5 (and of the PH4 reference runs): python3 -I summarize_e5.py <metrics comma separated> <case regex> <file.tsv...>
Rows are R <TAB> experiment <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes (the rows without a run column are skipped).
Prints one line per case: median [min-max] (runs) of every metric asked for, text values as they are (the first)."""
import re, statistics, sys
from collections import defaultdict

md = sys.argv[1] == '--md'
args = sys.argv[2:] if md else sys.argv[1:]
metrics = args[0].split(',')
pattern = re.compile(args[1])
data = defaultdict(lambda: defaultdict(list))
order = []
for path in args[2:]:
    for line in open(path):
        c = line.rstrip('\n').split('\t')
        if len(c) < 7 or c[0] != 'R' or not c[3].isdigit():
            continue
        case, run, metric, value = c[2], c[3], c[4], c[5]
        if not pattern.search(case) or metric not in metrics:
            continue
        if case not in order:
            order.append(case)
        data[case][metric].append(value)
if md:
    print('| case | ' + ' | '.join(metrics) + ' |')
    print('|---|' + '---|' * len(metrics))
for case in order:
    out = []
    for m in metrics:
        v = data[case].get(m, [])
        try:
            f = [float(x) for x in v]
        except ValueError:
            out.append(v[0] if v else '-')
            continue
        if not f:
            out.append('-')
        elif len(f) == 1:
            out.append('%.4g' % f[0])
        else:
            out.append('%.4g [%.4g-%.4g]' % (statistics.median(f), min(f), max(f)))
    if md:
        print('| %s | %s |' % (case, ' | '.join(out)))
    else:
        print('%-28s %s' % (case, '  '.join('%s=%s' % (m, o) for m, o in zip(metrics, out))))
