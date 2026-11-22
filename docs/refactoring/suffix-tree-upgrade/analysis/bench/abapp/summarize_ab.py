#!/usr/bin/env python3
"""Summary of the rows of run_ab.sh: per file and operation the median of the per-process ratios B / A, their range and the number of processes.
usage: summarize_ab.py <rows file>..."""
import re, sys, statistics
pat = re.compile(r'^AB (\S+)\s+A\s+(\d+)\s+B\s+(\d+)\s+B/A median-of-paired ([\d.]+)')
for path in sys.argv[1:]:
    ops = {}
    for line in open(path):
        m = pat.match(line)
        if m: ops.setdefault(m.group(1), []).append((float(m.group(4)), int(m.group(2)), int(m.group(3))))
    name = path.rsplit('/', 1)[-1]
    cells = []
    for op in ('count_present', 'locate', 'count_absent'):
        r = ops.get(op)
        if not r: cells.append(f'{op} -'); continue
        v = [x[0] for x in r]
        cells.append(f'{op} {statistics.median(v):.3f} [{min(v):.3f}-{max(v):.3f}] n={len(v)}')
    print(f'{name}: ' + ' | '.join(cells))
