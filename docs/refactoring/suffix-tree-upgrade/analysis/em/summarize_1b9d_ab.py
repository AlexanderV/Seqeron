#!/usr/bin/env python3
"""Table of the in-process A/B rows of run_1b9d_ab.sh (`sesa ab`): per input and cache state, the seconds of A (the growth is sparse: the behaviour before 1b.9d) and B (reserved by writing zeros) as medians [min-max] over
the pairs, the median and the geometric mean of the ratios B / A of the pairs [min-max] for the whole build and for PH4 + PH5, the seconds that were spent inside the reservation (median) and their share of the build.
The warm-up pair (pair 0) is left out. A cell whose growth is zero (dna) is the null cell: A and B run the same code there, so its ratios are the noise of the instrument.
usage: summarize_1b9d_ab.py rows..."""
import collections, math, statistics, sys

data = collections.defaultdict(lambda: collections.defaultdict(list))  # (input.cache) -> key -> values, in run order
for path in sys.argv[1:]:
    for line in open(path):
        p = line.rstrip('\n').split('\t')
        if len(p) < 6 or p[0] != 'R' or p[1] != 'sesa':
            continue
        metric, value, notes = p[2], p[3], p[5]
        parts = notes.split()[0].split('.')          # input.cache.pair.variant
        if len(parts) != 4 or parts[2] == '0':
            continue
        key = '.'.join(parts[:2])
        if metric in ('total_s', 'ph45_s', 'grown_bytes', 'file_size'):
            data[key][metric + '.' + parts[3]].append(float(value))
        elif metric == 'reserve_s':
            data[key]['reserve_s'].append(float(value))

def med(v, fmt='%.3f'):
    return '-' if not v else (fmt + ' [' + fmt + '-' + fmt + ']') % (statistics.median(v), min(v), max(v))

def geo(v):
    return math.exp(sum(math.log(x) for x in v) / len(v))

print('| input.cache | file MB | grown MB | A total s | B total s | ratio B / A: median [min-max] | geometric mean | A PH4+PH5 s | B PH4+PH5 s | ratio PH4+PH5 | reservation s (inside the call) | share of the build |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|')
for key in sorted(data):
    d = data[key]
    a, b = d.get('total_s.A', []), d.get('total_s.B', [])
    if not a or len(a) != len(b):
        continue
    ratios = [y / x for x, y in zip(a, b)]
    a45, b45 = d.get('ph45_s.A', []), d.get('ph45_s.B', [])
    r45 = [y / x for x, y in zip(a45, b45)]
    res = d.get('reserve_s', [])
    share = statistics.median(res) / statistics.median(b) * 100 if res else 0
    print('| %s (%d pairs) | %.0f | %.0f | %s | %s | %s | %.3f | %s | %s | %s | %s | %.2f %% |' % (key, len(ratios), statistics.median(d.get('file_size.B', [0])) / 1e6, statistics.median(d.get('grown_bytes.B', [0])) / 1e6, med(a, '%.2f'), med(b, '%.2f'), med(ratios), geo(ratios), med(a45, '%.2f'), med(b45, '%.2f'), med(r45), med(res, '%.4f'), share))
