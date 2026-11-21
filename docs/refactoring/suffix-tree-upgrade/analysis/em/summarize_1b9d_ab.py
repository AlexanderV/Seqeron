#!/usr/bin/env python3
"""Table of the in-process A/B rows of run_1b9d_ab.sh (`sesa ab`): per input and cache state, the seconds of A (the growth is sparse: the behaviour before 1b.9d) and B (reserved by writing zeros) as medians [min-max] over
the pairs, the median of the ratios B / A of the pairs [min-max] for the whole build and for PH4 + PH5, the seconds that were spent inside the reservation (median) and the bytes it wrote, and the share of the build.
The warm-up pair (pair 0) is left out. usage: summarize_1b9d_ab.py rows..."""
import collections, statistics, sys

data = collections.defaultdict(lambda: collections.defaultdict(list))  # (input.cache) -> key -> values
for path in sys.argv[1:]:
    for line in open(path):
        p = line.rstrip('\n').split('\t')
        if len(p) < 6 or p[0] != 'R' or p[1] != 'sesa':
            continue
        metric, value, notes = p[2], p[3], p[5]
        label = notes.split()[0]
        parts = label.split('.')
        if metric in ('total_s', 'ph45_s', 'counted_reserved_write', 'grown_bytes', 'file_size') and len(parts) == 4 and parts[2] != '0':
            data['.'.join(parts[:2])][metric + '.' + parts[3]].append(float(value))
        elif metric == 'reserve_s' and len(parts) == 4 and parts[2] != '0':
            data['.'.join(parts[:2])]['reserve_s'].append(float(value))
        elif metric in ('pair_ratio_total', 'pair_ratio_ph45'):
            name = notes.split()[0]
            parts = name.split('.')
            if parts[-1] == 'pair':
                pass
            data[notes.split()[0]]  # noqa

# pair ratios are re-derived from the per-run totals (the order of the rows is the order of the runs)
def med(v, fmt='%.3f'):
    if not v:
        return '-'
    return (fmt + ' [' + fmt + '-' + fmt + ']') % (statistics.median(v), min(v), max(v))

print('| input.cache | file MB | grown MB | A total s | B total s | ratio B / A of the pairs | A PH4+PH5 s | B PH4+PH5 s | ratio PH4+PH5 | reservation s (inside the call) | share of the build |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for key in sorted(data):
    d = data[key]
    a, b = d.get('total_s.A', []), d.get('total_s.B', [])
    if not a or len(a) != len(b):
        continue
    ratios = [y / x for x, y in zip(a, b)]
    a45, b45 = d.get('ph45_s.A', []), d.get('ph45_s.B', [])
    r45 = [y / x for x, y in zip(a45, b45)]
    res = d.get('reserve_s', [])
    grown = d.get('grown_bytes.B', [0])
    share = statistics.median(res) / statistics.median(b) * 100 if res else 0
    print('| %s | %.0f | %.0f | %s | %s | %s | %s | %s | %s | %s | %.2f %% |' % (key, statistics.median(d.get('file_size.B', [0])) / 1e6, statistics.median(grown) / 1e6, med(a, '%.2f'), med(b, '%.2f'), med(ratios, '%.3f'), med(a45, '%.2f'), med(b45, '%.2f'), med(r45, '%.3f'), med(res, '%.4f'), share))
