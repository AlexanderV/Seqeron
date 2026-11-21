#!/usr/bin/env python3
"""Tables of rows_1b9d_tail.tsv (run_1b9d_tail.sh): per size, cache state and variant, the median [min-max] of the seconds of the growth, of the write, of the sync and of the total, and of the bytes that the process wrote
(wchar), the pages it dirtied and the bytes that reached the device until the end of the sync; the cost of 'zero' and 'worst' against 'sparse' (the difference of the medians of the totals)."""
import collections, statistics, sys

rows = collections.defaultdict(lambda: collections.defaultdict(list))  # (size, cache, variant) -> metric -> values
for line in open(sys.argv[1] if len(sys.argv) > 1 else 'rows_1b9d_tail.tsv'):
    p = line.rstrip('\n').split('\t')
    if len(p) < 6 or p[0] != 'R':
        continue
    label = p[5].split()          # r1 sparse 100 MiB cold ...
    size, variant, cache = int(label[2]), label[1], 'cold' if 'cold' in label else 'warm'
    rows[(size, cache, variant)][p[2]].append(float(p[3]))

def med(v):
    return '%.3f [%.3f-%.3f]' % (statistics.median(v), min(v), max(v))

def mb(v):
    return '%.0f' % (statistics.median(v) / 1e6)

print('| size MiB | cache | variant | growth s | write s | sync s | total s | wchar MB | pages dirtied MB | device MB (to the end of the sync) |')
print('|---|---|---|---|---|---|---|---|---|---|')
for size in sorted({k[0] for k in rows}):
    for cache in ('cold', 'warm'):
        for variant in ('sparse', 'zero', 'worst'):
            m = rows.get((size, cache, variant))
            if not m:
                continue
            print('| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (size, cache, variant, med(m['grow_s']), med(m['write_s']), med(m['sync_s']), med(m['total_s']), mb(m['proc_wchar']), mb(m['proc_write_bytes']), mb(m['dev_written_bytes'])))
print()
print('| size MiB | cache | zero - sparse, total s (median of the three differences by round) | worst - sparse | zero / sparse | device MB zero - sparse |')
print('|---|---|---|---|---|---|')
for size in sorted({k[0] for k in rows}):
    for cache in ('cold', 'warm'):
        s, z, w = (rows.get((size, cache, v)) for v in ('sparse', 'zero', 'worst'))
        if not (s and z and w):
            continue
        dz = [a - b for a, b in zip(z['total_s'], s['total_s'])]
        dw = [a - b for a, b in zip(w['total_s'], s['total_s'])]
        print('| %d | %s | %+.3f [%+.3f..%+.3f] | %+.3f [%+.3f..%+.3f] | %.2f | %+.0f |' % (size, cache, statistics.median(dz), min(dz), max(dz), statistics.median(dw), min(dw), max(dw), statistics.median(z['total_s']) / statistics.median(s['total_s']), (statistics.median(z['dev_written_bytes']) - statistics.median(s['dev_written_bytes'])) / 1e6))
