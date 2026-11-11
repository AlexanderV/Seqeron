#!/usr/bin/env python3
"""Markdown tables of RESULTS-E3B.md from the raw rows (R <TAB> E3b <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes).
usage: tables_e3b.py <rows dir>   (e3b_g100.tsv e3b_small.tsv e3b_slope.tsv e3b_slopetouch.tsv e3b_fam.tsv e3b_extra.tsv)"""
import sys, os, math, statistics, collections, re

D = sys.argv[1]

def load(name):
    data = collections.OrderedDict()
    for line in open(os.path.join(D, name)):
        c = line.rstrip('\n').split('\t')
        if len(c) < 7 or c[0] != 'R' or c[1] != 'E3b':
            continue
        # rows of `extra` have no run column when --run was not given: detect by the metric position
        if len(c) >= 8:
            case, run, metric, value, unit, notes = c[2], c[3], c[4], c[5], c[6], c[7]
        else:
            case, run, metric, value, unit, notes = c[2], '', c[3], c[4], c[5], c[6]
        data.setdefault(case, collections.OrderedDict()).setdefault(metric, []).append((value, notes))
    return data

def merge(*ds):
    out = collections.OrderedDict()
    for d in ds:
        for case, m in d.items():
            for metric, vals in m.items():
                out.setdefault(case, collections.OrderedDict()).setdefault(metric, []).extend(vals)
    return out

def nums(d, case, metric):
    return [float(v) for v, _ in d[case][metric]]

def med(d, case, metric, fmt='{:.2f}', scale=1.0):
    v = [x * scale for x in nums(d, case, metric)]
    if len(v) == 1:
        return fmt.format(v[0])
    return (fmt + ' [' + fmt + '-' + fmt + ']').format(statistics.median(v), min(v), max(v))

def first(d, case, metric):
    return d[case][metric][0][0]

def table_inputs(d, cases, labels, nsym):
    print('| input | n | runs | PH1 s (4 threads) | ns/symbol | RssAnon MB | levels | shared names at level 0 | explicit I/O B/char (read+write) | mapped write B/char | block-layer read B/char (cold) | peak scratch B/char |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|')
    for c, lab in zip(cases, labels):
        n = float(first(d, c, 'n'))
        runs = len(d[c]['ph1_s'])
        er = nums(d, c, 'io_explicit_read')[0] / n
        ew = nums(d, c, 'io_explicit_write')[0] / n
        mw = nums(d, c, 'io_mapped_write')[0] / n
        br = statistics.median(nums(d, c, 'blocklayer_read')) / n
        sp = statistics.median(nums(d, c, 'scratch_peak')) / n
        print(f"| {lab} | {int(n):,} | {runs} | {med(d, c, 'ph1_s')} | {med(d, c, 'ns_per_symbol', '{:.1f}')} | {med(d, c, 'rssanon_peak', '{:.1f}')} | {int(nums(d, c, 'levels')[0])} | {nums(d, c, 'nonunique_share')[0]:.3f} % | {er + ew:.3f} ({er:.3f}+{ew:.3f}) | {mw:.3f} | {br:.3f} | {sp:.3f} |")

def slope(d, cases, label):
    xs, ys = [], []
    for c in cases:
        n = float(first(d, c, 'n'))
        xs.append(math.log2(n))
        ys.append(statistics.median(nums(d, c, 'rssanon_peak')))
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    rows = []
    for c, y in zip(cases, ys):
        rows.append(f"{int(float(first(d, c, 'n'))/1e6)} M: {y:.1f} MB [{min(nums(d, c, 'rssanon_peak')):.1f}-{max(nums(d, c, 'rssanon_peak')):.1f}] ({len(nums(d, c, 'rssanon_peak'))} runs)")
    # the same fit through every single run
    px = [math.log2(float(first(d, c, 'n'))) for c in cases for _ in nums(d, c, 'rssanon_peak')]
    py = [v for c in cases for v in nums(d, c, 'rssanon_peak')]
    mpx, mpy = sum(px) / len(px), sum(py) / len(py)
    allfit = sum((x - mpx) * (y - mpy) for x, y in zip(px, py)) / sum((x - mpx) ** 2 for x in px)
    pair = [(ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]) for i in range(len(ys) - 1)]
    print(f"| {label} | " + ' ; '.join(rows) + f" | {b:.2f} (all runs: {allfit:.2f}) | " + ' / '.join(f'{p:.2f}' for p in pair) + ' |')

if __name__ == '__main__':
    g = load('e3b_g100.tsv')
    sm = load('e3b_small.tsv')
    print('## inputs at 100 Mbp and small inputs (4 threads, scaled, cold cache, plain buffers)')
    table_inputs({**g, **sm}, ['dna_100m.t4', 'genome_100m.t4', 'dna_10m.t4', 'rep10.t4', 'a10.t4', 'fib10.t4', 'nrun10.t4', 'genome10.t4', 'tandem.t4'],
                 ['dna_100m (random)', 'genome-like 100 M', 'dna_10m', 'rep_10m', 'a^n 10 M', 'Fibonacci 10 M', 'N run 1 M in 10 M', 'genome-like 10 M', '171-bp tandem x 10^4'], None)
    print()
    sl = merge(load('e3b_slope.tsv'), load('e3b_slope2.tsv'))
    st = load('e3b_slopetouch.tsv')
    print('| series | RssAnon peak per n (median [min-max], runs per n in the text) | least-squares slope MB per doubling | step slopes MB per doubling |')
    print('|---|---|---|---|')
    slope(sl, ['dna10.slope', 'dna50.slope', 'dna_100m.slope', 'dna200.slope'], 'random DNA, plain buffers, 10-200 M')
    slope(sl, ['dna10.slope', 'dna50.slope', 'dna_100m.slope'], 'random DNA, plain buffers, 10-100 M')
    slope(sl, ['genome10.slope', 'genome50.slope', 'genome_100m.slope'], 'genome-like, plain buffers, 10-100 M')
    slope(st, ['dna10.slopetouch', 'dna50.slopetouch', 'dna_100m.slopetouch', 'dna200.slopetouch'], 'random DNA, all buffers touched, 10-200 M')
    slope(st, ['dna10.slopetouch', 'dna50.slopetouch', 'dna_100m.slopetouch'], 'random DNA, all buffers touched, 10-100 M')
    slope(st, ['genome10.slopetouch', 'genome50.slopetouch', 'genome_100m.slopetouch'], 'genome-like, all buffers touched, 10-100 M')
    print()
    ex = load('e3b_extra.tsv')
    for c in ex:
        print(c, {m: ex[c][m][0][0] for m in ('ph1_s', 'rssanon_peak', 'arena_peak', 'levels', 'nonunique_share', 'rank_equals_oracle', 'scratch_peak') if m in ex[c]})
