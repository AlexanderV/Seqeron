#!/usr/bin/env python3
"""Medians and ratios of the G3 / G5 rows of work package 1b.7e (rows_1b7e_g3.tsv, rows_1b7e_g5.tsv and any extra rows file given): per cell (input, cell, geometry) the median of the runs of every variant, PH6 production against the prototype, the
container against the sum of the separate phases (+ the checksum pass), the counted I/O, the scratch peak and RssAnon (in-process and external sampler). usage: summarize_1b7e_g3.py <rows file>... """
import statistics as st, sys, collections

CONT = 'label variant scan create pack ph1 ph23 ph45 bwt psi crc seal total read write mapped scratch_peak arena_peak rss_ph0 rss_ph1 rss_ph23 rss_ph45 rss_ph6 rss_seal rss_whole verdict file_size grown reservation overflow lcpkind rss_before_kb rss_max_kb samples'.split()
SEP = 'label variant ph0 ph1 ph23 ph45 bwt psi crc total read write mapped scratch_peak bytes copy_s copied rss_before_kb rss_max_kb samples'.split()
P6 = 'label variant bwt_s psi_s read1 write1 mapped1 read2 write2 mapped2 windows primary rss_before_kb rss_max_kb samples'.split()
P6PROD = 'label variant bwt_s psi_s read1 write1 mapped1 read2 write2 mapped2 windows primary rss_before_kb rss_max_kb samples'.split()


def num(x):
    try:
        return float(x)
    except ValueError:
        return x


rows = collections.defaultdict(list)  # (input, cell, variant) -> [dict]
ext = {}
for path in sys.argv[1:]:
    for line in open(path):
        p = line.rstrip('\n').split('\t')
        if line.startswith('#ext'):
            ext[p[1]] = (num(p[3]), num(p[4]))  # rss_max_kb (the test host), samples
            continue
        if line.startswith('#') or len(p) < 3:
            continue
        label, variant = p[0], p[1]
        inp, cell, name, _, _ = (label.split('.') + [''] * 5)[:5]
        variant = name  # the third part of the label names the variant of the run (the column repeats the family: containernopsi is a container)
        if variant.startswith('container'):
            d = dict(zip(CONT, map(num, p)))
        elif variant.startswith('separate'):
            d = dict(zip(SEP, map(num, p)))
        else:
            d = dict(zip(P6, map(num, p)))
        d['label'] = label
        d['ext_kb'] = ext.get(label, (None, None))[0]
        rows[(inp, cell, variant)].append(d)
# the #ext line follows its row: attach afterwards
for lst in rows.values():
    for d in lst:
        d['ext_kb'] = ext.get(d['label'], (None, None))[0]


def med(lst, key):
    v = [d[key] for d in lst if key in d and isinstance(d[key], float)]
    return st.median(v) if v else float('nan')


def rng(lst, key):
    v = [d[key] for d in lst if key in d and isinstance(d[key], float)]
    return (min(v), max(v)) if v else (float('nan'), float('nan'))


def pct(a, b):
    return 100.0 * (a - b) / b if b else float('nan')


for (inp, cell, variant), lst in sorted(rows.items()):
    n = len(lst)
    if variant.startswith('ph6'):
        print('%-8s %-9s %-14s n=%d  bwt %.3f [%.3f-%.3f]  psi %.3f [%.3f-%.3f]  io r/w/m %s' % (inp, cell, variant, n, med(lst, 'bwt_s'), *rng(lst, 'bwt_s'), med(lst, 'psi_s'), *rng(lst, 'psi_s'), '/'.join('%d' % med(lst, k) for k in ('read1', 'write1', 'mapped1'))))
    elif variant.startswith('container'):
        print('%-8s %-9s %-14s n=%d  total %.2f [%.2f-%.2f]  scan %.3f create %.4f pack %.3f ph1 %.2f ph23 %.2f ph45 %.2f bwt %.2f psi %.2f crc %.3f seal %.4f | io r/w/m %d/%d/%d scratch %d arena %d | rss whole %.1f ext %.1f MB' % (
            inp, cell, variant, n, med(lst, 'total'), *rng(lst, 'total'), med(lst, 'scan'), med(lst, 'create'), med(lst, 'pack'), med(lst, 'ph1'), med(lst, 'ph23'), med(lst, 'ph45'), med(lst, 'bwt'), med(lst, 'psi'), med(lst, 'crc'), med(lst, 'seal'),
            med(lst, 'read'), med(lst, 'write'), med(lst, 'mapped'), med(lst, 'scratch_peak'), med(lst, 'arena_peak'), med(lst, 'rss_whole'), (med(lst, 'ext_kb') or float('nan')) * 1024 / 1e6))
    else:
        print('%-8s %-9s %-14s n=%d  total %.2f [%.2f-%.2f]  ph0 %.3f ph1 %.2f ph23 %.2f ph45 %.2f bwt %.2f psi %.2f crc %.3f copy %.3f | io r/w/m %d/%d/%d scratch %d | ext %.1f MB' % (
            inp, cell, variant, n, med(lst, 'total'), *rng(lst, 'total'), med(lst, 'ph0'), med(lst, 'ph1'), med(lst, 'ph23'), med(lst, 'ph45'), med(lst, 'bwt'), med(lst, 'psi'), med(lst, 'crc'), med(lst, 'copy_s'),
            med(lst, 'read'), med(lst, 'write'), med(lst, 'mapped'), med(lst, 'scratch_peak'), (med(lst, 'ext_kb') or float('nan')) * 1024 / 1e6))
print()
print('--- comparisons (medians) ---')
for (inp, cell, variant), lst in sorted(rows.items()):
    if variant == 'ph6prod':
        pr = rows.get((inp, cell, 'ph6proto'))
        if pr:
            for k in ('bwt_s', 'psi_s'):
                print('%-8s %-9s PH6 %-5s production %.3f vs prototype %.3f: %+.1f %%' % (inp, cell, k, med(lst, k), med(pr, k), pct(med(lst, k), med(pr, k))))
            for k in ('read1', 'write1', 'mapped1'):
                print('%-8s %-9s PH6 %-8s production %d vs prototype %d: %+.3f %%' % (inp, cell, k, med(lst, k), med(pr, k), pct(med(lst, k), med(pr, k))))
    if variant == 'container':
        sp = rows.get((inp, cell, 'separate'))
        if sp:
            print('%-8s %-9s container total %.2f vs separate %.2f (incl. crc): %+.2f %%;  scratch peak %d vs %d;  write %d vs %d (%+.4f %%)' % (inp, cell, med(lst, 'total'), med(sp, 'total'), pct(med(lst, 'total'), med(sp, 'total')), med(lst, 'scratch_peak'), med(sp, 'scratch_peak'), med(lst, 'write'), med(sp, 'write'), pct(med(lst, 'write'), med(sp, 'write'))))
for (inp, cell, variant), lst in sorted(rows.items()):
    if variant in ('containernopsi', 'containersparse', 'containernosync'):
        base = rows.get((inp, cell, 'container'))
        if base:
            print('%-8s %-9s %s vs container: total %.2f vs %.2f (%+.2f %%), write %d vs %d, mapped %d vs %d, file %d vs %d' % (inp, cell, variant, med(lst, 'total'), med(base, 'total'), pct(med(lst, 'total'), med(base, 'total')), med(lst, 'write'), med(base, 'write'), med(lst, 'mapped'), med(base, 'mapped'), med(lst, 'file_size'), med(base, 'file_size')))
    if variant == 'separatecopy':
        base = rows.get((inp, cell, 'separate'))
        if base:
            print('%-8s %-9s separate + copy: total %.2f vs separate %.2f; copy %.3f s, %d bytes' % (inp, cell, med(lst, 'total'), med(base, 'total'), med(lst, 'copy_s'), med(lst, 'copied')))
