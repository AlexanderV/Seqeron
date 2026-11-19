#!/usr/bin/env python3
"""Tables of gate G3 / G5 of work package 1b.7f from rows_1b7f_g3.tsv (and rows_1b7f_g5.tsv, rows_1b7f_g2.tsv): the cold Header open (10 runs), the levels at 10 Mbp, the throughput of the file against the engine in memory
(medians over the rounds; one tree per process, alternating orders), and, for G5, the time and the memory of the structure level against the chunk of its Phi pass. usage: summarize_1b7f_g3.py [directory of the rows files]"""
import os, statistics as st, sys

D = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))


def rows(name, kind):
    for line in open(os.path.join(D, name)):
        p = line.rstrip('\n').split('\t')
        if p[0] == kind:
            yield p


op, warm, hdr, lv, rss = [], [], [], {'Sections': [], 'Structure': []}, []
for p in rows('rows_1b7f_g3.tsv', 'open'):
    if p[1].startswith('cold-open'):
        op.append(float(p[4]) * 1000)
        hdr.append(float(p[5]) * 1000)
    elif p[1].startswith('levels'):
        lv[p[2]].append(float(p[4]))
        if p[2] == 'Structure':
            rss.append((int(p[13]) - int(p[12])) / 1024)
for p in rows('rows_1b7f_g3.tsv', 'open-warm'):
    warm.append(float(p[4]) * 1000)
print('cold Header open (ms, 10 runs): median %.2f, min %.2f, max %.2f, %d of %d at most 10 ms; header level alone %.2f; warm second open %.2f-%.2f' % (st.median(op), min(op), max(op), sum(1 for x in op if x <= 10), len(op), st.median(hdr), min(warm), max(warm)))
for k, v in lv.items():
    print('%s at 10 Mbp (s): %s median %.3f' % (k, ' '.join('%.3f' % x for x in v), st.median(v)))
print('RssAnon of Structure above the process baseline (MB): %s' % ' '.join('%.1f' % x for x in rss))
t = {}
for p in rows('rows_1b7f_g3.tsv', 'throughput'):
    t.setdefault(p[2], []).append([float(x) for x in p[3:7]])
names = ['count20/s', 'locate8/s', 'count24 absent/s', 'child walk nodes/s']
for tree, v in t.items():
    print('%-8s (%d processes) median %s' % (tree, len(v), ' '.join('%s %d' % (n, st.median(c[i] for c in v)) for i, n in enumerate(names))))
for tree in ('checked', 'trusted'):
    print('%s against memory (median, %%): %s' % (tree, ' '.join('%.1f' % (100 * (st.median(c[i] for c in t[tree]) / st.median(c[i] for c in t['memory']) - 1)) for i in range(4))))
print('checked against trusted (median, %%): %s' % ' '.join('%.1f' % (100 * (st.median(c[i] for c in t['checked']) / st.median(c[i] for c in t['trusted']) - 1)) for i in range(4)))
if os.path.exists(os.path.join(D, 'rows_1b7f_g5.tsv')):
    for p in rows('rows_1b7f_g5.tsv', 'open'):
        print('G5 %s: %.1f s, RssAnon above baseline %.1f MB, buffers %.1f MB' % (p[1], float(p[4]), (int(p[13]) - int(p[12])) / 1024, int(p[11]) / 1e6))
