#!/usr/bin/env python3
"""Tables of gate G3 / G5 of work package 1b.7d from the raw rows (rows_1b7d_g3.tsv): per cell (input, geometry, cache) and variant the medians of the time of PH1, of PH2/PH3 and of the whole call, the counted I/O,
   the scratch peaks, the in-process RssAnon peak (maximum over the runs of the cell) and the independent sampler's, and the ratio production / prototype per cell. usage: tables_1b7d_g3.py <rows file>"""
import sys, statistics, collections

rows = [l.rstrip('\n').split('\t') for l in open(sys.argv[1]) if l.strip()]
data = collections.defaultdict(list)
ext = {}
cmp_ = collections.defaultdict(list)
for r in rows:
    if r[0] == '#ext':
        ext[(r[1], r[2])] = int(r[3])
    elif r[0] == '#cmp':
        cmp_[r[1]].append(r[3])
    elif r[0].startswith('#'):
        continue
    else:
        data[(r[0].rsplit('.', 1)[0], r[1])].append(r)  # label without the run number


def med(rs, i):
    return statistics.median(float(r[i]) for r in rs)


cells = sorted({k[0] for k in data})
print('cell\tvariant\tn\tph1_s\tph23_s\ttotal_s\tph1_rd\tph1_wr\tph1_map\tph23_rd\tph23_wr\tph23_map\tscr_ph1\tscr_all\trss_in_MB\trss_ext_MB\tverdicts')
for c in cells:
    base = {}
    for v in ('proto', 'prod', 'prodnr', 'prodsparse'):
        rs = data.get((c, v))
        if not rs:
            continue
        rss_in = max(int(r[18]) for r in rs) * 1024 / 1e6
        ex = [x for (l, vv), x in ext.items() if vv == v and l.startswith(c + '.')]
        rss_ex = max(ex) * 1024 / 1e6 if ex else float('nan')
        vd = sorted({r[16] for r in rs})
        t = [med(rs, 2), med(rs, 3), med(rs, 4)]
        base[v] = t
        print('\t'.join([c, v, str(len(rs)), '%.2f' % t[0], '%.2f' % t[1], '%.2f' % t[2]] + ['%d' % med(rs, i) for i in (5, 6, 7, 8, 9, 10, 11, 12)] + ['%.1f' % rss_in, '%.1f' % rss_ex, ','.join(vd)]))
    if 'proto' in base:
        for v in ('prod', 'prodnr', 'prodsparse'):
            if v in base:
                print('  ratio %s / proto: PH1 %.3f  PH2/PH3 %.3f  total %.3f' % (v, base[v][0] / base['proto'][0], base[v][1] / base['proto'][1], base[v][2] / base['proto'][2]))
    elif 'prod' in base:
        for v in ('prodnr', 'prodsparse'):
            if v in base:
                print('  ratio %s / prod: PH1 %.3f  PH2/PH3 %.3f  total %.3f' % (v, base[v][0] / base['prod'][0], base[v][1] / base['prod'][1], base[v][2] / base['prod'][2]))
for k, v in sorted(cmp_.items()):
    pass
bad = [(k, x) for k, v in cmp_.items() for x in v if x != 'equal-to-reference']
print('cmp rows', sum(len(v) for v in cmp_.values()), 'not equal:', bad)
