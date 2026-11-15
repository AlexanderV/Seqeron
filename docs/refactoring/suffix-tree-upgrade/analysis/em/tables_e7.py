#!/usr/bin/env python3
"""Tables of the raw rows of E7: python3 -I tables_e7.py <mode> <rows tsv>... 
Rows are R <TAB> E7 <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes (rows without a numeric run column are skipped).
Cases are  <family>_<size>m.e7.<series>[.t<threads>]  with series A (fixed parameters, --rule-n 10^8) or B (ratio emulation, --scale-n 3x10^9); a case without a thread suffix ran with 4 threads.
Modes:
  fit      per (family, series): total wall time per size (median [min-max]), the fit t = a n^b by least squares in log-log on the medians (and on the per-size min and on the per-size max), per-phase exponents
           (PH0, PH1, PH2/PH3, PH4+PH5, and the fused total), local exponents between neighbouring sizes, ns per symbol, and the ESTIMATE for 3x10^9 (a n^b with the three fits)
  rss      peak RssAnon per size (program sampler and independent sampler), per-phase peaks, and the slope of the median peak against log2(n) (least squares over the sizes, MB per doubling)
  io       I/O B/char (program counters, rchar + wchar), mapped writes, scratch peak, output footprint, block-layer reads, dirty pages per size
  cmp      per case: how many identity comparisons (cmp_*, exit_code, scratch_left) are yes / NO
  threads  total time and RssAnon of the 1 / 2 thread runs next to the 4 thread runs
  cache    the cache sweep (cases <name>.e7.C.<cell> and <name>.e7.calib.<cell>): per cell the runs, total time, slowdown against the unlimited cell, per-phase times, the effective cache (resident fraction of text + rank
           file, time-weighted over PH3), block-layer reads, memory.failcnt, kills, RssAnon, identity
Options --max-size <M>  (fit, rss): only the sizes up to M million symbols;  --kind C|calib  (cache): only the sweep cells or only the calibration cells.
"""
import math, re, statistics, sys
from collections import defaultdict

mode = sys.argv[1]
args = sys.argv[2:]
MAX_SIZE = None   # --max-size <M>: the fits use only the sizes up to M million symbols (genome-like: 200, the sizes the task fixed; 500 M was added on top)
if '--max-size' in args:
    i = args.index('--max-size')
    MAX_SIZE = int(args[i + 1])
    del args[i:i + 2]
KIND = None       # --kind C | calib (cache): only the cells of that kind are printed (the baselines still come from all the files)
if '--kind' in args:
    i = args.index('--kind')
    KIND = args[i + 1]
    del args[i:i + 2]
paths = args
data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))   # case -> run -> metric -> [values]
order = []
for path in paths:
    for line in open(path):
        c = line.rstrip('\n').split('\t')
        if len(c) < 7 or c[0] != 'R' or c[1] != 'E7' or not c[3].isdigit():
            continue
        case, run, metric, value = c[2], c[3], c[4], c[5]
        if case not in order:
            order.append(case)
        data[case][run][metric].append(value)

CASE = re.compile(r'^(?P<name>(?P<fam>[a-z]+)_(?P<size>\d+)m)\.e7\.(?P<series>[AB])(\.t(?P<th>\d))?$')

def info(case):
    m = CASE.match(case)
    if not m:
        return None
    return m.group('fam'), int(m.group('size')), m.group('series'), int(m.group('th') or 4)

def runs(case):
    return [r for r in sorted(data[case], key=int) if 'run_done' in data[case][r] or 'total_s' in data[case][r]]

def vals(case, metric, scale=1.0):
    out = []
    for r in runs(case):
        for x in data[case][r].get(metric, []):
            try:
                out.append(float(x) * scale)
            except ValueError:
                pass
    return out

def n_of(case):
    v = vals(case, 'n')
    return v[0] if v else None

def med(v):
    return statistics.median(v) if v else float('nan')

def fmt(v, d=4):
    if not v:
        return '-'
    if len(v) == 1:
        return ('%.' + str(d) + 'g') % v[0]
    return ('%.' + str(d) + 'g [%.' + str(d) + 'g-%.' + str(d) + 'g]') % (med(v), min(v), max(v))

def lsq(xs, ys):
    """least squares y = a + b x: (a, b)"""
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx if sxx > 0 else float('nan')
    return my - b * mx, b

def powfit(ns, ts):
    """t = a n^b by least squares in log-log: (a, b)"""
    a, b = lsq([math.log(n) for n in ns], [math.log(t) for t in ts])
    return math.exp(a), b

def groups():
    g = defaultdict(list)
    for case in order:
        i = info(case)
        if i and i[3] == 4 and runs(case) and (MAX_SIZE is None or i[1] <= MAX_SIZE):
            g[(i[0], i[2])].append((i[1], case))
    for k in g:
        g[k].sort()
    return g

PH = [('PH0', 'ph0_prepare_s'), ('PH1', 'ph1_dcs_s'), ('PH2/PH3', 'ph23_sa_s'), ('PH4+PH5', 'ph45_lcp_child_s'), ('total', 'total_s')]

if mode == 'fit':
    for (fam, ser), cases in sorted(groups().items()):
        ns = [n_of(c) for _, c in cases]
        print('\n### %s series %s (4 threads, cold cache, runs per size: %s)\n' % (fam, ser, ', '.join(str(len(runs(c))) for _, c in cases)))
        print('| size | n | total s | PH0 s | PH1 s | PH2/PH3 s | PH4+PH5 s | ns / symbol | local exponent b (previous size) |')
        print('|---|---|---|---|---|---|---|---|---|')
        prev = None
        for (sz, c), n in zip(cases, ns):
            tot = vals(c, 'total_s')
            loc = '-' if prev is None else '%.3f' % (math.log(med(tot) / prev[1]) / math.log(n / prev[0]))
            print('| %d M | %d | %s | %s | %s | %s | %s | %s | %s |' % (sz, n, fmt(tot), fmt(vals(c, 'ph0_prepare_s')), fmt(vals(c, 'ph1_dcs_s')), fmt(vals(c, 'ph23_sa_s')), fmt(vals(c, 'ph45_lcp_child_s')), '%.0f' % (med(tot) * 1e9 / n), loc))
            prev = (n, med(tot))
        if len(cases) >= 2:
            print('\nExponents b of t = a n^b (least squares in log-log over %d sizes):\n' % len(cases))
            print('| series | b (medians) | b (per-size min) | b (per-size max) | a (medians) | ESTIMATE t(3x10^9) from the medians, min fit, max fit |')
            print('|---|---|---|---|---|---|')
            rows = []
            fits = {}
            for label, f in (('medians', med), ('min', min), ('max', max)):
                ts = [f(vals(c, 'total_s')) for _, c in cases]
                fits[label] = powfit(ns, ts)
            est = ['%.0f s' % (fits[k][0] * 3e9 ** fits[k][1]) for k in ('medians', 'min', 'max')]
            print('| %s %s | **%.3f** | %.3f | %.3f | %.3g | %s (ESTIMATE, not measured) |' % (fam, ser, fits['medians'][1], fits['min'][1], fits['max'][1], fits['medians'][0], ' / '.join(est)))
            print('\nPer-phase exponents (medians):\n')
            print('| phase | b | share of the total at the largest size |')
            print('|---|---|---|')
            tot_last = med(vals(cases[-1][1], 'total_s'))
            for label, key in PH:
                ts = [med(vals(c, key)) for _, c in cases]
                if any(not (t > 0) for t in ts):
                    continue
                a, b = powfit(ns, ts)
                print('| %s | %.3f | %.1f %% |' % (label, b, 100 * ts[-1] / tot_last))
            a, b = powfit(ns, [med(vals(c, 'ph23_sa_s')) + med(vals(c, 'ph1_dcs_s')) for _, c in cases])
            print('| PH1+PH2/PH3 | %.3f | |' % b)
elif mode == 'rss':
    for (fam, ser), cases in sorted(groups().items()):
        print('\n### %s series %s: RssAnon (MB, median [min-max])\n' % (fam, ser))
        print('| size | program sampler (10 ms) | independent sampler (20 ms) | PH0 | PH1 | PH2/PH3 | PH4+PH5 | arena peak MB |')
        print('|---|---|---|---|---|---|---|---|')
        xs, ys, ye = [], [], []
        for sz, c in cases:
            p = vals(c, 'rssanon_peak')
            e = vals(c, 'ext_rssanon_peak')
            xs.append(math.log2(n_of(c)))
            ys.append(med(p))
            ye.append(med(e) if e else float('nan'))
            print('| %d M | %s | %s | %s | %s | %s | %s | %s |' % (sz, fmt(p, 4), fmt(e, 4), fmt(vals(c, 'ph0_prepare_rssanon_peak'), 3), fmt(vals(c, 'ph1_dcs_rssanon_peak'), 3), fmt(vals(c, 'ph23_sa_rssanon_peak'), 3), fmt(vals(c, 'ph45_lcp_child_rssanon_peak'), 3), fmt(vals(c, 'arena_peak'), 3)))
        if len(xs) >= 2:
            a, b = lsq(xs, ys)
            print('\nSlope of the median peak (program sampler) against log2(n): **%.2f MB per doubling** (least squares over %d sizes; the last doubling alone: %.2f MB)' % (b, len(xs), ys[-1] - ys[-2] if len(xs) >= 2 else float('nan')))
            if not any(math.isnan(v) for v in ye):
                a2, b2 = lsq(xs, ye)
                print('Slope with the independent sampler: %.2f MB per doubling' % b2)
            if len(xs) >= 3:
                print('Doublings step by step: ' + ', '.join('%d M -> %d M: %+.1f MB' % (cases[i][0], cases[i + 1][0], (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])) for i in range(len(xs) - 1)))
elif mode == 'io':
    print('| case | n | I/O B/char (program counters) | rchar + wchar B/char (program / independent whole process) | mapped write B/char | scratch peak B/char (tracked / apparent / allocated) | output footprint B/char | block-layer read B/char | dirty peak MB (system / mapped outputs) | disk free at start GB |')
    print('|---|---|---|---|---|---|---|---|---|---|')
    for case in order:
        i = info(case)
        if not i or not runs(case):
            continue
        n = n_of(case)
        pc = lambda m: [x / n for x in vals(case, m)]
        ext = [(a + b) / n for a, b in zip(vals(case, 'ext_proc_rchar'), vals(case, 'ext_proc_wchar'))]
        print('| %s | %d | %s | %s / %s | %s | %s / %s / %s | %s | %s | %s / %s | %s |' % (case, n, fmt(pc('io_total_counted'), 4), fmt(pc('proc_rw_total'), 4), fmt(ext, 4), fmt(pc('io_mapped_write'), 3), fmt(pc('scratch_peak'), 4),
              fmt(pc('ext_scratch_apparent_peak'), 4), fmt(pc('ext_scratch_allocated_peak'), 4), fmt(pc('output_bytes'), 4), fmt(pc('blocklayer_read'), 3), fmt(vals(case, 'dirty_peak'), 3), fmt(vals(case, 'proc_dirty_peak'), 3), fmt(vals(case, 'disk_free_at_start'), 3)))
elif mode == 'cmp':
    for case in order:
        yes = no = 0
        names = []
        for r in data[case]:
            for metric, vs in data[case][r].items():
                if metric.startswith('cmp_') or metric in ('exit_code', 'scratch_left'):
                    for v in vs:
                        if v in ('yes', '0'):
                            yes += 1
                        else:
                            no += 1
                            names.append('run%s:%s=%s' % (r, metric, v))
        print('%-32s runs=%d yes=%d NO=%d %s' % (case, len(runs(case)), yes, no, ' '.join(names)))
elif mode == 'threads':
    print('| input | series | threads | n | total s | PH2/PH3 s | PH4+PH5 s | RssAnon MB (program) | runs |')
    print('|---|---|---|---|---|---|---|---|---|')
    rowsout = []
    for case in order:
        i = info(case)
        if not i or not runs(case):
            continue
        rowsout.append((i[0], i[2], i[1], i[3], case))
    for fam, ser, sz, th, case in sorted(rowsout, key=lambda x: (x[0], x[1], x[2], -x[3])):
        print('| %s_%dm | %s | %d | %d | %s | %s | %s | %s | %d |' % (fam, sz, ser, th, n_of(case), fmt(vals(case, 'total_s')), fmt(vals(case, 'ph23_sa_s')), fmt(vals(case, 'ph45_lcp_child_s')), fmt(vals(case, 'rssanon_peak')), len(runs(case))))
elif mode == 'cache':
    CC = re.compile(r'^(?P<name>[a-z]+_\d+m)\.e7\.(?P<kind>C|calib)\.(?P<cell>.+)$')
    base = {}
    base3 = {}
    for case in order:
        m = CC.match(case)
        if m and m.group('kind') == 'C' and m.group('cell') == 'U':
            t = [x for x in vals(case, 'total_s')]
            if t:
                base[m.group('name')] = med(t)
                base3[m.group('name')] = med(vals(case, 'ph23_sa_s'))
    for case in order:
        m = CC.match(case)
        if m and m.group('kind') == 'calib' and m.group('cell') == 'C1' and m.group('name') not in base:
            base[m.group('name')] = med(vals(case, 'total_s'))
            base3[m.group('name')] = med(vals(case, 'ph23_sa_s'))
    print('| input | cell | limit MB (raised to MB after PH3) | runs finished / started | total s | slowdown (median) | PH1 s | PH2/PH3 s | PH2/PH3 slowdown | PH4+PH5 s | effective cache (PH3 resident fraction of text + rank) | block-layer read MB (whole process) | memory.failcnt | killed | RssAnon peak MB | identity |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for case in order:
        m = CC.match(case)
        if not m or (KIND is not None and m.group('kind') != KIND):
            continue
        rs = [r for r in sorted(data[case], key=int) if 'run_done' in data[case][r]]
        if not rs:
            continue
        ok = [r for r in rs if data[case][r].get('exit_code', ['?'])[0] == '0' and 'total_s' in data[case][r]]
        def v(metric, only=None, scale=1.0):
            out = []
            for r in (only if only is not None else rs):
                for x in data[case][r].get(metric, []):
                    try:
                        out.append(float(x) * scale)
                    except ValueError:
                        pass
            return out
        st = data[case][rs[0]].get('settings', [''])[0]
        lim = re.search(r'limit=(\S+)', st).group(1)
        rzm = re.search(r'raise=(\S+)', st)
        rz = rzm.group(1) if rzm else 'none'
        limtxt = 'unlimited' if lim == 'max' else '%.0f' % (int(lim) / 1048576)
        if rz not in ('none', ''):
            limtxt += ' (%.0f)' % (int(rz) / 1048576)
        tot = v('total_s', ok)
        slow = [t / base[m.group('name')] for t in tot] if m.group('name') in base else []
        killed = sorted(set(data[case][r].get('killed', ['?'])[0] for r in rs))
        idn = []
        for r in ok:
            idn.append(all(data[case][r].get(k, ['?'])[0] == 'yes' for k in ('cmp_pk', 'cmp_sa', 'cmp_lcp', 'cmp_ovf', 'cmp_cld', 'cmp_stats')))
        idtxt = ('%d of %d runs equal to the references' % (sum(idn), len(idn))) if idn else 'no finished run'
        tottxt = fmt(tot) if tot else '> watchdog'
        slow3 = [t / base3[m.group('name')] for t in v('ph23_sa_s', ok)] if m.group('name') in base3 else []
        print('| %s | %s | %s | %d / %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
            m.group('name'), m.group('cell') if m.group('kind') == 'C' else 'calib ' + m.group('cell'), limtxt, len(ok), len(rs), tottxt, fmt(slow, 3), fmt(v('ph1_dcs_s', ok), 3), fmt(v('ph23_sa_s', ok), 3), fmt(slow3, 3), fmt(v('ph45_lcp_child_s', ok), 3),
            fmt(v('eff_resident_ph3'), 3), fmt(v('ext_proc_read_bytes', None, 1e-6), 4), fmt(v('cg_failcnt_final'), 4), '/'.join(k for k in killed), fmt(v('ext_rssanon_peak'), 3), idtxt))
else:
    print(__doc__)
