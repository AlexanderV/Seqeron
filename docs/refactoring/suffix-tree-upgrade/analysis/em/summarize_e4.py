#!/usr/bin/env python3
"""Summaries of the E4 raw rows (analysis/em/e4_*.tsv): one markdown table per question of the gates G1-G4 of the work package.
usage: summarize_e4.py <analysis/em dir> [section]   (sections: g1 phases slope setup big all)
Row format: R <TAB> E4 <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes (the run column exists when the command was given --run)."""
import collections, glob, math, os, re, statistics, sys

D = sys.argv[1] if len(sys.argv) > 1 else '.'
SECTION = sys.argv[2] if len(sys.argv) > 2 else 'all'


def load(pattern):
    rows = collections.defaultdict(lambda: collections.defaultdict(dict))  # case -> run -> metric -> (value, notes)
    for path in sorted(glob.glob(os.path.join(D, pattern))):
        for line in open(path):
            c = line.rstrip('\n').split('\t')
            if len(c) < 7 or c[0] != 'R':
                continue
            if c[1] not in ('E4', 'E3'):
                continue
            if len(c) >= 8:
                _, exp, case, run, metric, value, unit, notes = c[:8]
            else:
                _, exp, case, metric, value, unit, notes = c[:7]
                run = '1'
            rows[case][run][metric] = (value, unit, notes)
    return rows


def num(x):
    try:
        return float(x)
    except ValueError:
        return None


def series(rows, case, metric):
    out = []
    for run, m in sorted(rows[case].items()):
        if metric in m:
            v = num(m[metric][0])
            if v is not None:
                out.append(v)
    return out


def med(xs, fmt='{:.2f}'):
    if not xs:
        return '-'
    s = fmt.format(statistics.median(xs))
    if len(xs) > 1:
        s += ' [' + fmt.format(min(xs)) + '-' + fmt.format(max(xs)) + ']'
    return s


def per_char(rows, case, metric):
    """metric (bytes) divided by n, per run"""
    out = []
    for run, m in sorted(rows[case].items()):
        if metric in m and 'n' in m:
            out.append(float(m[metric][0]) / float(m['n'][0]))
    return out


def ratio_series(rows, case, metric):
    out = []
    for run, m in sorted(rows[case].items()):
        if metric in m and 'n' in m:
            out.append(float(m[metric][0]) * 1e9 / float(m['n'][0]))
    return out


def io_per_char(rows, case):
    """explicit I/O + mapped writes per character, plus the 4 B/char the bounded mode's gather reads through a read-only mapping"""
    out = []
    for run, m in sorted(rows[case].items()):
        if 'io_total_counted' in m and 'n' in m:
            v = float(m['io_total_counted'][0]) / float(m['n'][0])
            if 'mode=Bounded' in m['n'][2]:
                v += 4.0
            out.append(v)
    return out


def width_of(rows, case):
    for run, m in sorted(rows[case].items()):
        if 'width' in m:
            return m['width'][0]
    return '-'


def verdicts(rows, case):
    y = n = 0
    for run, m in rows[case].items():
        if 'lcp_equals_oracle' in m:
            if m['lcp_equals_oracle'][0] == 'yes':
                y += 1
            else:
                n += 1
    return y, n


def g1():
    allrows = load('e4_*.tsv')
    yes = no = 0
    cases = 0
    for case in allrows:
        y, n = verdicts(allrows, case)
        yes += y
        no += n
        cases += 1 if (y + n) else 0
    print(f'LCP equals the LcpKernel oracle (byte-identical base and overflow files): {yes} of {yes + no} protocol runs, {no} different, {cases} case labels')


def phases(group):
    rows = load(f'e4_{group}.tsv')
    print(f'\n### {group}\n')
    print('| case | runs | equal | PH4 s | B1 s | B2 s | B3 s | ns/symbol | Msym/s | RssAnon MB | RssFile MB | I/O B/char (explicit + mapped + gather) | scratch B/char | width | escapes >= 255 % | restart symbols |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for case in sorted(rows):
        if not series(rows, case, 'ph4_s'):
            continue
        y, n = verdicts(rows, case)
        esc = []
        restart = []
        for run, m in sorted(rows[case].items()):
            if 'width' in m:
                mm = re.search(r'escapes >= 255: \d+ \(([\d.]+) %\)', m['width'][2])
                if mm:
                    esc.append(float(mm.group(1)))
            if 'compare_steps' in m:
                mm = re.search(r'no carry\): (\d+)', m['compare_steps'][2])
                if mm:
                    restart.append(float(mm.group(1)))
        print('| ' + ' | '.join([
            case, str(len(series(rows, case, 'ph4_s'))), f'{y}/{y + n}',
            med(series(rows, case, 'ph4_s')), med(series(rows, case, 'b1_s')), med(series(rows, case, 'b2_s')), med(series(rows, case, 'b3_s')),
            med(series(rows, case, 'ns_per_symbol'), '{:.1f}'), med(series(rows, case, 'throughput'), '{:.1f}'),
            med(series(rows, case, 'rssanon_peak'), '{:.1f}'), med(series(rows, case, 'rssfile_peak'), '{:.1f}'),
            med(io_per_char(rows, case), '{:.2f}'),
            med(per_char(rows, case, 'scratch_peak'), '{:.3f}'), width_of(rows, case),
            med(esc[:1], '{:.4f}'), med(restart[:1], '{:.0f}'),
        ]) + ' |')


def slope(label_re, sizes_of):
    rows = load('e4_slope.tsv')
    pts = collections.defaultdict(list)
    for case in rows:
        mm = re.fullmatch(label_re, case)
        if not mm:
            continue
        for run, m in rows[case].items():
            if 'rssanon_peak' in m and 'n' in m:
                pts[int(float(m['n'][0]))].append(float(m['rssanon_peak'][0]))
    if len(pts) < 2:
        return None
    ns = sorted(pts)
    xs = [math.log2(n) for n in ns]
    ys = [statistics.median(pts[n]) for n in ns]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    allx = [math.log2(n) for n in ns for _ in pts[n]]
    ally = [v for n in ns for v in pts[n]]
    ax, ay = statistics.mean(allx), statistics.mean(ally)
    ball = sum((x - ax) * (y - ay) for x, y in zip(allx, ally)) / sum((x - ax) ** 2 for x in allx)
    steps = [(ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]) for i in range(len(ns) - 1)]
    return ns, pts, b, ball, steps


def slopes():
    print('\n### slope of RssAnon against n (scale-down emulation, 4 threads, cold cache)\n')
    print('| series | RssAnon peak per n (median [min-max], runs) | least-squares slope MB per doubling (all runs) | step slopes |')
    print('|---|---|---|---|')
    for name, rx in [
        ('random DNA, base mode, plain, 10-100 M', r'(dna10r|dna50r|dna_100m)\.dist\.slope'),
        ('random DNA, base mode, plain, 10-200 M', r'(dna10r|dna50r|dna_100m|dna200r)\.dist\.slope'),
        ('random DNA, base mode, touched, 10-100 M', r'(dna10r|dna50r|dna_100m)\.dist\.slopetouch'),
        ('random DNA, base mode, touched, 10-200 M', r'(dna10r|dna50r|dna_100m|dna200r)\.dist\.slopetouch'),
        ('random DNA, bounded (2 rounds), plain, 10-100 M', r'(dna10r|dna50r|dna_100m)\.b2\.slope'),
        ('random DNA, bounded (2 rounds), plain, 10-200 M', r'(dna10r|dna50r|dna_100m|dna200r)\.b2\.slope'),
        ('random DNA, bounded (2 rounds), touched, 10-100 M', r'(dna10r|dna50r|dna_100m)\.b2\.slopetouch'),
        ('random DNA, bounded (2 rounds), touched, 10-200 M', r'(dna10r|dna50r|dna_100m|dna200r)\.b2\.slopetouch'),
        ('genome-like, base mode, plain, 10-100 M', r'(genome10|genome50|genome_100m)\.dist\.slope'),
        ('genome-like, base mode, touched, 10-100 M', r'(genome10|genome50|genome_100m)\.dist\.slopetouch'),
        ('genome-like, bounded (2 rounds), plain, 10-100 M', r'(genome10|genome50|genome_100m)\.b2\.slope'),
        ('genome-like, bounded (2 rounds), touched, 10-100 M', r'(genome10|genome50|genome_100m)\.b2\.slopetouch'),
    ]:
        r = slope(rx, None)
        if r is None:
            continue
        ns, pts, b, ball, steps = r
        cell = ' ; '.join(f'{n // 1_000_000} M: {statistics.median(pts[n]):.1f} [{min(pts[n]):.1f}-{max(pts[n]):.1f}] ({len(pts[n])})' for n in ns)
        print(f'| {name} | {cell} | {b:.2f} ({ball:.2f}) | {" / ".join(f"{s:.2f}" for s in steps)} |')


def setup():
    rows = load('e4_setup.tsv')
    print('\n### SA phase (E3 builder, 4 threads, cold cache, same session)\n')
    print('| input | runs | SA phase s | ns/symbol |')
    print('|---|---|---|---|')
    for case in sorted(rows):
        if case.endswith('.sa_phase'):
            print(f'| {case[:-9]} | {len(series(rows, case, "sa_s"))} | {med(series(rows, case, "sa_s"))} | {med(series(rows, case, "ns_per_symbol"), "{:.1f}")} |')
    print('\n### in-RAM kernel (LcpKernel, one thread)\n')
    print('| input | kernel s |')
    print('|---|---|')
    for case in sorted(rows):
        if series(rows, case, 'kernel_s'):
            print(f'| {case} | {med(series(rows, case, "kernel_s"))} |')


def ratio():
    setup_rows = load('e4_setup.tsv')
    sa = {}
    for case in setup_rows:
        if case.endswith('.sa_phase'):
            sa[case[:-9]] = statistics.median(series(setup_rows, case, 'sa_s'))
    print('\n### LCP phase against SA phase (abandon rule: LCP phase > 3 x SA phase)\n')
    print('| input | SA phase s (median) | PH4 base mode s | ratio | PH4 bounded 2 rounds s | ratio | bounded 1 round s | ratio |')
    print('|---|---|---|---|---|---|---|---|')
    allrows = {}
    for g in ('g100', 'small'):
        allrows.update(load(f'e4_{g}.tsv'))
    for name in ['dna_100m', 'genome_100m', 'dna10', 'rep10', 'a50', 'fib10', 'tandem', 'nrun10', 'genome10']:
        if name not in sa:
            continue
        cells = [name, f'{sa[name]:.2f}']
        for suffix in ('dist', 'b2', 'b1'):
            xs = series(allrows, f'{name}.{suffix}.t4', 'ph4_s')
            if xs:
                m = statistics.median(xs)
                cells += [f'{m:.2f}', f'{m / sa[name]:.2f}']
            else:
                cells += ['-', '-']
        print('| ' + ' | '.join(cells) + ' |')


def big():
    rows = load('e4_big*.tsv')
    print('\n### real-size runs (unscaled parameters)\n')
    print('| case | runs | equal | PH4 s | B1 | B2 | B3 | ns/symbol | RssAnon MB | RssFile MB | I/O B/char | scratch B/char | width |')
    print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for case in sorted(rows):
        if not series(rows, case, 'ph4_s'):
            continue
        y, n = verdicts(rows, case)
        print('| ' + ' | '.join([case, str(len(series(rows, case, 'ph4_s'))), f'{y}/{y + n}', med(series(rows, case, 'ph4_s')), med(series(rows, case, 'b1_s')), med(series(rows, case, 'b2_s')), med(series(rows, case, 'b3_s')),
                                 med(series(rows, case, 'ns_per_symbol'), '{:.1f}'), med(series(rows, case, 'rssanon_peak'), '{:.1f}'), med(series(rows, case, 'rssfile_peak'), '{:.1f}'),
                                 med(io_per_char(rows, case), '{:.2f}'), med(per_char(rows, case, 'scratch_peak'), '{:.3f}'), width_of(rows, case)]) + ' |')


if SECTION in ('g1', 'all'):
    g1()
if SECTION in ('phases', 'all'):
    for g in ('g100', 'small', 'fam'):
        if glob.glob(os.path.join(D, f'e4_{g}.tsv')):
            phases(g)
if SECTION in ('slope', 'all'):
    slopes()
if SECTION in ('setup', 'all'):
    setup()
    ratio()
if SECTION in ('big', 'all') and glob.glob(os.path.join(D, 'e4_big*.tsv')):
    big()
