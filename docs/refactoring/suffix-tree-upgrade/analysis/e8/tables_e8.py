#!/usr/bin/env python3
"""E8: turns the rows file of the run (R<TAB>E8<TAB>case<TAB>run<TAB>metric<TAB>value<TAB>unit<TAB>notes) into the markdown tables of RESULTS-E8 and judges the pre-registered gates G1-G4.
usage: python3 -I tables_e8.py <rows file> [--out tables.md] [--json gates.json]
Prints the markdown to stdout (and to --out). Exit code 0 always (the gates are reported, not enforced). UNBUILT AND UNTESTED ON A REAL E8 RUN (tested on synthetic rows and on small local runs).

The gates are fixed by the orchestrator before the run and are quoted verbatim:"""
import json, statistics, sys
from collections import defaultdict

GATES = {
    'G1': 'the verifier passes at 1, 2.2 and 3 Gbp and the output SHA-256 are identical between the generous and every limited-cache run of the same text and between thread counts',
    'G2': 'independently sampled RssAnon of the dotnet process at most 256 MB at every size (E7: 124 MB at 500 Mbp)',
    'G3': 'peak scratch at most 12.5 B/char and peak total disk at most 70 GB at 3 Gbp',
    'G4': 'reported, not gated: time per phase and total (ESTIMATE from E7 for 3 Gbp: 11-15 min on the prototype machine, linear), I/O in B/char, and the cache table of (d) with the answer to the open question',
    'abandon': "a run that takes more than 6 times the generous run of the same text is stopped and reported as 'does not finish in practice at this memory'",
}
__doc__ += '\n' + '\n'.join('  %s: %s' % kv for kv in GATES.items())
SIZES = {'1000m': 1_000_000_000, '2200m': 2_200_000_000, '3000m': 3_000_000_000}   # the sizes of the gates (labels n<millions>m)
LIMIT_RSS_MB = 256.0
LIMIT_SCRATCH = 12.5
LIMIT_DISK_BYTES = 70e9
SLOWDOWN = 3.0
PHASES = [('ph0_prepare_s', 'PH0'), ('ph1_dcs_s', 'PH1'), ('ph23_sa_s', 'PH2/3'), ('ph45_lcp_child_s', 'PH4+5')]


def load(path):
    data = defaultdict(dict)       # (case, run) -> metric -> (value, unit, notes)
    for line in open(path, encoding='utf-8', errors='replace'):
        p = line.rstrip('\n').split('\t')
        if len(p) < 7 or p[0] != 'R' or p[1] != 'E8':
            continue
        if len(p) == 7:
            p.insert(3, '')
        case, run, metric, value, unit = p[2], p[3], p[4], p[5], p[6]
        notes = p[7] if len(p) > 7 else ''
        data[(case, run)][metric] = (value, unit, notes)
    return data


def num(d, metric):
    try:
        return float(d[metric][0])
    except (KeyError, ValueError):
        return None


def done_runs(data):
    return {k: v for k, v in data.items() if v.get('run_done', ('',))[0] == 'yes'}


def label_n(case):
    head = case.split('.')[0]
    if head.startswith('n') and head.endswith('m') and head[1:-1].isdigit():
        return int(head[1:-1]) * 1_000_000
    if head.startswith('n') and head[1:].isdigit():
        return int(head[1:])
    return None


def fmt(x, digits=2):
    if x is None:
        return '-'
    if isinstance(x, str):
        return x
    return ('%.' + str(digits) + 'f') % x


def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join('---' for _ in header) + '|']
    out += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
    return '\n'.join(out)


def run_numbers(d, n):
    """the numbers of one finished run, per symbol where the unit says so"""
    r = {}
    r['finished'] = d.get('exit_code', ('',))[0] == '0'
    r['killed'] = d.get('killed', ('no',))[0]
    r['total_s'] = num(d, 'total_s')
    r['wall'] = num(d, 'process_wall_s')
    for k, _ in PHASES:
        r[k] = num(d, k)
    r['ns'] = num(d, 'ns_per_symbol')
    r['rss_prog'] = num(d, 'rssanon_peak')
    r['rss_ext'] = num(d, 'ext_rssanon_peak')
    io = num(d, 'io_total_counted')
    r['io_counted'] = io / n if io is not None else None
    rw = num(d, 'proc_rw_total')
    r['io_rw'] = rw / n if rw is not None else None
    br = num(d, 'ext_proc_read_bytes')
    r['dev_read'] = br / n if br is not None else None
    bw = num(d, 'ext_proc_write_bytes')
    r['dev_write'] = bw / n if bw is not None else None
    sp = num(d, 'scratch_peak')
    r['scratch_prog'] = sp / n if sp is not None else None
    se = num(d, 'ext_scratch_apparent_peak')
    r['scratch_ext'] = se / n if se is not None else None
    ob = num(d, 'output_bytes')
    r['out'] = ob / n if ob is not None else None
    r['disk_peak'] = num(d, 'ext_disk_total_apparent_peak')
    r['eff'] = num(d, 'eff_resident_ph3')
    r['eff_budget'] = num(d, 'eff_budget_ph3')
    r['failcnt'] = num(d, 'cg_failcnt')
    r['oom'] = num(d, 'cg_oom_kill')
    r['cold'] = d.get('cold_start', ('?',))[0]
    r['psi_some_s'] = (num(d, 'psi_memory_some_us') or 0) / 1e6 if 'psi_memory_some_us' in d else None
    return r


def main():
    path = sys.argv[1]
    out_md = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None
    out_json = sys.argv[sys.argv.index('--json') + 1] if '--json' in sys.argv else None
    data = load(path)
    done = done_runs(data)
    md = []
    gates = {}

    # ---- environment
    pre = {}
    for (case, run), d in data.items():
        if case == 'preflight':
            pre.update({k: v[0] for k, v in d.items()})
    md.append('### Environment (preflight rows)\n')
    keys = [('uname_m', 'machine'), ('image_arch', 'image built for'), ('cpuinfo_format', '/proc/cpuinfo format'), ('emulation', 'emulation'), ('cores', 'cores'), ('cpu_limit', 'CPU quota'), ('mem_total', 'memory visible (bytes)'),
            ('cgroup_version', 'cgroup version'), ('cg_memory_limit', 'container memory limit'), ('cg_swap_limit', 'container swap allowance'), ('disk_size', 'disk size (bytes)'), ('disk_avail', 'disk free at preflight (bytes)'),
            ('fs_type', 'file system'), ('io_test', 'write / mmap / msync test'), ('cold_start_possible', 'drop_caches possible'), ('env_framework', '.NET'), ('env_process_arch', 'process architecture'),
            ('env_lce_path', 'suffix comparison path'), ('env_page_size', 'page size'), ('preflight_result', 'preflight result')]
    md.append(table(['item', 'value'], [(label, pre.get(k, '-')) for k, label in keys]))
    md.append('')

    # ---- generous runs
    gen_rows = []
    verify_rows = []
    generous = {}       # size label -> numbers of run 1 (the baseline of step d)
    for (case, run), d in sorted(done.items(), key=lambda kv: (label_n(kv[0][0]) or 0, kv[0][0], kv[0][1])):
        if not case.endswith('.gen'):
            continue
        n = label_n(case) or num(d, 'n')
        r = run_numbers(d, n)
        generous.setdefault(case.split('.')[0], {})[run] = r
        gen_rows.append([case, run, 'yes' if r['finished'] else 'NO (' + r['killed'] + ')', fmt(r['total_s'], 1), fmt(r['wall'], 1)] + [fmt(r[k], 1) for k, _ in PHASES]
                         + [fmt(r['ns'], 0), fmt(r['rss_prog'], 1), fmt(r['rss_ext'], 1), fmt(r['io_counted'], 1), fmt(r['io_rw'], 1), fmt(r['dev_read'], 2), fmt(r['scratch_prog'], 2), fmt(r['scratch_ext'], 2), fmt(r['out'], 2),
                            fmt((r['disk_peak'] or 0) / 1e9 if r['disk_peak'] else None, 1), r['cold']])
        v = {k[len('verify_'):]: val[0] for k, val in d.items() if k.startswith('verify_') and not k.endswith('_s') and not k.startswith('verify_fact') and k not in ('verify_all', 'verify_exit_code', 'verify_process_wall_s')}
        verify_rows.append([case, run, d.get('verify_all', ('not run',))[0], fmt(num(d, 'verify_wall_s'), 1)] + [v.get(name, '-') for name in ('files', 'meta', 'text', 'input', 'sa_permutation', 'sa_order', 'lcp_values', 'lcp_width', 'stats', 'child', 'sha256')]
                           + [d.get('verify_fact_long_pairs', ('-',))[0], d.get('verify_fact_lce_cap', ('-',))[0], d.get('verify_fact_stack_peak', ('-',))[0]])
    md.append('### Generous runs: time, memory, I/O and disk (steps a, b, c; G4 reported, G2 and G3 judged here)\n')
    md.append(table(['case', 'run', 'finished', 'total s', 'process s'] + [p for _, p in PHASES] + ['ns/sym', 'RssAnon prog MB', 'RssAnon sampler MB', 'I/O counted B/char', 'rchar+wchar B/char', 'device read B/char', 'scratch peak B/char (prog)', 'scratch peak B/char (sampler)',
                                                                                               'outputs B/char', 'peak disk GB', 'cold start'], gen_rows))
    md.append('')
    md.append('Phase seconds are the program\'s own; "process s" includes the runtime start. ESTIMATE from E7 for comparison: 250 ns per symbol (205-300), i.e. 4.2 / 9.2 / 12.5 min at 1 / 2.2 / 3 Gbp on the prototype machine (4 vCPU x86, virtio disk).\n')
    md.append('### Verifier (e8-verify) after the generous runs (G1)\n')
    md.append(table(['case', 'run', 'all', 'verify s', 'files', 'meta', 'text', 'input', 'permutation', 'order', 'lcp values', 'lcp width', 'stats', 'child', 'sha256', 'long pairs', 'LCE cap', 'stack peak'], verify_rows))
    md.append('')

    # ---- identity
    ident = []
    ident_bad = 0
    ident_n = 0
    thread_cells = defaultdict(list)
    for (case, run), d in sorted(done.items()):
        if case.endswith('.gen') and run == '1' and '.t' not in case:
            continue
        if 'identity_vs_generous' in d:
            val = d['identity_vs_generous'][0]
            ident.append([case, run, d.get('exit_code', ('?',))[0], val])
            if val == 'NO':
                ident_bad += 1
            if val in ('yes', 'NO'):
                ident_n += 1
            if '.gen.t' in case:
                thread_cells[case.split('.')[0]].append(val)
    md.append('### Identity of the outputs against the generous run of the same text (G1)\n')
    md.append(table(['case', 'run', 'exit code', 'seven files SHA-256 equal'], ident) if ident else 'No limited-memory or thread run finished yet.')
    md.append('')

    # ---- cache table
    cache_rows = []
    verdicts = {}
    for size in sorted({c.split('.')[0] for (c, _) in done if '.m' in c}, key=lambda s: label_n(s) or 0):
        n = label_n(size)
        base_runs = generous.get(size, {})
        base = base_runs.get('1')
        base_total = base['total_s'] if base and base['finished'] else None
        cells = defaultdict(list)
        for (case, run), d in done.items():
            if case.startswith(size + '.m'):
                cells[int(case.split('.m')[1])].append((run, d))
        best = None
        for mem in sorted(cells, reverse=True):
            rs = [(run, run_numbers(d, n)) for run, d in sorted(cells[mem])]
            fin = [r for _, r in rs if r['finished']]
            if fin:
                tot = statistics.median(r['total_s'] for r in fin)
                slow = tot / base_total if base_total else None
                eff = statistics.median(r['eff'] for r in fin if r['eff'] is not None) if any(r['eff'] is not None for r in fin) else None
                ph = [statistics.median(r[k] for r in fin if r[k] is not None) if any(r[k] is not None for r in fin) else None for k, _ in PHASES]
                status = 'finished (%d run%s)' % (len(fin), '' if len(fin) == 1 else 's')
                if slow is not None and slow <= SLOWDOWN and (best is None or mem < best[0]):
                    best = (mem, slow, eff)
                row = [size, mem, status, fmt(eff, 3), fmt(tot, 1), fmt(slow, 2) + 'x' if slow else '-'] + [fmt(x, 1) for x in ph] + [fmt(statistics.median(r['dev_read'] for r in fin if r['dev_read'] is not None), 2) if any(r['dev_read'] is not None for r in fin) else '-',
                                                                                                                                 fmt(max(r['failcnt'] or 0 for r in fin), 0), fmt(max((r['rss_ext'] or 0) for r in fin), 1), 'yes' if all(r['eff'] is not None for r in fin) else '-']
            else:
                r = rs[0][1]
                status = 'does not finish in practice at this memory (%s)' % r['killed'] if r['killed'] in ('watchdog',) else 'did not finish (%s, exit code %s)' % (r['killed'], '137' if r['killed'] == 'oom-or-signal' else '?')
                row = [size, mem, status, fmt(r['eff'], 3), '-', '> %dx' % int(6) if r['killed'] == 'watchdog' else '-'] + ['-'] * len(PHASES) + [fmt(r['dev_read'], 2), fmt(r['failcnt'], 0), fmt(r['rss_ext'], 1), '-']
            cache_rows.append(row)
        verdicts[size] = best
    md.append('### The limited-cache sweep (step d; G4 reported): the container limit, the effective cache, the time against the generous run\n')
    md.append(table(['text', 'container memory MiB', 'result', 'effective cache of text + rank in PH3', 'total s', 'slowdown', 'PH0 s', 'PH1 s', 'PH2/3 s', 'PH4+5 s', 'device read B/char', 'cgroup failcnt / events', 'RssAnon sampler MB', 'residency sampled'], cache_rows) if cache_rows else 'No limited-memory cell finished or was abandoned yet.')
    md.append('')
    md.append("Slowdown = total_s of the cell (median of its runs) / total_s of run 1 of the generous run of the same text. 'Effective cache' = time-weighted resident fraction of the packed text and the DCS rank file over PH3 (mincore, 200 ms), as in E7. Abandon rule: " + GATES['abandon'] + '.\n')
    md.append('#### Verdict: the smallest memory with a slowdown of at most %.0fx (the E7 criterion)\n' % SLOWDOWN)
    vrows = []
    for size, best in verdicts.items():
        n = label_n(size)
        base = (generous.get(size, {}) or {}).get('1')
        text_rank = None
        for (case, run), d in done.items():
            if case.startswith(size + '.m'):
                tb, rb = num(d, 'text_bytes'), num(d, 'rank_bytes')
                if tb and rb:
                    text_rank = tb + rb
        if best:
            mem, slow, eff = best
            vrows.append([size, '%d MiB' % mem, fmt(slow, 2) + 'x', fmt(eff, 3), '%.2f x (text + rank file = %s MiB)' % (mem * 1048576 / text_rank, fmt(text_rank / 1048576, 0)) if text_rank else '-', '%.2f B/char' % (mem * 1048576 / n)])
        else:
            vrows.append([size, 'none of the finished cells', '-', '-', '-', '-'])
    md.append(table(['text', 'smallest container memory with slowdown <= 3x', 'slowdown', 'effective cache', 'as a multiple of text + rank file', 'per symbol'], vrows) if vrows else 'No cache cells.')
    md.append('')
    md.append('The open question of ADR-004 6c (3): how much page cache PH4\'s 12 B/char of scratch needs beyond the 1x (text + rank file). Read it from the table above: the container limit includes the anonymous memory (about 0.1-0.15 GB), the page cache of the text, the rank file and the scratch pages in flight; the first cell whose PH4+5 time grows against the generous run is the one where the scratch no longer fits next to the text. The interpretation is for the orchestrator; the numbers are above.\n')

    # ---- gates
    # G1
    ver_state = {}
    for size, label in (('1000m', '1 Gbp'), ('2200m', '2.2 Gbp'), ('3000m', '3 Gbp')):
        runs = generous.get('n' + size, {})
        if not runs:
            ver_state[size] = 'not run'
            continue
        vs = [done[('n%s.gen' % size, r)].get('verify_all', ('not run',))[0] for r in runs]
        ver_state[size] = 'pass' if all(v == 'pass' for v in vs) else ('FAIL' if any(v == 'FAIL' for v in vs) else 'not run')
    if any(v == 'FAIL' for v in ver_state.values()) or ident_bad:
        g1 = 'FAIL'
    elif all(v == 'pass' for v in ver_state.values()) and ident_n and thread_cells:
        g1 = 'PASS'
    else:
        g1 = 'NOT RUN (partly)' if any(v == 'pass' for v in ver_state.values()) or ident_n else 'NOT RUN'
    gates['G1'] = {'state': g1, 'verifier': ver_state, 'identity_compared': ident_n, 'identity_differences': ident_bad, 'thread_comparisons': {k: v for k, v in thread_cells.items()}}
    # G2
    rss = {}
    for (case, run), d in done.items():
        v = num(d, 'ext_rssanon_peak')
        if v is not None:
            rss[case.split('.')[0]] = max(rss.get(case.split('.')[0], 0), v)
    g2 = 'NOT RUN' if not rss else ('PASS' if all(v <= LIMIT_RSS_MB for v in rss.values()) else 'FAIL')
    gates['G2'] = {'state': g2, 'max_rssanon_mb_by_size': rss}
    # G3
    c3 = generous.get('n3000m', {}).get('1')
    if c3 and c3['finished']:
        sp = [x for x in (c3['scratch_ext'], c3['scratch_prog']) if x is not None]
        dk = c3['disk_peak']
        g3 = 'PASS' if sp and max(sp) <= LIMIT_SCRATCH and dk is not None and dk <= LIMIT_DISK_BYTES else ('NOT RUN (a measurement is missing)' if not sp or dk is None else 'FAIL')
        gates['G3'] = {'state': g3, 'scratch_peak_b_per_char': sp, 'peak_total_disk_gb': (dk or 0) / 1e9}
    else:
        gates['G3'] = {'state': 'NOT RUN'}
    gates['G4'] = {'state': 'reported', 'verdict_by_size': {k: (list(v) if v else None) for k, v in verdicts.items()}}
    md.append('### Gates\n')
    grows = [
        ['G1', GATES['G1'], gates['G1']['state'], 'verifier by size: ' + ', '.join('%s %s' % kv for kv in ver_state.items()) + '; identity comparisons %d, differences %d; thread comparisons: %s' % (ident_n, ident_bad, ', '.join('%s %s' % (k, v) for k, v in thread_cells.items()) or 'none')],
        ['G2', GATES['G2'], gates['G2']['state'], ', '.join('%s %.1f MB' % (k, v) for k, v in sorted(rss.items())) or '-'],
        ['G3', GATES['G3'], gates['G3']['state'], ('scratch peak %s B/char, peak total disk %.1f GB (3 Gbp generous run)' % ('/'.join(fmt(x, 2) for x in gates['G3'].get('scratch_peak_b_per_char', [])), gates['G3'].get('peak_total_disk_gb', 0))) if 'scratch_peak_b_per_char' in gates['G3'] else 'no finished 3 Gbp generous run'],
        ['G4', GATES['G4'], 'reported', 'see the tables above'],
    ]
    md.append(table(['gate', 'text (verbatim)', 'result', 'numbers'], grows))
    md.append('')
    abandoned = [(c, r) for (c, r), d in done.items() if d.get('killed', ('no',))[0] == 'watchdog']
    md.append('Abandon rule: ' + GATES['abandon'] + '. Cells stopped by it: ' + (', '.join('%s run %s' % cr for cr in sorted(abandoned)) or 'none') + '.\n')
    pending = sorted({k for k in data if k not in done and k[0] not in ('preflight',) and not k[0].startswith('gen.')})
    if pending:
        md.append('Runs with rows but no run_done row (interrupted; the runner purges and repeats them): ' + ', '.join('%s run %s' % k for k in pending) + '\n')
    text = '\n'.join(md)
    print(text)
    if out_md:
        open(out_md, 'w').write(text + '\n')
    if out_json:
        json.dump(gates, open(out_json, 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
