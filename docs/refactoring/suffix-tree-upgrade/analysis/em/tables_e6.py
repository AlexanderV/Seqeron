#!/usr/bin/env python3
"""Tables of the raw rows of E6: python3 -I tables_e6.py <rows tsv> phases <case regex> | totals <case regex>... | cmp
Rows are R <TAB> E6 <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes (rows without a numeric run column are skipped).
  phases <regex>   the per-phase table (median [min-max] over the runs of every case that matches): wall time, program I/O (explicit read, explicit write, mapped write), /proc/<pid>/io of the phase (rchar + wchar), RssAnon peak, scratch peak so far, dirty pages
  totals <regex>   the gate numbers per case: total_s, sa_stage_s, io_total_counted per char, rssanon (program, external), dirty (program, external, per-process mapped), scratch
  cmp              every cmp_* and *_equals_* row: how many yes and NO per case"""
import re, statistics, sys
from collections import defaultdict

path, mode = sys.argv[1], sys.argv[2]
pattern = re.compile(sys.argv[3]) if len(sys.argv) > 3 else re.compile('.')
data = defaultdict(lambda: defaultdict(list))
order = []
n_of = {}
for line in open(path):
    c = line.rstrip('\n').split('\t')
    if len(c) < 7 or c[0] != 'R' or c[1] != 'E6' or not c[3].isdigit():
        continue
    case, run, metric, value = c[2], c[3], c[4], c[5]
    if not pattern.search(case):
        continue
    if case not in order:
        order.append(case)
    data[case][metric].append(value)
    if metric == 'n':
        n_of[case] = float(value)

def nums(case, metric):
    out = []
    for x in data[case].get(metric, []):
        try:
            out.append(float(x))
        except ValueError:
            pass
    return out

def fmt(f, digits=4):
    if not f:
        return '-'
    if len(f) == 1:
        return ('%.' + str(digits) + 'g') % f[0]
    return ('%.' + str(digits) + 'g [%.' + str(digits) + 'g-%.' + str(digits) + 'g]') % (statistics.median(f), min(f), max(f))

def per_char(case, f):
    n = n_of.get(case, 1.0)
    return [x / n for x in f]

if mode == 'phases':
    phases = ['ph0_prepare', 'ph1_dcs', 'ph23_sa', 'ph45_lcp_child', 'sync']
    for case in order:
        print('\n### %s (%d runs)\n' % (case, len(data[case].get('total_s', []))))
        print('| phase | wall s | explicit read B/char | explicit write B/char | mapped write B/char | rchar + wchar B/char | RssAnon peak MB | scratch peak so far B/char | dirty peak MB (system / mapped outputs) |')
        print('|---|---|---|---|---|---|---|---|---|')
        for ph in phases:
            if not nums(case, ph + '_s'):
                continue
            rw = [a + b for a, b in zip(nums(case, ph + '_proc_rchar'), nums(case, ph + '_proc_wchar'))]
            print('| %s | %s | %s | %s | %s | %s | %s | %s | %s / %s |' % (
                ph, fmt(nums(case, ph + '_s')), fmt(per_char(case, nums(case, ph + '_io_explicit_read')), 3), fmt(per_char(case, nums(case, ph + '_io_explicit_write')), 3),
                fmt(per_char(case, nums(case, ph + '_io_mapped_write')), 3), fmt(per_char(case, rw), 3), fmt(nums(case, ph + '_rssanon_peak')), fmt(per_char(case, nums(case, ph + '_scratch_peak_cum')), 3),
                fmt(nums(case, ph + '_dirty_peak'), 3), fmt(nums(case, ph + '_proc_dirty_peak'), 3)))
elif mode == 'totals':
    print('| case | total s | SA stage s (PH0-PH3) | I/O B/char (program) | rchar + wchar B/char (program run / whole process) | RssAnon MB (program / independent) | dirty MB (program system / independent / mapped outputs) | scratch peak B/char (tracked / apparent / allocated) | footprint B/char |')
    print('|---|---|---|---|---|---|---|---|---|')
    for case in order:
        n = n_of.get(case, 1.0)
        ext_rw = [a + b for a, b in zip(nums(case, 'ext_proc_rchar'), nums(case, 'ext_proc_wchar'))]
        print('| %s | %s | %s | %s | %s / %s | %s / %s | %s / %s / %s | %s / %s / %s | %s |' % (
            case, fmt(nums(case, 'total_s')), fmt(nums(case, 'sa_stage_s')), fmt(per_char(case, nums(case, 'io_total_counted')), 3), fmt(per_char(case, nums(case, 'proc_rw_total')), 3), fmt(per_char(case, ext_rw), 3),
            fmt(nums(case, 'rssanon_peak')), fmt(nums(case, 'ext_rssanon_peak')), fmt(nums(case, 'dirty_peak'), 3), fmt(nums(case, 'ext_dirty_peak'), 3), fmt(nums(case, 'proc_dirty_peak'), 3),
            fmt(per_char(case, nums(case, 'scratch_peak')), 3), fmt(per_char(case, nums(case, 'ext_scratch_apparent_peak')), 3), fmt(per_char(case, nums(case, 'ext_scratch_allocated_peak')), 3), fmt(per_char(case, nums(case, 'output_bytes')), 3)))
elif mode == 'cmp':
    for case in order:
        yes = no = 0
        names = []
        for metric, vals in data[case].items():
            if metric.startswith('cmp_') or metric.endswith('_equals_ref') or metric.endswith('_equals_libsais') or metric == 'exit_code' or metric == 'clean' or metric == 'scratch_left':
                for v in vals:
                    if v in ('yes', '0'):
                        yes += 1
                    else:
                        no += 1
                        names.append('%s=%s' % (metric, v))
        print('%-36s yes=%d NO=%d %s' % (case, yes, no, ' '.join(names)))
