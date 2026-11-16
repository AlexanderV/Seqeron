#!/usr/bin/env python3
# Tables of the E9 timed runs (e9_protocol.tsv, rows  R E9 label run metric value unit notes): per label the number of runs, the median [min-max] of the time, RssAnon (program sampler and independent sampler), and the
# I/O in bytes per character by the three counters the gate names: the program's counters (explicit reads + explicit writes + mapped writes), /proc/<pid>/io rchar + wchar (+ the mapped writes the syscalls do not see),
# and the block layer (/sys/block/vda/stat sectors, at exit and after sync, with the cost of starting the runtime from the e9.base runs subtracted).
# usage: python3 -I tables_e9.py <protocol rows file> [--identity <identity rows file>]
import statistics, sys
from collections import defaultdict

N = {'dna_10m': 10_000_000, 'prot_1m': 1_000_000, 'rep_10m': 10_000_000, 'a_1m': 1_000_000, 'fib_1m': 1_000_000, 'tm_1m': 1_000_000, 'dna_100m': 100_000_000, 'rnd256_10m': 10_000_000, 'tiny': 1000}
PH = {'psi': 'psi', 'psi256': 'psi', 'bwt': 'bwt', 'isa': 'isa', 'psiisa': 'psi_from_isa', 'base': 'psi'}

rows = defaultdict(lambda: defaultdict(dict))   # label -> run -> metric -> value
for line in open(sys.argv[1]):
    f = line.rstrip('\n').split('\t')
    if len(f) < 6 or f[0] != 'R' or f[1] != 'E9' or f[2] in ('protocol',):
        continue
    label, run, metric, value = f[2], f[3], f[4], f[5]
    if run in ('head',):
        continue
    rows[label][run][metric] = value

def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

def med(label, getter):
    xs = []
    for run, m in rows[label].items():
        if m.get('run_done') != 'yes':
            continue
        v = getter(m)
        if v is not None:
            xs.append(v)
    if not xs:
        return None
    return statistics.median(xs), min(xs), max(xs), len(xs)

def fmt(t, digits=3):
    if t is None:
        return '-'
    m, lo, hi, n = t
    return ('%.' + str(digits) + 'g [%.' + str(digits) + 'g-%.' + str(digits) + 'g]') % (m, lo, hi)

# the start-up cost of the runtime (e9.base): block layer bytes read
base = {}
for key in ('blk_read_at_exit', 'blk_write_after_sync', 'ext_proc_rchar', 'ext_proc_wchar'):
    t = med('e9.base' if 'e9.base' in rows else 'tiny.e9.base', lambda m, key=key: num(m.get(key)))
    base[key] = t[0] if t else 0.0
print('# start-up cost of the runtime, e9.base (bytes): ' + ', '.join('%s=%.0f' % kv for kv in base.items()))
print()
print('| label | runs | phase s (program) | process wall s | RssAnon MB program / independent | counted B/char (explicit r + w + mapped) | rchar+wchar B/char (phase / whole process) | block read B/char at exit (net of start-up) | block write B/char after sync | cmp |')
print('|---|---|---|---|---|---|---|---|---|---|')
for label in sorted(rows):
    name = label.split('.')[0]
    grp = label.split('.')[-1]
    ph = PH.get(grp, 'psi')
    n = N.get(name, N.get('tiny'))
    t_s = med(label, lambda m: num(m.get(ph + '_s')))
    t_w = med(label, lambda m: num(m.get('process_wall_s')))
    r_p = med(label, lambda m: num(m.get(ph + '_rssanon_peak')))
    r_e = med(label, lambda m: num(m.get('ext_rssanon_peak')))
    cnt = med(label, lambda m: num(m.get(ph + '_io_counted_total')) / n if num(m.get(ph + '_io_counted_total')) is not None else None)
    rw_phase = med(label, lambda m: (num(m.get(ph + '_proc_rchar')) + num(m.get(ph + '_proc_wchar'))) / n if num(m.get(ph + '_proc_rchar')) is not None else None)
    rw_proc = med(label, lambda m: (num(m.get('ext_proc_rchar')) + num(m.get('ext_proc_wchar'))) / n if num(m.get('ext_proc_rchar')) is not None else None)
    br = med(label, lambda m: (num(m.get('blk_read_at_exit')) - base['blk_read_at_exit']) / n if num(m.get('blk_read_at_exit')) is not None else None)
    bw = med(label, lambda m: (num(m.get('blk_write_after_sync')) - base['blk_write_after_sync']) / n if num(m.get('blk_write_after_sync')) is not None else None)
    cmps = sorted({m.get('cmp_output', '-') for m in rows[label].values() if m.get('run_done') == 'yes'})
    r_cell = '%s / %s' % (fmt(r_p, 4), fmt(r_e, 4))
    rw_cell = '%s / %s' % (fmt(rw_phase, 3), fmt(rw_proc, 3))
    print('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (label, t_s[3] if t_s else 0, fmt(t_s), fmt(t_w), r_cell, fmt(cnt, 4), rw_cell, fmt(br, 3), fmt(bw, 3), ','.join(cmps)))
