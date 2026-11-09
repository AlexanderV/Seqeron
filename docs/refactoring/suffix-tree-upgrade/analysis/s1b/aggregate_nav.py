#!/usr/bin/env python3
"""Aggregation of nav_raw.tsv (work package 1b.4, run_nav.sh / run_nav_all.sh): per (input, op, group, metric) the median over the fresh processes with [min .. max] and the process count.
usage: aggregate_nav.py <nav_raw.tsv> [summary|rows]      summary (default): the headline table; rows: every metric."""
import re, statistics, sys
from collections import defaultdict

def parse(path):
    data = defaultdict(list)       # (input, op, group, metric) -> [values]
    notes = {}
    exits = []
    cur = None
    for line in open(path, encoding='utf-8'):
        line = line.rstrip('\n')
        if line.startswith('# process'):
            m = re.match(r'# process input=(\S+) op=(\S+) run=(\d+).*options: ?(.*)$', line)
            if m and 'EXIT' not in line:
                opts = m.group(4).strip()
                cur = (m.group(1), m.group(2) + (' ' + opts if opts else ''))
            elif 'EXIT' in line:
                exits.append(line)
            continue
        if line.startswith('R\t') and cur:
            f = line.split('\t')
            if len(f) < 5: continue
            try: v = float(f[3])
            except ValueError: continue
            data[(cur[0], cur[1], f[1], f[2])].append(v)
            if len(f) > 5: notes[(cur[0], cur[1], f[1], f[2])] = f[5]
    return data, notes, exits

def fmt(vals, digits=4):
    vals = sorted(vals)
    med = statistics.median(vals)
    return f"{med:.{digits}f} [{vals[0]:.{digits}f} .. {vals[-1]:.{digits}f}] n={len(vals)}"

def main():
    path = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else 'summary'
    data, notes, exits = parse(path)
    if mode == 'rows':
        for k in sorted(data):
            print('\t'.join(k), fmt(data[k], 6), sep='\t')
    else:
        inputs = sorted({k[0] for k in data})
        for inp in inputs:
            for op in ('mem', 'mum', 'anchors'):
                if (inp, op, 'nav_esa', f'{op}_first_call_s') not in data: continue
                print(f"== {inp} {op} (100 kbp query, minLength 20)")
                for group, label in (('nav_esa', 'ESA'), ('nav_ukkonen', 'Ukkonen')):
                    for metric in ('first_call_s', 'warm_median_s', 'steady_median_s'):
                        k = (inp, op, group, f'{op}_{metric}')
                        if k in data: print(f"  {label:8s} {metric:16s} {fmt(data[k])}")
                for metric in ('results',):
                    k = (inp, op, 'nav_esa', f'{op}_{metric}')
                    if k in data: print(f"  results {fmt(data[k], 0)}")
                eq = data.get((inp, op, 'nav_check', f'{op}_equal'))
                if eq: print(f"  answers equal to Ukkonen: {sum(eq):.0f} of {len(eq)}")
                for metric in ('build_esa_s', 'pack_build_s', 'pack_bytes_per_char', 'pyramid_build_s'):
                    k = (inp, op, 'nav', metric)
                    if k in data: print(f"  {metric:20s} {fmt(data[k])}")
                k = (inp, op, 'nav_ukkonen', 'build_s')
                if k in data: print(f"  ukkonen build_s      {fmt(data[k], 2)}")
                k = (inp, op, 'nav', 'peak_rss_hwm')
                if k in data: print(f"  peak_rss_hwm         {fmt([v / 1e9 for v in data[k]], 2)} GB")
            for op in sorted({k[1] for k in data if k[0] == inp and k[1].startswith('lcs')}):
                keys = sorted(k for k in data if k[0] == inp and k[1] == op and k[2] in ('nav_esa', 'nav_ukkonen'))
                if keys: print(f"== {inp} {op}")
                for k in keys: print(f"  {k[2]:12s} {k[3]:28s} {fmt(data[k])}")
            keys = sorted(k for k in data if k[0] == inp and k[1] == 'links' and k[2] == 'nav_links' and k[3] in ('production_median_ns', 'prototype_c2_isa_f32_pos_median_ns', 'ratio_production_over_prototype', 'differences'))
            if keys: print(f"== {inp} links")
            for k in keys: print(f"  {k[3]:40s} {fmt(data[k], 3)}")
            keys = sorted(k for k in data if k[0] == inp and k[1] == 'split')
            if keys: print(f"== {inp} split")
            for k in keys:
                if k[2] in ('split_calls', 'split_time', 'replay', 'split_plain') and not k[3].endswith(('checksum',)): print(f"  {k[2]:12s} {k[3]:28s} {fmt(data[k], 5)}")
    for e in exits: print(e)

main()
