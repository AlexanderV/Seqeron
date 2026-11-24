#!/usr/bin/env python3
"""Summarise results/width_synthetic.tsv (and descent/text rows) as tables: python3 summarize_width.py <tsv> [metric-substring]"""
import sys, csv, collections
path = sys.argv[1]
want = sys.argv[2] if len(sys.argv) > 2 else None
rows = list(csv.DictReader(open(path), delimiter='\t'))
cells = collections.OrderedDict()
for r in rows:
    if want and want not in r['metric']:
        continue
    cells.setdefault((r['kind'], r['metric'], r['threads']), []).append(r)
for (kind, metric, thr), rs in cells.items():
    ns = sorted({int(r['n']) for r in rs})
    layouts = list(collections.OrderedDict.fromkeys(r['layout'] for r in rs))
    print(f"\n== {kind} {metric} threads={thr}  (median [min-max] {rs[0]['unit']})")
    print("layout".ljust(8) + "".join(f"n={n:,}".rjust(26) for n in ns))
    for l in layouts:
        line = l.ljust(8)
        for n in ns:
            m = [r for r in rs if r['layout'] == l and int(r['n']) == n]
            line += (f"{float(m[0]['median']):8.2f} [{float(m[0]['min']):.2f}-{float(m[0]['max']):.2f}]".rjust(26) if m else "-".rjust(26))
        print(line)
