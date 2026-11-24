#!/usr/bin/env python3
"""Descent rows as a table with the ratio to a reference variant: python3 summarize_descent.py <tsv> [label-substring] [reference-variant]"""
import sys, csv, collections
path = sys.argv[1]
lab = sys.argv[2] if len(sys.argv) > 2 else ''
ref = sys.argv[3] if len(sys.argv) > 3 else 'v1'
rows = [r for r in csv.DictReader(open(path), delimiter='\t') if lab in r['notes']]
cells = collections.OrderedDict()
for r in rows:
    cells.setdefault((r['n'], r['metric'], r['threads']), {})[r['layout']] = r
for (n, metric, thr), vs in cells.items():
    base = float(vs[ref]['median']) if ref in vs else None
    print(f"\n== n={int(n):,} {metric} threads={thr} (ns per query per thread, median [min-max]; ratio to {ref})")
    for v, r in vs.items():
        m, lo, hi = float(r['median']), float(r['min']), float(r['max'])
        print(f"  {v:4} {m:9.1f} [{lo:8.1f}-{hi:8.1f}]  x{(m / base if base else float('nan')):.3f}   {r['notes'].split('bytes(SA+child)=')[1][:90] if 'bytes(SA+child)=' in r['notes'] else ''}")
