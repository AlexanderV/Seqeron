#!/usr/bin/env python3
"""The tables of ADR-005 section 5 from results/*.tsv (markdown on stdout): python3 adr_tables.py synthetic|descent <label> [group]"""
import csv, sys, collections
mode = sys.argv[1]
def load(p): return list(csv.DictReader(open(p), delimiter='\t'))
if mode == 'synthetic':
    rows = [r for r in load('results/width_synthetic.tsv') if r['kind'] in ('chase', 'thr', 'scan')]
    layouts = ['w4', 'w5', 'w8', 'hyb', 'p34']
    for metric in ['chase_ns_per_read', 'thr_ns_per_read', 'scan1500_ns_per_entry']:
        for thr in ['1', '4']:
            print(f"\n{metric}, {thr} thread(s): median [min-max] ns (ratio to w4)")
            print("| n | " + " | ".join(layouts) + " |")
            print("|---|" + "---|" * len(layouts))
            for n in sorted({int(r['n']) for r in rows}):
                cells = {}
                for r in rows:
                    if r['metric'] == metric and r['threads'] == thr and int(r['n']) == n:
                        cells[r['layout']] = r
                if 'w4' not in cells: continue
                b = float(cells['w4']['median'])
                out = []
                for l in layouts:
                    r = cells.get(l)
                    out.append(f"{float(r['median']):.1f} [{float(r['min']):.1f}-{float(r['max']):.1f}] ({float(r['median'])/b:.2f})" if r else '-')
                print(f"| {n:,} | " + " | ".join(out) + " |")
elif mode == 'descent':
    label = sys.argv[2]
    group = sys.argv[3] if len(sys.argv) > 3 else None
    rows = [r for r in load('results/width_descent.tsv') if label in r['notes'] and (group is None or f'group {group}' in r['notes'])]
    ref = 'v1'
    for metric in ['count20_present', 'count20_absent', 'locate8']:
        for thr in ['4', '1']:
            rs = [r for r in rows if r['metric'].startswith(metric) and r['threads'] == thr]
            if not rs: continue
            base = float(next(r for r in rs if r['layout'] == ref)['median']) if any(r['layout'] == ref for r in rs) else None
            print(f"\n{metric}, {thr} thread(s), {label}{' ' + group if group else ''}: ns per query per thread, median [min-max] (ratio to {ref})")
            for r in rs:
                m = float(r['median'])
                print(f"  {r['layout']:5} {m:8.1f} [{float(r['min']):8.1f}-{float(r['max']):8.1f}]" + (f" ({m/base:.3f})" if base else ''))
