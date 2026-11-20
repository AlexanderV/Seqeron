#!/usr/bin/env python3
"""Work package 1b.9a: medians over the processes of a rows file written by run_1b9a_matrix.sh.
For each variant: the median, the minimum and the maximum of count20 / locate8 / absent24 per second over its processes, and the slowdown against the reference variant
(the median of the throughputs, and the median of the per-round paired ratios: round r of the variant against round r of the reference).
usage: summarize_1b9a.py <rows file> [reference variant, default memory] [more rows files]"""
import statistics
import sys
from collections import defaultdict

def load(path):
    rows = []
    for line in open(path, encoding="utf-8", errors="replace"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "row" and len(f) >= 7:
            rows.append(dict(round=f[1], variant=f[2], c=float(f[3]), l=float(f[4]), a=float(f[5]), open=float(f[6]), extra=f[8:]))
    return rows

def main():
    paths = [p for p in sys.argv[1:] if not p.startswith("--")]
    ref = "memory"
    if len(paths) > 1 and not paths[1].endswith(".tsv"):
        ref = paths[1]
        paths = [paths[0]] + paths[2:]
    rows = [r for p in paths for r in load(p)]
    by = defaultdict(list)
    for r in rows:
        by[r["variant"]].append(r)
    rounds = defaultdict(dict)
    for r in rows:
        rounds[r["variant"]][r["round"]] = r
    print(f"{'variant':34s} {'n':>3s} | {'count20/s':>10s} {'min':>8s} {'max':>8s} {'slower%':>8s} {'paired%':>8s} | {'locate8/s':>10s} {'slower%':>8s} {'paired%':>8s} | {'absent24/s':>10s} {'slower%':>8s} {'paired%':>8s}")
    refrows = by.get(ref, [])
    def med(rs, k):
        return statistics.median(x[k] for x in rs)
    for v in sorted(by, key=lambda v: (v != ref, v)):
        rs = by[v]
        out = f"{v:34s} {len(rs):3d} |"
        for k in ("c", "l", "a"):
            m = med(rs, k)
            slow = (1 - m / med(refrows, k)) * 100 if refrows else float("nan")
            pr = []
            for rd, x in rounds[v].items():
                y = rounds[ref].get(rd)
                if y:
                    pr.append((1 - x[k] / y[k]) * 100)
            pm = statistics.median(pr) if pr else float("nan")
            if k == "c":
                out += f" {m:10.0f} {min(x[k] for x in rs):8.0f} {max(x[k] for x in rs):8.0f} {slow:8.1f} {pm:8.1f} |"
            else:
                out += f" {m:10.0f} {slow:8.1f} {pm:8.1f} |"
        print(out.rstrip("|"))
    print("\nslower% = (1 - median(variant) / median(reference)) * 100 over all processes; paired% = the median over rounds of the same quantity for the processes of one round")

main()
