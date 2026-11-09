#!/usr/bin/env python3
"""E3: for every cache-limited run print the sub-group limit, the mean anonymous memory of the sub-group during the run, the effective page-cache budget
(limit - mean rss) as a multiple of (text + rank), and the time. Usage: python3 -I e3_effective_cache.py <rows.tsv> <text bytes + rank bytes> [baseline_case]
Limits are recomputed as the protocol script does: (85 + 8) MiB + rho * (text + rank), rho taken from the case name 'cache_rho<rho>'."""
import re
import statistics
import sys

rows = {}
for line in open(sys.argv[1], encoding="utf-8"):
    c = line.rstrip("\n").split("\t")
    if c[0] != "R" or len(c) != 8:
        continue
    _, _, case, run, metric, value, unit, notes = c
    rows.setdefault((case, run), {})[metric] = (value, notes)
t = float(sys.argv[2])
MIB = 1048576
out = []
for (case, run), m in sorted(rows.items()):
    mo = re.search(r"cache_rho([0-9.]+)", case)
    if not mo or "sa_s" not in m:
        continue
    rho = float(mo.group(1))
    limit = (85 + 8) * MIB + rho * t
    mean_rss = float(re.search(r"mean ([0-9.]+) MB", m["cgroup_rss"][1]).group(1)) * 1e6 if "cgroup_rss" in m else float("nan")
    eff = (limit - mean_rss) / t
    out.append((rho, case, run, limit / 1e6, mean_rss / 1e6, eff, float(m["sa_s"][0])))
print("| rho (label) | run | limit (MB) | mean rss (MB) | effective cache / (text+rank) | sa_s (s) |")
print("|---|---|---|---|---|---|")
for rho, case, run, lim, rss, eff, s in sorted(out, key=lambda x: (-x[0], x[2])):
    print(f"| {rho} | {run} | {lim:.0f} | {rss:.0f} | {eff:.2f} | {s:.1f} |")
