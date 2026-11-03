#!/usr/bin/env python3
"""Print markdown comparison tables from bench.tsv.  usage: summarize.py <bench.tsv> [input ...]"""
import csv
import sys
from collections import defaultdict

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
inputs = sys.argv[2:] or sorted({r["input"] for r in rows})
idx = defaultdict(dict)
for r in rows:
    idx[(r["tool"], r["input"])][(r["operation"], r["metric"])] = r["value"]


def g(tool, inp, op, met, fmt="{:,.3g}", scale=1.0):
    v = idx[(tool, inp)].get((op, met))
    if v is None and met == "peak_rss_hwm_after_build":
        v = idx[(tool, inp)].get(("build_process", "ext_peak_rss"))
    if v is None and op == "build" and met == "time":
        v = idx[(tool, inp)].get(("index_process", "ext_wall"))
    if v is None and met == "peak_rss_hwm_after_build":
        v = idx[(tool, inp)].get(("index_process", "ext_peak_rss"))
    if v is None and op == "build" and met == "bytes_per_char":
        v = idx[(tool, inp)].get(("index", "bytes_per_char"))
    if v is None:
        return "-"
    try:
        return fmt.format(float(v) * scale)
    except ValueError:
        return v


COLS = [
    ("build s", "build", "time", "{:,.3g}", 1),
    ("B/char idx", "build", "bytes_per_char", "{:,.3g}", 1),
    ("peak RSS MB (build)", "build", "peak_rss_hwm_after_build", "{:,.0f}", 1e-6),
    ("count q/s", "count_present20", "throughput", "{:,.0f}", 1),
    ("absent q/s", "count_absent20", "throughput", "{:,.0f}", 1),
    ("locate8 pat/s", "locate8", "throughput", "{:,.0f}", 1),
    ("count 4T q/s", "count_present20_4t", "throughput", "{:,.0f}", 1),
]
for inp in inputs:
    tools = sorted({t for (t, i) in idx if i == inp})
    if not tools:
        continue
    print(f"\n### {inp}\n")
    print("| tool | " + " | ".join(c[0] for c in COLS) + " |")
    print("|---|" + "---|" * len(COLS))
    for t in tools:
        print(f"| {t} | " + " | ".join(g(t, inp, op, m, f, s) for _, op, m, f, s in COLS) + " |")
