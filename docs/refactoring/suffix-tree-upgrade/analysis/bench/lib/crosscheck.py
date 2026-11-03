#!/usr/bin/env python3
"""Cross-tool oracle check: every subject must return identical count/locate checksums on the same input.
usage: crosscheck.py <bench.tsv>   -> prints a markdown table and exits 1 on any disagreement"""
import csv
import sys
from collections import defaultdict

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
KEYS = [("count_present20", "checksum"), ("count_absent20", "checksum"), ("contains_present20", "checksum"),
        ("locate8", "total_hits"), ("locate8", "checksum"), ("locate4", "checksum")]
by = defaultdict(lambda: defaultdict(dict))
for r in rows:
    if (r["operation"], r["metric"]) in KEYS:
        by[r["input"]][(r["operation"], r["metric"])][r["tool"]] = r["value"]
bad = 0
print("| input | check | agreeing tools | value | DISAGREE |\n|---|---|---|---|---|")
for inp in sorted(by):
    for k in KEYS:
        d = by[inp].get(k)
        if not d:
            continue
        vals = defaultdict(list)
        for t, v in d.items():
            vals[v].append(t)
        dis = "" if len(vals) == 1 else "; ".join(f"{v}: {','.join(ts)}" for v, ts in vals.items())
        if dis:
            bad += 1
        main = max(vals.items(), key=lambda kv: len(kv[1]))
        print(f"| {inp} | {k[0]}/{k[1]} | {len(main[1])}/{len(d)} | {main[0]} | {dis} |")
sys.exit(1 if bad else 0)
