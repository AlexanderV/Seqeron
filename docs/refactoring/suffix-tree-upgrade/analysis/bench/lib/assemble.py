#!/usr/bin/env python3
"""Merge every runs/raw*.tsv into bench.tsv (median over repeats) + append the NA rows of runs/na_rows.tsv.
usage: assemble.py <bench-dir>"""
import glob
import os
import subprocess
import sys

B = sys.argv[1]
subprocess.check_call([sys.executable, f"{B}/lib/make_na.py", B])
raws = sorted(glob.glob(f"{B}/runs/raw*.tsv"))
raws = [r for r in raws if not r.endswith("raw_all.tsv")]
all_raw = f"{B}/runs/raw_all.tsv"
with open(all_raw, "w") as o:
    first = True
    for r in raws:
        with open(r) as f:
            h = next(f)
            if first:
                o.write(h)
                first = False
            for line in f:
                o.write(line)
out = f"{B}/bench.tsv"
subprocess.check_call([sys.executable, f"{B}/lib/benchlib.py", "aggregate", "--raw", all_raw, "--out", out])
na = f"{B}/runs/na_rows.tsv"
if os.path.exists(na):
    with open(out, "a") as o, open(na) as f:
        for line in f:
            if line.strip() and not line.startswith("tool\t"):
                o.write(line)
print("bench.tsv rows:", sum(1 for _ in open(out)) - 1)
