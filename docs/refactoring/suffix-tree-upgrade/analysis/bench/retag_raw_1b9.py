#!/usr/bin/env python3
"""Work package 1b.9: the first session wrote its `open` rows (warm file cache and dropped page cache) and the rows of the nav run with the in-memory Ukkonen tree under the same operation names; the runner
(run_1b9_gate.sh) now gives them the operation prefixes through lib/benchlib.py --op-prefix: `pwarm_` / `pcold_` for the first session's `open` rows (it read the pattern file before the first Open: the command `openp`; the
rows and their per-process rows were renamed from `warm_` / `cold_` (`open_*_process`) to these names once the later runs, with nothing before the first Open, took the shorter ones) and `ukk_`. This script gives the rows already written the same prefixes, from the phase of the
per-process rows that close every repeat (a repeat is the run of rows up to its `ext_cpu_user_sys`): `open_warmcache_process` -> `warm_`, `open_cold_process` -> `cold_`, `nav_*_ukk_process` -> `ukk_` on every
row of that repeat. Rows with a prefix already are left alone; nothing else is changed. usage: retag_raw_1b9.py <raw file> (in place)"""
import sys

PREFIX = {"popen_warmcache_process": "pwarm_", "popen_cold_process": "pcold_", "open_warmcache_process": "warm_", "open_cold_process": "cold_"}
KNOWN = ("pwarm_", "pcold_", "warm_", "cold_", "ukk_", "lean_")

path = sys.argv[1]
with open(path) as f:
    lines = f.read().split("\n")
out, block = [], []
header = lines[0]
changed = 0


def flush(block):
    global changed
    if not block:
        return
    proc = next((p[4] for p in block if p[4].endswith("_process")), None)
    prefix = PREFIX.get(proc, "ukk_" if proc and proc.startswith("nav_") and proc.endswith("_ukk_process") else "")
    for p in block:
        if prefix and not p[4].startswith(KNOWN) and not p[4].endswith("_process"):
            p[4] = prefix + p[4]
            changed += 1
        out.append(p)


for line in lines[1:]:
    if not line:
        continue
    p = line.split("\t")
    block.append(p)
    if p[5] == "ext_cpu_user_sys":
        flush(block)
        block = []
flush(block)
with open(path, "w") as f:
    f.write(header + "\n" + "\n".join("\t".join(p) for p in out) + "\n")
print("retagged", changed, "rows")
