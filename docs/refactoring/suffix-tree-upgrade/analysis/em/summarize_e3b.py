#!/usr/bin/env python3
"""Summarise the raw rows of the E3b protocol: median [min-max] over the runs of every (case, metric); rows are
R <TAB> E3b <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes.
usage: summarize_e3b.py rows.tsv [rows2.tsv ...] [--metrics a,b,c]"""
import sys, statistics, collections

def main(argv):
    metrics = None
    files = []
    i = 0
    while i < len(argv):
        if argv[i] == "--metrics":
            metrics = argv[i + 1].split(","); i += 2
        else:
            files.append(argv[i]); i += 1
    if metrics is None:
        metrics = ["ph1_s", "ns_per_symbol", "rssanon_peak", "rssfile_peak", "arena_peak", "io_total_counted", "scratch_peak", "levels", "nonunique_share", "rank_equals_oracle"]
    data = collections.OrderedDict()
    for f in files:
        for line in open(f):
            c = line.rstrip("\n").split("\t")
            if len(c) < 7 or c[0] != "R" or c[1] != "E3b":
                continue
            case, run, metric, value = c[2], c[3], c[4], c[5]
            data.setdefault(case, collections.OrderedDict()).setdefault(metric, []).append(value)
    for case, m in data.items():
        out = []
        for name in metrics:
            if name not in m:
                continue
            vals = m[name]
            try:
                nums = [float(v) for v in vals]
                med = statistics.median(nums)
                out.append(f"{name} {med:.4g} [{min(nums):.4g}-{max(nums):.4g}]" if len(nums) > 1 else f"{name} {med:.4g}")
            except ValueError:
                out.append(f"{name} {'/'.join(vals)}")
        print(f"{case:34s} runs={len(m.get('ph1_s', []))}  " + "; ".join(out))

if __name__ == "__main__":
    main(sys.argv[1:])
