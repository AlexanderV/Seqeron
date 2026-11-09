#!/usr/bin/env python3
"""Summarise raw rows of the ADR-004 prototypes (E1/E2/E3): median [min-max] over the runs of every (experiment, case, metric).

Usage: python3 -I summarize.py <rows.tsv>... [--metrics sa_s,tie_share,...] [--ref <case>]
Rows: R <TAB> experiment <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes (the run column is optional).
Prints one markdown table row per case; --ref names the case whose sa_s median is the denominator of the ratio column.
"""
import statistics
import sys


def load(paths):
    rows = []
    for p in paths:
        with open(p, encoding="utf-8") as f:
            for line in f:
                if not line.startswith("R\t"):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) == 7:  # no run column
                    c = c[:3] + ["1"] + c[3:]
                rows.append(c[1:8])  # exp, case, run, metric, value, unit, notes
    return rows


def main():
    args = sys.argv[1:]
    metrics = ["sa_s", "tie_share", "lce_symbols", "compares"]
    ref = None
    paths = []
    i = 0
    while i < len(args):
        if args[i] == "--metrics":
            metrics = args[i + 1].split(",")
            i += 2
        elif args[i] == "--ref":
            ref = args[i + 1]
            i += 2
        else:
            paths.append(args[i])
            i += 1
    rows = load(paths)
    cases = []
    data = {}
    status = {}
    for exp, case, run, metric, value, unit, notes in rows:
        key = (exp, case)
        if key not in cases:
            cases.append(key)
        if metric == "status":
            status.setdefault(key, []).append(value)
            continue
        try:
            v = float(value)
        except ValueError:
            continue
        data.setdefault((key, metric), []).append(v)

    def fmt(v):
        return f"{v:.4g}" if abs(v) < 1e5 else f"{v:.3e}"

    refmed = None
    if ref:
        for key in cases:
            if key[1] == ref and (key, "sa_s") in data:
                refmed = statistics.median(data[(key, "sa_s")])
    print("| exp | case | runs | " + " | ".join(metrics) + (" | ratio to " + ref if refmed else "") + " |")
    print("|---|---|---|" + "---|" * (len(metrics) + (1 if refmed else 0)))
    for key in cases:
        cells = []
        runs = 0
        for m in metrics:
            vs = data.get((key, m))
            if vs:
                runs = max(runs, len(vs))
                cells.append(f"{fmt(statistics.median(vs))} [{fmt(min(vs))}-{fmt(max(vs))}]" if len(vs) > 1 else fmt(vs[0]))
            else:
                cells.append("-")
        extra = ""
        if refmed:
            vs = data.get((key, "sa_s"))
            extra = f" | {fmt(statistics.median(vs) / refmed)}x" if vs else " | -"
        st = status.get(key)
        label = f"{runs}" if not st else f"{runs} ({', '.join(sorted(set(st)))})"
        print(f"| {key[0]} | {key[1]} | {label} | " + " | ".join(cells) + extra + " |")


if __name__ == "__main__":
    main()
