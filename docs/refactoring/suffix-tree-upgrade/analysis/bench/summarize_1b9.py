#!/usr/bin/env python3
"""Work package 1b.9: median [min-max] over the fresh processes of a raw rows file of lib/benchlib.py (runs_esa/raw_1b9*.tsv).
usage: summarize_1b9.py <raw file> [input] [op-substring] [metric-substring]      one line per (input, op, metric): n, median, min, max, spread = max / min, unit
       summarize_1b9.py --abapp <rows file>...                                   the lines of abapp (rows files of run_ab.sh): per operation the median of the processes of A, B and B / A
       summarize_1b9.py --gate <rows file>                                       the table of the gate of 1b.9: per text and thread count the null experiment, file-backed (Checked) / in-memory in both orders and
                                                                                 balanced, Trusted / Checked, facade / index (the labels of run_1b9_session2.sh and run_1b9_session2b.sh)
"""
import os
import re
import statistics
import sys

VERSION = os.environ.get("VERSION_FILTER", "")   # keeps the rows whose version column contains this text (the libsais rows carry "session start|mid|end")


def raw(path, inp=None, op=None, metric=None):
    groups, order = {}, []
    with open(path) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            tool, ver, x, n, o, m, v, unit, rep = p[:9]
            if VERSION and VERSION not in ver:
                continue
            if inp and x != inp:
                continue
            if op and op not in o:
                continue
            if metric and metric not in m:
                continue
            try:
                val = float(v)
            except ValueError:
                continue
            k = (x, o, m, unit)
            if k not in groups:
                groups[k] = []
                order.append(k)
            groups[k].append(val)
    return groups, order


def fmt(v):
    if abs(v) >= 1e6:
        return f"{v:,.0f}"
    if abs(v) >= 100:
        return f"{v:.1f}"
    return f"{v:.4g}"


def main_raw(a):
    groups, order = raw(a[0], *(a[1:4] + [None] * (3 - len(a[1:4]))))
    print("input\top\tmetric\tn\tmedian\tmin\tmax\tspread\tunit")
    for k in order:
        v = groups[k]
        lo, hi = min(v), max(v)
        print("\t".join([k[0], k[1], k[2], str(len(v)), fmt(statistics.median(v)), fmt(lo), fmt(hi), f"{hi / lo:.3f}" if lo > 0 else "-", k[3]]))


AB = re.compile(r"^AB\s+(\S+)\s+A\s+(\S+)\s+B\s+(\S+)\s+B/A median-of-paired\s+(\S+)\s+p25\s+(\S+)\s+p75\s+(\S+)\s+n\s+(\d+)")


def main_ab(files):
    for path in files:
        per = {}
        label = None
        with open(path) as f:
            for line in f:
                if line.startswith("# "):
                    label = line.strip()[2:].rsplit(" proc ", 1)[0]
                    continue
                m = AB.match(line)
                if m:
                    per.setdefault((label, m.group(1)), []).append(tuple(float(m.group(i)) for i in (2, 3, 4)))
        print(f"# {path}")
        print("label\toperation\tprocesses\tA median\tB median\tB/A median [min-max]")
        for (label, op), rows in per.items():
            a = [r[0] for r in rows]
            b = [r[1] for r in rows]
            q = [r[2] for r in rows]
            print(f"{label}\t{op}\t{len(rows)}\t{fmt(statistics.median(a))}\t{fmt(statistics.median(b))}\t{statistics.median(q):.3f} [{min(q):.3f}-{max(q):.3f}]")


def ab_by_label(path):
    per = {}
    label = None
    with open(path) as f:
        for line in f:
            if line.startswith("# "):
                label = line.strip()[2:].rsplit(" proc ", 1)[0]
                continue
            m = AB.match(line)
            if m:
                per.setdefault((label, m.group(1)), []).append(tuple(float(m.group(i)) for i in (2, 3, 4)))
    return per


def main_gate(path):
    per = ab_by_label(path)

    def med(label, op, k):
        rows = per.get((label, op))
        return statistics.median([r[k] for r in rows]) if rows else None

    def rng(label, op):
        rows = per.get((label, op))
        return (min(r[2] for r in rows), max(r[2] for r in rows), len(rows)) if rows else None

    print("text\tthreads\toperation\tnull B/A [min-max]\tfile/mem, A = mem\tfile/mem, A = file (1 / B/A)\tbalanced file/mem\tTrusted/Checked\tfacade/index\tmem q/s (all mem passes, median)\tfile q/s (median)")
    for text in ("dna_10m", "prot_1m"):
        for t in ("4t", "1t"):
            for op in ("count_present", "locate", "count_absent"):
                base = f"{text}.{t}."
                null = rng(base + "null_mem_mem", op)
                o1 = med(base + "mem_file", op, 2)
                o2 = med(base + "file_mem", op, 2)
                o2i = 1 / o2 if o2 else None
                bal = (o1 * o2i) ** 0.5 if o1 and o2i else None
                tr = med(base + "file_filetrusted", op, 2)
                ix = med(base + "idx_file", op, 2)
                memq = [r[0] for r in per.get((base + "mem_file", op), [])] + [r[0] for r in per.get((base + "null_mem_mem", op), [])] + [r[1] for r in per.get((base + "null_mem_mem", op), [])] + [r[1] for r in per.get((base + "file_mem", op), [])]
                fileq = [r[1] for r in per.get((base + "mem_file", op), [])] + [r[0] for r in per.get((base + "file_mem", op), [])]

                def f(x, nd=3):
                    return "-" if x is None else f"{x:.{nd}f}"
                nn = f"{null[2]}: " + f"{f(med(base + 'null_mem_mem', op, 2))} [{null[0]:.3f}-{null[1]:.3f}]" if null else "-"
                print("\t".join([text, t[0], op, nn, f(o1), f(o2i), f(bal), f(tr), f(ix), fmt(statistics.median(memq)) if memq else "-", fmt(statistics.median(fileq)) if fileq else "-"]))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    if sys.argv[1] == "--abapp":
        main_ab(sys.argv[2:])
    elif sys.argv[1] == "--gate":
        main_gate(sys.argv[2])
    else:
        main_raw(sys.argv[1:])
