#!/usr/bin/env python3
"""Work package 1b.9b: the table of both instruments over the rows files that run_1b9b_bench.sh wrote.
For every text (the files cross_<text><tag>.tsv and inproc_<text><tag>.tsv of the rows directory):
  cross-process: the throughput ratio file-backed (Checked) / in-memory as the ratio of the medians over the processes, as the median over the rounds of the per-round ratio (paired), and the 95 % bootstrap interval of the first
                 (the processes of both variants resampled); also the absolute medians (queries per second) of both;
  in process:    the median over the processes of the per-process median ratio (Checked / memory, and Trusted / memory) of the passes.
usage: summarize_1b9b.py <rows directory> [tag]"""
import os
import random
import statistics
import sys
from collections import defaultdict

TEXTS = ["dna", "sigma5", "protu", "protr", "text90", "bytes256"]
LABEL = {"dna": "DNA (4)", "sigma5": "5 symbols", "protu": "protein uniform (20)", "protr": "protein realistic (20)", "text90": "text (90)", "bytes256": "bytes (256)"}


def load_cross(path):
    rows = []
    for line in open(path, encoding="utf-8", errors="replace"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "row" and len(f) >= 7:
            rows.append(dict(round=f[1], variant=f[2], v=(float(f[3]), float(f[4]), float(f[5]))))
    return rows


def load_inproc(path):
    by = defaultdict(list)
    for line in open(path, encoding="utf-8", errors="replace"):
        f = line.split()
        if len(f) >= 10 and f[0] == "pair":
            by[f[2]].append((float(f[7]), float(f[8]), float(f[9])))
    return by


def boot(a, b, k, n=3000, seed=1):
    rng = random.Random(seed)
    xs = [x["v"][k] for x in a]
    ys = [y["v"][k] for y in b]
    out = []
    for _ in range(n):
        out.append(statistics.median(rng.choices(xs, k=len(xs))) / statistics.median(rng.choices(ys, k=len(ys))))
    out.sort()
    return out[int(0.025 * n)], out[int(0.975 * n) - 1]


def stratified(a, b, k):
    """The ratio of the medians inside each of the two modes of the throughput (the processes of both variants are split at the threshold that maximises the between-group variance of the log throughput, Otsu), averaged over the two
    modes: the processes of a variant fall into a slow and a fast group independently of the variant, so the ratio of the medians over all processes depends on how many of each fell into either variant; inside a mode it does not.
    Exploratory (it needs both modes in both variants); returns None when a mode is empty in a variant."""
    import math
    vals = sorted(math.log(x["v"][k]) for x in a + b)
    best, cut = -1.0, None
    for i in range(3, len(vals) - 3):
        lo, hi = vals[:i], vals[i:]
        w = len(lo) * len(hi) * (statistics.mean(lo) - statistics.mean(hi)) ** 2
        if w > best:
            best, cut = w, (vals[i - 1] + vals[i]) / 2
    ratios = []
    for low in (True, False):
        ga = [x["v"][k] for x in a if (math.log(x["v"][k]) < cut) == low]
        gb = [x["v"][k] for x in b if (math.log(x["v"][k]) < cut) == low]
        if len(ga) < 2 or len(gb) < 2:
            return None
        ratios.append(statistics.median(ga) / statistics.median(gb))
    return statistics.mean(ratios)


def cross(path, ref="memory", var="checked"):
    rows = load_cross(path)
    by = defaultdict(list)
    rounds = defaultdict(dict)
    for r in rows:
        by[r["variant"]].append(r)
        rounds[r["variant"]][r["round"]] = r
    a, b = by.get(var, []), by.get(ref, [])
    if not a or not b:
        return None
    res = []
    for k in range(3):
        ratio = statistics.median(x["v"][k] for x in a) / statistics.median(x["v"][k] for x in b)
        pairs = [rounds[var][rd]["v"][k] / rounds[ref][rd]["v"][k] for rd in rounds[var] if rd in rounds[ref]]
        lo, hi = boot(a, b, k)
        res.append((ratio, statistics.median(pairs), lo, hi, statistics.median(x["v"][k] for x in b), statistics.median(x["v"][k] for x in a), stratified(a, b, k)))
    return len(a), len(b), res


def pairs(path):
    """The rows of the PairBench / ThreadBench output: per variant the list of (abs count, abs locate, abs absent, ratio count, ratio locate, ratio absent) over the processes."""
    by = defaultdict(list)
    for line in open(path, encoding="utf-8", errors="replace"):
        f = line.split()
        if len(f) >= 10 and f[0] in ("pair", "mt"):
            off = 3 if f[0] == "pair" else 4
            ratio = f.index("ratio_to_first:") + 1
            by[f[off - 1]].append((float(f[off]), float(f[off + 1]), float(f[off + 2]), float(f[ratio]), float(f[ratio + 1]), float(f[ratio + 2])))
    return by


def short(v):
    return os.path.basename(os.path.dirname(v)) if v.startswith("alc=") else v


def families(d, tag):
    """The other instruments of 1b.9b: the in-memory tree against e4a0b2da (old_*), the ablation (abl_*), the four threads (threads_*), the null experiment (null_*), 100 Mbp (cross_dna100 / inproc_dna100)."""
    for kind, title, ref in (("old", "in-memory tree against e4a0b2da (in process, throughput relative to the first variant)", None), ("abl", "ablation at the engine level (in process, relative to the facade of the in-memory tree)", None),
                             ("threads", "four threads on one tree (in process, relative to the first variant)", None)):
        for t in TEXTS + ["dna100"]:
            p = os.path.join(d, "%s_%s%s.tsv" % (kind, t, tag))
            if not os.path.exists(p):
                continue
            by = pairs(p)
            if not by:
                continue
            print("\n%s, %s" % (title, LABEL.get(t, t)))
            for v, rows in by.items():
                cols = list(zip(*rows))
                print("  %-34s n=%d  abs k/s %s  ratio %s  [count range %.3f-%.3f]" % (short(v), len(rows), " ".join("%5d" % (statistics.median(c) / 1000) for c in cols[:3]), " ".join("%.3f" % statistics.median(c) for c in cols[3:]), min(cols[3]), max(cols[3])))
    for t in TEXTS + ["dna100"]:
        p = os.path.join(d, "null_%s%s.tsv" % (t, tag))
        if os.path.exists(p):
            r = cross(p, ref="memory", var="memory@head2")
            if r:
                n, m, res = r
                print("\nnull experiment (the same build twice, cross-process), %s: n=%d, ratios %s" % (LABEL.get(t, t), n, " | ".join("%.3f paired %.3f [%.3f-%.3f]" % x[:4] for x in res)))
    p = os.path.join(d, "cross_dna100%s.tsv" % tag)
    if os.path.exists(p):
        r = cross(p)
        if r:
            n, m, res = r
            print("\n100 Mbp DNA cross-process, n=%d: %s" % (n, " | ".join("%.3f paired %.3f [%.3f-%.3f] mem %.0fk" % (x[0], x[1], x[2], x[3], x[4] / 1000) for x in res)))
    p = os.path.join(d, "inproc_dna100%s.tsv" % tag)
    if os.path.exists(p):
        by = load_inproc(p)
        for v in ("checked", "trusted"):
            if v in by:
                cols = list(zip(*by[v]))
                print("100 Mbp DNA in process (%d processes), %s / memory: %s" % (len(by[v]), v, " ".join("%.3f [%.3f-%.3f]" % (statistics.median(c), min(c), max(c)) for c in cols)))


def main():
    d = sys.argv[1]
    tag = sys.argv[2] if len(sys.argv) > 2 else ""
    if "--families" in sys.argv:
        families(d, tag if tag != "--families" else "")
        return
    names = ["count (present long)", "locate (present short)", "absent"]
    print("cross-process: Checked / in-memory throughput ratio (>1: the file-backed tree is faster); n = processes per variant")
    print("%-24s %4s | %s" % ("text", "n", " | ".join("%-44s" % nm for nm in names)))
    for t in TEXTS:
        p = os.path.join(d, "cross_%s%s.tsv" % (t, tag))
        if not os.path.exists(p):
            continue
        r = cross(p)
        if not r:
            continue
        n, m, res = r
        cells = " | ".join("%.3f paired %.3f [%.3f-%.3f] mem %6.0fk strat %s" % (x[0], x[1], x[2], x[3], x[4] / 1000, "%.3f" % x[6] if x[6] else "n/a") for x in res)
        print("%-24s %4d | %s" % (LABEL[t], n, cells))
    print()
    print("in process: throughput relative to the in-memory tree (median over the processes of the median over the passes of the per-pass ratio)")
    print("%-24s %4s | %s" % ("text", "n", " | ".join("%-30s" % nm for nm in names)))
    for t in TEXTS:
        p = os.path.join(d, "inproc_%s%s.tsv" % (t, tag))
        if not os.path.exists(p):
            continue
        by = load_inproc(p)
        for v in ("checked", "trusted"):
            if v not in by:
                continue
            cols = list(zip(*by[v]))
            cells = " | ".join("%.3f [%.3f-%.3f]%s" % (statistics.median(c), min(c), max(c), " " * 8) for c in cols)
            print("%-24s %4d | %s   (%s)" % (LABEL[t], len(by[v]), cells, v))


main()
