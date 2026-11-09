#!/usr/bin/env python3
"""Aggregates the raw rows of prototypes P1 and P2 (work package 1b.4a): median over the fresh processes, min/max in brackets.

  aggregate.py p1 p1_links_raw.tsv [--out p1_summary.tsv]      per input and (kernel, source): ns per link, steps, build costs, the decision rule of ADR-001 section 5
  aggregate.py p2 p2_mem_raw.tsv                               per input and variant: cold / warm time, counters, Ukkonen reference, MEM-set equality

Raw files hold `# process input=... run=... options: ...` header lines and `R<TAB>op<TAB>metric<TAB>value<TAB>unit<TAB>notes` rows.
"""
import re
import statistics
import sys
from collections import defaultdict


def read(path):
    """-> list of (header dict, [(op, metric, value-string, unit, notes)])"""
    blocks = []
    cur = None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith("# process"):
            m = dict(re.findall(r"(\w+)=(\S+)", line.split("options:")[0]))
            m["options"] = line.split("options:", 1)[1].strip() if "options:" in line else ""
            cur = (m, [])
            blocks.append(cur)
        elif line.startswith("R\t") and cur is not None:
            f = line.split("\t")
            if len(f) >= 5:
                cur[1].append((f[1], f[2], f[3], f[4], f[5] if len(f) > 5 else ""))
    return blocks


def num(s):
    try:
        return float(s)
    except ValueError:
        return None


def med(vals):
    return statistics.median(vals) if vals else float("nan")


def span(vals, fmt="{:.0f}"):
    if not vals:
        return "-"
    return f"{fmt.format(med(vals))} [{fmt.format(min(vals))}..{fmt.format(max(vals))}]"


def collect(blocks):
    """-> {input: {(op, metric): [values per process]}}, {input: n processes}"""
    data = defaultdict(lambda: defaultdict(list))
    runs = defaultdict(int)
    for hdr, rows in blocks:
        x = hdr["input"] + (f"  [{hdr['options']}]" if hdr.get("options") else "")
        if rows:
            runs[x] += 1
        for op, metric, value, unit, notes in rows:
            v = num(value)
            data[x][(op, metric)].append(v if v is not None else value)
    return data, runs


COMBOS = [("a_scan", "isa"), ("a_scan", "psi"), ("b_descent", "isa"), ("b_descent", "psi"),
          ("c_rmq", "isa"), ("c_rmq", "psi"), ("c2_psv", "isa"), ("c2_psv", "psi")]


def p1(path, out):
    blocks = read(path)
    data, runs = collect(blocks)
    def sort_key(x):
        base = ["dna_10m", "rep_10m", "prot_1m", "a^n", "fibonacci", "thue-morse"]
        head = x.split("  [")[0]
        name, _, size = head.partition(":")
        return (base.index(name) if name in base else len(base), int(size) if size.isdigit() else 0, x)
    order = sorted(data, key=sort_key)
    rows_out = []
    summary = {}
    for x in order:
        d = data[x]
        print(f"\n## {x}  ({runs[x]} fresh processes)")
        print("nav: " + "; ".join(f"{m} {span(d[('nav', m)], '{:.4g}')}" for m in ["build_isa_s", "build_psi_s", "build_rmq_s", "isa_bytes_per_char", "psi_bytes_per_char", "rmq_bytes_per_char"] if ('nav', m) in d))
        print("timer floor ns: " + span(d.get(("timer", "per_link_timer_floor_ns"), []), "{:.1f}") + "; internal nodes: " + span(d.get(("nodes", "internal_nodes"), []), "{:.0f}"))
        print("| kernel/source | median ns (net) | p99 ns (net) | max ns (raw, max over runs) | mean batch ns | steps median | steps p99 | steps max | budget exceeded | first link ns (cold) |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for k, s in COMBOS:
            op = f"sample_{k}_{s}"
            if (op, "median_ns") not in d and (op, "steps_median") not in d:
                continue
            g = lambda m: d.get((op, m), [])
            mx = max(g("max_ns")) if g("max_ns") else float("nan")
            smx = max(g("steps_max")) if g("steps_max") else float("nan")
            print(f"| {k} {s} | {span(g('median_net_ns'))} | {span(g('p99_net_ns'))} | {mx:.0f} | {span(g('mean_batch_ns'))} | {span(g('steps_median'), '{:.1f}')} | {span(g('steps_p99'), '{:.0f}')} | {smx:.0f} | {span(g('steps_budget_exceeded'), '{:.0f}')} of {span(g('steps_links'), '{:.0f}')} | {span(g('first_link_ns'), '{:.0f}')} |")
            summary[(x, k, s)] = dict(median=med(g("median_net_ns")), p99=med(g("p99_net_ns")), maxns=mx, steps_max=smx,
                                      exceeded=med(g("steps_budget_exceeded")), counted=med(g("steps_links")), steps_median=med(g("steps_median")))
            rows_out.append((x, k, s, med(g("median_net_ns")), med(g("p99_net_ns")), mx, med(g("mean_batch_ns")), med(g("steps_median")), med(g("steps_p99")), smx, med(g("steps_budget_exceeded")), med(g("steps_links")), med(g("first_link_ns"))))
        for key in sorted({o for (o, m) in d if o.startswith("worst_")} ):
            pass
        worst_ops = sorted({o for (o, m) in d if o.startswith("worst_") and m == "steps_max"})
        if worst_ops:
            print("\nworst-node sweeps (steps max over runs / median ns net, median of runs):")
            print("| list / kernel / source | nodes timed | steps max | steps median | median ns (net) | p99 ns (net) | max ns (raw) | budget exceeded |")
            print("|---|---|---|---|---|---|---|---|")
            for op in worst_ops:
                g = lambda m: d.get((op, m), [])
                print(f"| {op[6:]} | {span(g('timed_links'), '{:.0f}')} | {max(g('steps_max')):.0f} | {span(g('steps_median'), '{:.1f}')} | {span(g('median_net_ns'))} | {span(g('p99_net_ns'))} | {(max(g('max_ns')) if g('max_ns') else float('nan')):.0f} | {span(g('steps_budget_exceeded'), '{:.0f}')} |")
        reps = sorted({o for (o, m) in d if o.startswith("repeat_") and m == "median_net_ns"})
        for op in reps:
            first = op.replace("repeat_", "sample_")
            print(f"\ndrift check: {op[7:]} median ns (net) first {span(d.get((first, 'median_net_ns'), []))} -> repeated after all other kernels {span(d.get((op, 'median_net_ns'), []))}")
        sw = [(m, d[("sweeps", m)]) for m in ["max_scan_cost", "max_target_width", "max_string_depth", "max_node_depth"] if ("sweeps", m) in d]
        if sw:
            print("\nsweeps: " + "; ".join(f"{m} {span(v)}" for m, v in sw))
        mism = [(op, m, v) for (op, m), v in d.items() if m in ("cross_check_different", "mismatches") and any(t for t in v if t)]
        print("cross-check mismatches: " + (str(mism) if mism else "none"))
    print("\n## decision rule (ADR-001 section 5, pre-registered): ratio of the median ns/link (net of the timer floor) to kernel (c)")
    print("| input | source | (a)/(c) | (b)/(c) | (c2)/(c) | (c) median ns | (a) steps max | (b) steps max | (a) budget | (b) budget |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for x in order:
        for s in ["isa", "psi"]:
            c = summary.get((x, "c_rmq", s))
            if not c:
                continue
            a = summary.get((x, "a_scan", s)); b = summary.get((x, "b_descent", s)); c2 = summary.get((x, "c2_psv", s))
            r = lambda t: f"{t['median'] / c['median']:.2f}" if t else "-"
            bud = lambda t: f"{t['exceeded']:.0f} of {t['counted']:.0f}" if t else "-"
            print(f"| {x} | {s} | {r(a)} | {r(b)} | {r(c2)} | {c['median']:.0f} | {a['steps_max'] if a else float('nan'):.0f} | {b['steps_max'] if b else float('nan'):.0f} | {bud(a)} | {bud(b)} |")
    # ISA against Psi: ratio per fresh process (both were measured in the same process), median of the ratios
    print("\n## ISA against Psi: median over the fresh processes of (Psi median ns / ISA median ns), per process; [min..max]")
    print("| input | (a) scan | (b) descent | (c) rmq | (c2) psv |")
    print("|---|---|---|---|---|")
    for x in order:
        d = data[x]
        cells = []
        for k in ["a_scan", "b_descent", "c_rmq", "c2_psv"]:
            ia = d.get((f"sample_{k}_isa", "median_net_ns"), [])
            ps = d.get((f"sample_{k}_psi", "median_net_ns"), [])
            ratios = [p / i for i, p in zip(ia, ps) if i and i > 0]
            cells.append(span(ratios, "{:.2f}") if ratios else "-")
        print(f"| {x} | " + " | ".join(cells) + " |")

    # boundedness: largest step counts of the sample and of the worst-node sweeps, and the exact worst scan cost over all nodes, by text size
    print("\n## boundedness by text size (steps; max over the fresh processes; the sweeps include the worst nodes by scan cost, target width, string depth, node depth)")
    print("| input | max scan cost over ALL nodes | max node depth | (a) steps max sample+worst | (b) steps max sample+worst | (c) steps max sample+worst | (c2) steps max sample+worst |")
    print("|---|---|---|---|---|---|---|")
    for x in order:
        d = data[x]
        def worst_steps(k, src="psi"):
            vals = []
            for (op, m), v in d.items():
                if m == "steps_max" and (op == f"sample_{k}_{src}" or (op.startswith("worst_") and op.endswith(f"_{k}_{src}"))):
                    vals += [t for t in v if isinstance(t, float)]
            return max(vals) if vals else float("nan")
        sc = d.get(("sweeps", "max_scan_cost"), [float("nan")])
        nd = d.get(("sweeps", "max_node_depth"), [float("nan")])
        print(f"| {x} | {max(sc):.0f} | {max(nd):.0f} | {worst_steps('a_scan'):.0f} | {worst_steps('b_descent'):.0f} | {worst_steps('c_rmq'):.0f} | {worst_steps('c2_psv'):.0f} |")

    if out:
        with open(out, "w", encoding="utf-8") as f:
            f.write("# P1 summary: medians over fresh processes; columns: input kernel source median_net_ns p99_net_ns max_ns_raw mean_batch_ns steps_median steps_p99 steps_max budget_exceeded counted first_link_ns\n")
            for r in rows_out:
                f.write("\t".join(str(v) if not isinstance(v, float) else f"{v:.6g}" for v in r) + "\n")


def short(x):
    """Summary key: the input name without the size of a synthetic text, except when sizes differ (scaling rows)."""
    return x


def p2(path):
    blocks = read(path)
    groups = defaultdict(list)
    for hdr, rows in blocks:
        groups[(hdr["input"], hdr["options"])].append(rows)
    for (x, opts), runs in groups.items():
        print(f"\n## {x}  options: {opts or '(defaults)'}  ({len(runs)} fresh processes)")
        def col(op, metric):
            out = []
            for rows in runs:
                for o, m, v, u, n in rows:
                    if o == op and m == metric:
                        t = num(v)
                        if t is not None:
                            out.append(t)
            return out
        for op, metric in [("mem", "build_esa_s"), ("mem", "build_nav_s"), ("mem_esa", "first_call_s"), ("mem_esa", "warm_median_s"), ("mem_esa", "steady_median_s"), ("mem_esa", "matches"), ("mem_esa", "links_followed"),
                           ("mem_esa", "extend_steps"), ("mem_esa", "child_lookups"), ("mem_esa", "rescan_lookups"), ("mem_esa", "link_kernel_share"), ("mem_esa", "link_kernel_ns_per_link"),
                           ("mem_ukkonen", "build_s"), ("mem_ukkonen", "first_call_s"), ("mem_ukkonen", "warm_median_s"), ("mem_ukkonen", "steady_median_s"), ("mem_ukkonen", "matches"), ("mem_check", "mem_set_equal"),
                           ("mem_esa", "peak_rss_hwm"), ("mem", "peak_rss_hwm_final")]:
            v = col(op, metric)
            if v:
                print(f"  {op}.{metric}: {span(v, '{:.4g}')}")


def p2_pairs(path):
    """Paired comparison of variants that differ only in one option, pairing the fresh processes of the same repetition (run number) and input."""
    blocks = read(path)
    by = defaultdict(dict)       # (input, run) -> {options: {metric: value}}
    sums = defaultdict(set)
    for hdr, rows in blocks:
        d = {}
        for op, metric, value, unit, notes in rows:
            t = num(value)
            if t is not None:
                d[(op, metric)] = t
        by[(hdr["input"], hdr["run"])][hdr["options"]] = d
        if ("mem_esa", "checksum") in d:
            sums[hdr["input"]].add(d[("mem_esa", "checksum")])
    print("\n## checksums of the MEM sets over all variants and runs (one value per input = identical sets)")
    for x, v in sums.items():
        print(f"  {x}: {sorted(v)}")
    def key(o, **kv):
        for k, v in kv.items():
            if f"--{k} {v}" not in o:
                return False
        return True
    pairs = [
        ("kernel rmq: ISA / Psi (warm)", dict(kernel="rmq", loci="2"), "source", "isa", "psi"),
        ("kernel psv: ISA / Psi (warm)", dict(kernel="psv", loci="2"), "source", "isa", "psi"),
        ("psv psi: loci 1 / loci 2 (warm)", dict(kernel="psv", source="psi"), "loci", "1", "2"),
        ("psi, loci 2: psv / rmq (warm)", dict(source="psi", loci="2"), "kernel", "psv", "rmq"),
        ("psi, loci 2: scan / rmq (warm)", dict(source="psi", loci="2"), "kernel", "scan", "rmq"),
        ("psi, loci 2: descent / rmq (warm)", dict(source="psi", loci="2"), "kernel", "descent", "rmq"),
    ]
    print("\n## paired ratios of the warm median (steady median of the last five calls), per repetition; median [min..max]")
    for title, fixed, opt, a, b in pairs:
        for inp in ("dna_10m", "rep_10m"):
            ratios, rc = [], []
            for (x, run), variants in by.items():
                if x != inp:
                    continue
                va = [d for o, d in variants.items() if key(o, **fixed) and f"--{opt} {a}" in o]
                vb = [d for o, d in variants.items() if key(o, **fixed) and f"--{opt} {b}" in o]
                if va and vb and ("mem_esa", "steady_median_s") in va[0] and ("mem_esa", "steady_median_s") in vb[0]:
                    ratios.append(va[0][("mem_esa", "steady_median_s")] / vb[0][("mem_esa", "steady_median_s")])
                    rc.append(va[0][("mem_esa", "first_call_s")] / vb[0][("mem_esa", "first_call_s")])
            if ratios:
                print(f"  {inp}: {title}: warm {span(ratios, '{:.2f}')}; first call {span(rc, '{:.2f}')}  (n = {len(ratios)})")


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("p1", "p2"):
        sys.exit(__doc__)
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
    if sys.argv[1] == "p1":
        p1(sys.argv[2], out)
    else:
        p2(sys.argv[2])
        p2_pairs(sys.argv[2])
