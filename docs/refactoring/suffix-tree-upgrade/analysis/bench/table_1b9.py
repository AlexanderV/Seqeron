#!/usr/bin/env python3
"""Work package 1b.9: the tables of the ADR-004 subsection and of the Log, made from the raw rows (runs_esa/raw_1b9.tsv): median [min-max] over the fresh processes.
usage: table_1b9.py <raw file> build|query|nav|open|entry"""
import statistics
import sys


def load(path):
    d = {}
    with open(path) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            try:
                v = float(p[6])
            except ValueError:
                continue
            d.setdefault((p[2], p[4], p[5]), []).append(v)
    return d


def cell(d, inp, op, metric, scale=1.0, nd=3, comma=False):
    v = d.get((inp, op, metric))
    if not v:
        return "-"
    v = [x * scale for x in v]
    f = (lambda x: f"{x:,.0f}") if comma else (lambda x: f"{x:.{nd}f}")
    return f"{f(statistics.median(v))} [{f(min(v))}-{f(max(v))}]"


def main():
    d = load(sys.argv[1])
    what = sys.argv[2]
    ins = ["dna_10m", "dna_1m", "rep_10m", "prot_1m", "dna_20m", "dna_50m"]
    if what == "build":
        print("| input | op | total s | SA s | LCP s | child s | VmHWM B/char | ru_maxrss MB | Core B/char (arrays) | n |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for op in ("build", "lean_build"):
            for x in ins:
                if (x, op, "build_total_s") not in d:
                    continue
                print(f"| {x} | {op} | {cell(d, x, op, 'build_total_s')} | {cell(d, x, op, 'build_sa_s')} | {cell(d, x, op, 'build_lcp_s')} | {cell(d, x, op, 'build_child_s')} | "
                      f"{cell(d, x, op, 'peak_rss_bytes_per_char', nd=2)} | {cell(d, x, op.replace('build', 'build') + '_process' if False else ('build_process' if op == 'build' else 'lean_process'), 'ext_peak_rss', 1e-6, 1)} | "
                      f"{cell(d, x, op, 'index_bytes_core_per_char', nd=2)} | {len(d[(x, op, 'build_total_s')])} |")
    elif what == "query":
        print("| input | count present q/s | count absent q/s | locate (FindAll) patterns/s | locate lazy patterns/s | n |")
        print("|---|---|---|---|---|---|")
        for x in ins:
            op = "layout_Tagged"
            if (x, op, "count_present_qps") in d:
                print(f"| {x} | {cell(d, x, op, 'count_present_qps', comma=True)} | {cell(d, x, op, 'count_absent_qps', comma=True)} | {cell(d, x, op, 'locate_findall_pps', comma=True)} | "
                      f"{cell(d, x, op, 'locate_lazy_pps', comma=True)} | {len(d[(x, op, 'count_present_qps')])} |")
    elif what == "nav":
        print("| input | op | first call s | steady s | results | peak RSS MB | n |")
        print("|---|---|---|---|---|---|---|")
        for x, op in (("dna_10m", "mem"), ("dna_10m", "mum"), ("dna_10m", "anchors"), ("rep_10m", "mem"), ("dna_10m", "ukk_mem")):
            o = "nav_esa" if not op.startswith("ukk_") else "ukk_nav_esa"
            if (x, o, f"{op.replace('ukk_', '')}_first_call_s") not in d:
                continue
            m = op.replace("ukk_", "")
            proc = f"nav_{m}_process" if not op.startswith("ukk_") else f"nav_{m}_ukk_process"
            print(f"| {x} | {op} | {cell(d, x, o, m + '_first_call_s')} | {cell(d, x, o, m + '_steady_median_s')} | {int(d[(x, o, m + '_results')][0])} | {cell(d, x, proc, 'ext_peak_rss', 1e-6, 0)} | {len(d[(x, o, m + '_first_call_s')])} |")
        ukk = ("dna_10m", "ukk_nav_ukkonen", "mem_first_call_s")
        if ukk in d:
            print(f"| dna_10m | Ukkonen, same process | {cell(d, 'dna_10m', 'ukk_nav_ukkonen', 'mem_first_call_s')} | {cell(d, 'dna_10m', 'ukk_nav_ukkonen', 'mem_steady_median_s')} | - | {cell(d, 'dna_10m', 'nav_mem_ukk_process', 'ext_peak_rss', 1e-6, 0)} | {len(d[ukk])} |")
    elif what == "open":
        print("| input | cache | first Open ms | first query ms | warm Open us | Sections open ms | Structure open s | n |")
        print("|---|---|---|---|---|---|---|---|")
        for x in ("dna_10m", "prot_1m"):
            for tag in ("warm", "cold", "pwarm", "pcold"):
                op = f"{tag}_open"
                if (x, op, "first_open_s") in d:
                    print(f"| {x} | {'file cached' if tag.endswith('warm') else 'page cache dropped'}{', pattern file read first' if tag.startswith('p') else ''} | {cell(d, x, op, 'first_open_s', 1e3, 2)} | {cell(d, x, op, 'first_query_s', 1e3, 2)} | "
                          f"{cell(d, x, op, 'warm_open_median_s', 1e6, 1)} | {cell(d, x, op, 'open_sections_s', 1e3, 1)} | {cell(d, x, op, 'open_structure_s', 1, 2)} | {len(d[(x, op, 'first_open_s')])} |")
    elif what == "entry":
        print(f"entry dna_10m: {cell(d, 'dna_10m', 'entry', 'entry_build_s')} s, first count {cell(d, 'dna_10m', 'entry', 'first_count_s', 1e3, 2)} ms")


main()
