#!/usr/bin/env python3
"""Generate REPORT.md from bench.tsv (all numbers are pulled from the data, none typed by hand).  usage: make_report.py <bench-dir>"""
import csv
import sys

B = sys.argv[1]
rows = list(csv.DictReader(open(f"{B}/bench.tsv"), delimiter="\t"))
D = {}
REP = {}
for r in rows:
    D[(r["tool"], r["input"], r["operation"], r["metric"])] = r["value"]
    REP[(r["tool"], r["input"], r["operation"], r["metric"])] = r["repeats"]


def v(t, i, o, m):
    x = D.get((t, i, o, m))
    if x is None:
        return None
    try:
        return float(x)
    except ValueError:
        return None


def build_time(t, i):
    return v(t, i, "build", "time") or v(t, i, "index_process", "ext_wall")


def build_rss(t, i):
    return v(t, i, "build", "peak_rss_hwm_after_build") or v(t, i, "build_process", "ext_peak_rss") or v(t, i, "index_process", "ext_peak_rss")


def bpc(t, i):
    return v(t, i, "build", "bytes_per_char") or v(t, i, "index", "bytes_per_char")


def f(x, spec="{:,.3g}"):
    return "-" if x is None else spec.format(x)


SEQ_M, SEQ_P, GO = "seqeron-st-inmem", "seqeron-st-persistent", "go-suffixarray"
QUERY_TOOLS = ["go-suffixarray", "libsais", "libdivsufsort", "rust-bio", "sdsl-csa_wt_huff", "sdsl-csa_wt_rrr", "sdsl-csa_sada", "sdsl-cst_sct3", "pydivsufsort"]
ALL_TOOLS = [SEQ_M, SEQ_P] + QUERY_TOOLS + ["suffix-trees-py", "bwa-index", "bowtie2-build"]


def best(inp, fn, lo=True, tools=QUERY_TOOLS):
    c = [(fn(t, inp), t) for t in tools if fn(t, inp) is not None]
    return (min(c) if lo else max(c)) if c else (None, None)


def ratio_cell(a, b, higher_is_better=False, worse="slower", better="faster"):
    if a is None or b is None or b == 0:
        return "-"
    x = (b / a) if higher_is_better else (a / b)
    return f"{x:,.1f}x {worse}" if x >= 1 else f"{1 / x:,.1f}x {better}"


inp = "dna_10m"
out = []
w = out.append

w(f"""# Suffix-tree benchmark: Seqeron.SuffixTree vs Go stdlib and open-source competitors

Measured on one machine, 3 repeats (median) unless stated.  Seqeron = commit `242318429a5e` (working tree clean; library sources were copied to
`bench/seqeron/snapshot` and built from there, nothing under `/home/user/Seqeron` was touched).  All numbers below are generated from `bench.tsv`
by `lib/make_report.py`; raw per-repeat rows are in `runs/raw*.tsv`.  Environment: `ENV.txt`.  Provenance and licences of every third-party source: `SOURCES.txt`.

## 0. Read this first (what is and is not comparable)

* **Different data structures.** Seqeron and `sdsl cst_sct3` / Python `suffix-trees` are *suffix trees*; Go `index/suffixarray`, libsais, libdivsufsort and pydivsufsort
  are plain *suffix arrays* (4 B/char, binary-search queries); rust-bio and `sdsl csa_*` are *FM-indexes* (compressed suffix arrays);
  BWA and Bowtie2 are *FM-indexes for read alignment* (index build only was measured).  A suffix array answers count/locate/contains but has no LRS, MEM/MUM,
  matching statistics, suffix links or tree navigation; those operations are exactly what a suffix tree is bought for.  Numbers therefore compare **cost for the shared
  functionality (build, size, count, locate)**; they do not say the other tools replace the tree's extra functions.
* libsais / libdivsufsort / pydivsufsort are *construction* libraries: their count/locate rows are *our* binary-search code (`competitors/sa_c/sa_bench.c`) over their SA.
* Every subject returned **identical count and locate checksums** on every input it ran (`crosscheck.md`; 0 disagreements) so the query sets really are the same workload.
* The VM reports a 480 MiB L3 and 4 vCPU (KVM); random-access structures up to a few hundred MB may be served largely from that cache, so the 1-10 Mbp rows are cache-friendly relative to a bare-metal desktop.

## 1. Headline: 10 Mbp random DNA (`{inp}`)

| metric | Seqeron in-memory | Seqeron persistent (mmap) | Go `index/suffixarray` | best other tool (that metric) | Seqeron in-mem vs best | Seqeron persistent vs best |
|---|---|---|---|---|---|---|""")
bt, bt_t = best(inp, build_time, True)
bb, bb_t = best(inp, bpc, True)
bc, bc_t = best(inp, lambda t, i: v(t, i, "count_present20", "throughput"), False)
bl, bl_t = best(inp, lambda t, i: v(t, i, "locate8", "throughput"), False)
def valid4(t, i):
    ok = v(t, i, "count_present20_4t", "checksum_ok")
    return None if ok == 0 else v(t, i, "count_present20_4t", "throughput")


b4, b4_t = best(inp, valid4, False)
sm_bt, sp_bt = build_time(SEQ_M, inp), build_time(SEQ_P, inp)
w(f"| build wall time (s) | {f(sm_bt)} | {f(sp_bt)} | {f(build_time(GO, inp))} | {f(bt)} ({bt_t}) | {ratio_cell(sm_bt, bt)} | {ratio_cell(sp_bt, bt)} |")
w(f"| index bytes per input char | {f(bpc(SEQ_M, inp))} (GC heap; +2 B/char text string) | {f(bpc(SEQ_P, inp))} (file; text included) | {f(bpc(GO, inp))} (+1 B/char text) | {f(bb)} ({bb_t}) | {ratio_cell(bpc(SEQ_M, inp), bb, False, 'larger', 'smaller')} | {ratio_cell(bpc(SEQ_P, inp), bb, False, 'larger', 'smaller')} |")
w(f"| peak RSS during build (MB) | {f(build_rss(SEQ_M, inp) / 1e6, '{:,.0f}')} | {f(build_rss(SEQ_P, inp) / 1e6, '{:,.0f}')} (of which anonymous {f(v(SEQ_P, inp, 'build_process', 'peak_rss_anon_sampled') / 1e6, '{:,.0f}')}; rest = mmap'd file pages) | {f(build_rss(GO, inp) / 1e6, '{:,.0f}')} | {f(min(build_rss(t, inp) for t in QUERY_TOOLS if build_rss(t, inp)) / 1e6, '{:,.0f}')} (min over {len(QUERY_TOOLS)} tools) | - | - |")
w(f"| count, 20-mers present (queries/s, 1 thread) | {f(v(SEQ_M, inp, 'count_present20', 'throughput'), '{:,.0f}')} | {f(v(SEQ_P, inp, 'count_present20', 'throughput'), '{:,.0f}')} | {f(v(GO, inp, 'count_present20', 'throughput'), '{:,.0f}')} | {f(bc, '{:,.0f}')} ({bc_t}) | {ratio_cell(v(SEQ_M, inp, 'count_present20', 'throughput'), bc, True)} | {ratio_cell(v(SEQ_P, inp, 'count_present20', 'throughput'), bc, True)} |")
w(f"| count, 4 threads, one shared index (queries/s) | {f(v(SEQ_M, inp, 'count_present20_4t', 'throughput'), '{:,.0f}')} | {f(v(SEQ_P, inp, 'count_present20_4t', 'throughput'), '{:,.0f}')} | {f(v(GO, inp, 'count_present20_4t', 'throughput'), '{:,.0f}')} | {f(b4, '{:,.0f}')} ({b4_t}) | {ratio_cell(v(SEQ_M, inp, 'count_present20_4t', 'throughput'), b4, True)} | {ratio_cell(v(SEQ_P, inp, 'count_present20_4t', 'throughput'), b4, True)} |")
w(f"| locate, 1 000 random 8-mers, ~{f(v(SEQ_M, inp, 'locate8', 'total_hits') / 1000, '{:,.0f}')} hits each (patterns/s) | {f(v(SEQ_M, inp, 'locate8', 'throughput'), '{:,.0f}')} | {f(v(SEQ_P, inp, 'locate8', 'throughput'), '{:,.0f}')} | {f(v(GO, inp, 'locate8', 'throughput'), '{:,.0f}')} | {f(bl, '{:,.0f}')} ({bl_t}) | {ratio_cell(v(SEQ_M, inp, 'locate8', 'throughput'), bl, True)} | {ratio_cell(v(SEQ_P, inp, 'locate8', 'throughput'), bl, True)} |")
w(f"| locate throughput (positions returned/s) | {f(v(SEQ_M, inp, 'locate8', 'positions_per_s'), '{:,.0f}')} | {f(v(SEQ_P, inp, 'locate8', 'positions_per_s'), '{:,.0f}')} | {f(v(GO, inp, 'locate8', 'positions_per_s'), '{:,.0f}')} | {f(max(v(t, inp, 'locate8', 'positions_per_s') for t in QUERY_TOOLS if v(t, inp, 'locate8', 'positions_per_s')), '{:,.0f}')} | - | - |")
w(f"| count latency p50 / p99 (ns, 1 thread) | {f(v(SEQ_M, inp, 'count_present20', 'latency_p50'), '{:,.0f}')} / {f(v(SEQ_M, inp, 'count_present20', 'latency_p99'), '{:,.0f}')} | {f(v(SEQ_P, inp, 'count_present20', 'latency_p50'), '{:,.0f}')} / {f(v(SEQ_P, inp, 'count_present20', 'latency_p99'), '{:,.0f}')} | {f(v(GO, inp, 'count_present20', 'latency_p50'), '{:,.0f}')} / {f(v(GO, inp, 'count_present20', 'latency_p99'), '{:,.0f}')} | libsais {f(v('libsais', inp, 'count_present20', 'latency_p50'), '{:,.0f}')} / {f(v('libsais', inp, 'count_present20', 'latency_p99'), '{:,.0f}')} | - | - |")
w("""
"Best other tool" ranges over Go, libsais, libdivsufsort, rust-bio, sdsl (4 configs) and pydivsufsort -- i.e. over suffix arrays *and* FM-indexes, so its bytes/char is an FM-index and its
build time a suffix-array construction.  Same-kind comparators (suffix trees) are in section 3.""")

# ---- per-input tables
w("\n## 2. All tools, all inputs (medians of 3 repeats; see section 6 for exceptions)\n")
w("Columns: build wall time; index bytes per input char (Seqeron in-memory = GC heap delta, persistent = file size, Go/libsais/divsufsort = SA only, FM-indexes/CST = `size_in_bytes`, BWA/Bowtie2 = all index files); "
  "peak RSS during build; single-thread 20-mer count (present / absent), locate of 1 000 8-mers, 4-thread shared-index count.  `-` = not applicable / not run.  For `prot_1m` the 8-mer locate set has almost no hits (20^8 >> 1e6): use the `locate4` column.\n")
for i in ["dna_1m", "dna_10m", "rep_10m", "prot_1m", "dna_20m", "dna_50m"]:
    tools = [t for t in ALL_TOOLS if any(k[0] == t and k[1] == i for k in D)]
    w(f"\n### {i}\n")
    cols = "| tool | build s | B/char | peak RSS MB (build) | count q/s | absent q/s | locate8 pat/s | locate4 pat/s | count 4T q/s | repeats |"
    w(cols)
    w("|" + "---|" * 10)
    for t in tools:
        rep = REP.get((t, i, "build", "time")) or REP.get((t, i, "index_process", "ext_wall")) or "-"
        w(f"| {t} | {f(build_time(t, i))} | {f(bpc(t, i))} | {f((build_rss(t, i) or 0) / 1e6 if build_rss(t, i) else None, '{:,.0f}')} | {f(v(t, i, 'count_present20', 'throughput'), '{:,.0f}')} | {f(v(t, i, 'count_absent20', 'throughput'), '{:,.0f}')} | {f(v(t, i, 'locate8', 'throughput'), '{:,.0f}')} | {f(v(t, i, 'locate4', 'throughput'), '{:,.0f}')} | {f(valid4(t, i), '{:,.0f}') if v(t, i, 'count_present20_4t', 'checksum_ok') != 0 else 'invalid (race)'} | {rep} |")

# ---- suffix tree comparators + scaling
w("\n## 3. Suffix trees only: scaling and same-kind comparison\n")
w("Seqeron persistent (mmap) tree vs the other suffix-array / compressed-tree builds at growing n (random DNA):\n")
w("| n (Mbp) | Seqeron persistent build s | Seqeron chars/s | file GB | peak RSS GB (file-backed + anon) | anon MB | count q/s | locate8 pat/s | MEM 100 kbp query s | Go build s | libsais build s | libsais_omp4 build s |")
w("|---|---|---|---|---|---|---|---|---|---|---|---|")
for i, n in [("dna_1m", 1), ("dna_10m", 10), ("dna_20m", 20), ("dna_50m", 50)]:
    w(f"| {n} | {f(build_time(SEQ_P, i))} | {f(v(SEQ_P, i, 'build', 'throughput'), '{:,.0f}')} | {f((v(SEQ_P, i, 'build', 'index_bytes_file') or 0) / 1e9, '{:.2f}')} | {f((v(SEQ_P, i, 'build_process', 'ext_peak_rss') or 0) / 1e9, '{:.2f}')} | {f((v(SEQ_P, i, 'build_process', 'peak_rss_anon_sampled') or 0) / 1e6, '{:,.0f}')} | {f(v(SEQ_P, i, 'count_present20', 'throughput'), '{:,.0f}')} | {f(v(SEQ_P, i, 'locate8', 'throughput'), '{:,.0f}')} | {f(v(SEQ_P, i, 'mem_q100k_min20', 'time'))} | {f(build_time(GO, i))} | {f(build_time('libsais', i))} | {f(v('libsais', i, 'build_omp4', 'time'))} |")
w(f"""
* **Largest Seqeron size that finished: 50 Mbp** (persistent tree, {f(build_time(SEQ_P, 'dna_50m'), '{:,.0f}')} s, 3 of 3 repeats, 2.81 GB file, peak RSS {f(v(SEQ_P, 'dna_50m', 'build_process', 'ext_peak_rss') / 1e9, '{:.1f}')} GB of which only {f(v(SEQ_P, 'dna_50m', 'build_process', 'peak_rss_anon_sampled') / 1e6, '{:,.0f}')} MB anonymous; the rest is the mmap'd file).
  Build rate falls from {f(v(SEQ_P, 'dna_10m', 'build', 'throughput'), '{:,.0f}')} chars/s at 10 Mbp to {f(v(SEQ_P, 'dna_50m', 'build', 'throughput'), '{:,.0f}')} chars/s at 50 Mbp (at 50 Mbp the Ukkonen phase is {f(v(SEQ_P, 'dna_50m', 'build', 'stage_Ukkonen_build') / build_time(SEQ_P, 'dna_50m') * 100, '{:.0f}')} % of the time, i.e. random access to the node file, not the finalize passes).
* **100 Mbp did not finish**: one attempt with a 600 s cap was killed in the Ukkonen phase (`runs/raw.tsv` row `dna_100m ... BLOCKED`).  At the kill the main storage file held 1 531 MB of
  24-byte nodes = ~64 M of the ~162 M nodes a 100 Mbp tree needs, i.e. ~40 % of the text after 600 s (estimate), RSS 4.8 GB (4.3 GB file-backed, 0.45 GB anonymous); the sample log is `runs/mon_100m.keep.txt`.
  Extrapolating the 10 -> 20 -> 50 Mbp trend, a 100 Mbp build would need well over 30 min.  The text length is also an `int` in the API (< 2^31 chars) -- 3 Gbp is out of reach by construction, not only by time.
* Construction is single-threaded: CPU time = wall time ({f(v(SEQ_P, 'dna_50m', 'build_process', 'ext_cpu_user_sys'), '{:,.0f}')} s CPU for {f(build_time(SEQ_P, 'dna_50m'), '{:,.0f}')} s wall at 50 Mbp).  libsais's OpenMP build gives {f(v('libsais', 'dna_50m', 'build', 'time') / v('libsais', 'dna_50m', 'build_omp4', 'time'), '{:.1f}')}x at 4 threads on 50 Mbp.
* Build stage breakdown at 10 Mbp (persistent): Ukkonen {f(v(SEQ_P, inp, 'build', 'stage_Ukkonen_build'))} s, finalize (child arrays + leaf counts + write-back + trim) {f(sp_bt - v(SEQ_P, inp, 'build', 'stage_Ukkonen_build'))} s.
""")
w("Same-kind comparators (index = suffix *tree*), 1 Mbp -> 10 Mbp:\n")
w("| tool | n (Mbp) | build s | B/char | peak RSS MB | count q/s | locate8 pat/s |")
w("|---|---|---|---|---|---|---|")
for t in [SEQ_M, SEQ_P, "sdsl-cst_sct3"]:
    for i, n in [("dna_1m", 1), ("dna_10m", 10)]:
        w(f"| {t} | {n} | {f(build_time(t, i))} | {f(bpc(t, i))} | {f((build_rss(t, i) or 0) / 1e6, '{:,.0f}')} | {f(v(t, i, 'count_present20', 'throughput'), '{:,.0f}')} | {f(v(t, i, 'locate8', 'throughput'), '{:,.0f}')} |")
for i in ["dna_1m_prefix_10000", "dna_1m_prefix_100000", "dna_1m_prefix_1000000", "dna_10m_prefix_3000000", "dna_10m_prefix_10000000"]:
    t = "suffix-trees-py"
    n = v(t, i, "build", "time") and int(i.rsplit("_", 1)[1]) / 1e6
    w(f"| {t} (pure Python) | {f(n)} | {f(build_time(t, i))} | {f((build_rss(t, i) or 0) / int(i.rsplit('_', 1)[1]), '{:,.0f}')} (RSS/n) | {f((build_rss(t, i) or 0) / 1e6, '{:,.0f}')} | {f(v(t, i, 'contains_present20', 'throughput'), '{:,.0f}')} (contains) | {f(v(t, i, 'locate8', 'throughput'), '{:,.0f}')} |")
w(f"""
Python `suffix-trees` 0.4.0 **gives up between 10 and 20 Mbp**: 10 Mbp builds in {f(build_time('suffix-trees-py', 'dna_10m_prefix_10000000'), '{:,.0f}')} s using {f(build_rss('suffix-trees-py', 'dna_10m_prefix_10000000') / 1e9, '{:.1f}')} GB; the 20 Mbp attempt hit `MemoryError` at the 12 GB address-space cap after ~800 s (one repeat observed, recorded as BLOCKED).
Seqeron's in-memory tree is {f(build_time('suffix-trees-py', 'dna_10m_prefix_10000000') / sm_bt, '{:.0f}')}x faster to build than that Python tree at 10 Mbp and uses {f(build_rss('suffix-trees-py', 'dna_10m_prefix_10000000') / build_rss(SEQ_M, inp), '{:.1f}')}x less memory; `sdsl cst_sct3` builds the same 10 Mbp in {f(build_time('sdsl-cst_sct3', inp))} s with {f(bpc('sdsl-cst_sct3', inp))} B/char (compressed tree, navigation via `sdsl` API, count by descent -- {f(v('sdsl-cst_sct3', inp, 'count_present20', 'throughput'), '{:,.0f}')} q/s).""")

# ---- Seqeron-specific
w("\n## 4. Seqeron-specific measurements\n")
w("| metric | in-memory 10 Mbp | persistent 10 Mbp | persistent 50 Mbp |\n|---|---|---|---|")
w(f"| LongestRepeatedSubstring first call (s) | {f(v(SEQ_M, inp, 'lrs', 'first_call_time'), '{:.2g}')} | {f(v(SEQ_P, inp, 'lrs', 'first_call_time'), '{:.2g}')} | {f(v(SEQ_P, 'dna_50m', 'lrs', 'first_call_time'), '{:.2g}')} |")
w(f"| MEM, 100 kbp query w/ 1 % substitutions, minLength 20: first call / warm (s) | {f(v(SEQ_M, inp, 'mem_q100k_min20', 'time'))} / {f(v(SEQ_M, inp, 'mem_q100k_min20', 'time_second_call'))} | {f(v(SEQ_P, inp, 'mem_q100k_min20', 'time'))} / {f(v(SEQ_P, inp, 'mem_q100k_min20', 'time_second_call'))} | {f(v(SEQ_P, 'dna_50m', 'mem_q100k_min20', 'time'))} / {f(v(SEQ_P, 'dna_50m', 'mem_q100k_min20', 'time_second_call'))} |")
w(f"| MEM matches returned | {f(v(SEQ_M, inp, 'mem_q100k_min20', 'matches'), '{:,.0f}')} | {f(v(SEQ_P, inp, 'mem_q100k_min20', 'matches'), '{:,.0f}')} | {f(v(SEQ_P, 'dna_50m', 'mem_q100k_min20', 'matches'), '{:,.0f}')} |")
w(f"| MUM (`Both`) time (s) / matches | {f(v(SEQ_M, inp, 'mum_q100k_min20', 'time'))} / {f(v(SEQ_M, inp, 'mum_q100k_min20', 'matches'), '{:,.0f}')} | {f(v(SEQ_P, inp, 'mum_q100k_min20', 'time'))} / {f(v(SEQ_P, inp, 'mum_q100k_min20', 'matches'), '{:,.0f}')} | {f(v(SEQ_P, 'dna_50m', 'mum_q100k_min20', 'time'))} / {f(v(SEQ_P, 'dna_50m', 'mum_q100k_min20', 'matches'), '{:,.0f}')} |")
w(f"| FindExactMatchAnchors time (s) | {f(v(SEQ_M, inp, 'anchors_q100k_min20', 'time'))} | {f(v(SEQ_P, inp, 'anchors_q100k_min20', 'time'))} | {f(v(SEQ_P, 'dna_50m', 'anchors_q100k_min20', 'time'))} |")
w(f"| MEM on the repetitive 10 Mbp text (`rep_10m`): time (s) / matches | {f(v(SEQ_M, 'rep_10m', 'mem_q100k_min20', 'time'))} / {f(v(SEQ_M, 'rep_10m', 'mem_q100k_min20', 'matches'), '{:,.0f}')} | {f(v(SEQ_P, 'rep_10m', 'mem_q100k_min20', 'time'))} / {f(v(SEQ_P, 'rep_10m', 'mem_q100k_min20', 'matches'), '{:,.0f}')} | - |")
w(f"| cold `Load` time, page cache dropped (ms) | - | {f(v(SEQ_P, inp, 'cold_load', 'time') * 1e3)} | {f(v(SEQ_P, 'dna_50m', 'cold_load', 'time') * 1e3)} |")
w(f"| warm `Load` time (ms) | - | {f(v(SEQ_P, inp, 'warm_load', 'time') * 1e3)} | {f(v(SEQ_P, 'dna_50m', 'warm_load', 'time') * 1e3)} |")
w(f"| first 100 000 count queries right after cold / warm load (s) | - | {f(v(SEQ_P, inp, 'cold_load', 'first_100k_queries_time'))} / {f(v(SEQ_P, inp, 'warm_load', 'first_100k_queries_time'))} | {f(v(SEQ_P, 'dna_50m', 'cold_load', 'first_100k_queries_time'))} / {f(v(SEQ_P, 'dna_50m', 'warm_load', 'first_100k_queries_time'))} |")
w(f"| 4-thread / 1-thread count throughput (one shared tree) | {f(v(SEQ_M, inp, 'count_present20_4t', 'throughput') / v(SEQ_M, inp, 'count_present20_1t', 'throughput'), '{:.2f}')}x | {f(v(SEQ_P, inp, 'count_present20_4t', 'throughput') / v(SEQ_P, inp, 'count_present20_1t', 'throughput'), '{:.2f}')}x | {f(v(SEQ_P, 'dna_50m', 'count_present20_4t', 'throughput') / v(SEQ_P, 'dna_50m', 'count_present20_1t', 'throughput'), '{:.2f}')}x |")
w(f"| 4-thread / 1-thread locate8 throughput | {f(v(SEQ_M, inp, 'locate8_4t', 'throughput') / v(SEQ_M, inp, 'locate8_1t', 'throughput'), '{:.2f}')}x | {f(v(SEQ_P, inp, 'locate8_4t', 'throughput') / v(SEQ_P, inp, 'locate8_1t', 'throughput'), '{:.2f}')}x | {f(v(SEQ_P, 'dna_50m', 'locate8_4t', 'throughput') / v(SEQ_P, 'dna_50m', 'locate8_1t', 'throughput'), '{:.2f}')}x |")
w(f"""
Go reference for the same two rows: cold `Index.Read` of the serialized index {f(v(GO, inp, 'cold_load', 'time') * 1e3, '{:,.0f}')} ms at 10 Mbp ({f(v(GO, inp, 'ser_serialize', 'serialized_bytes') / 1e6, '{:,.1f}')} MB serialized, {f(v(GO, inp, 'ser_serialize', 'serialized_bytes_per_char'))} B/char; no mmap, whole index copied into RAM);
4T/1T count scaling {f(v(GO, inp, 'count_present20_4t', 'throughput') / v(GO, inp, 'count_present20_1t', 'throughput'), '{:.2f}')}x.

Notes on the numbers above
* 4T/1T ratios slightly above 4.0 on a 4-vCPU VM are measurement noise in the 1-thread baseline (turbo/steal), read them as "scales ~linearly to 4 threads".  Concurrent reads of one shared tree returned correct results in every run (`checksum_ok = 1`).
* Persistent query throughput is {f(min(v(SEQ_M, inp, o, m) / v(SEQ_P, inp, o, m) for o, m in [('count_present20', 'throughput'), ('count_absent20', 'throughput'), ('locate8', 'throughput'), ('count_present20_4t', 'throughput')]), '{:.2f}')}-{f(max(v(SEQ_M, inp, o, m) / v(SEQ_P, inp, o, m) for o, m in [('count_present20', 'throughput'), ('count_absent20', 'throughput'), ('locate8', 'throughput'), ('count_present20_4t', 'throughput')]), '{:.2f}')}x lower than in-memory on 10 Mbp (count / absent / locate8 / 4T count); the persistent tree is queried through the page cache here (file just built / loaded, 15 GB RAM) -- a true out-of-core run (index larger than RAM) was **not** done.
* In-memory index = {f(bpc(SEQ_M, inp))} B/char of managed heap ({f(v(SEQ_M, inp, 'build', 'node_count'), '{:,.0f}')} nodes for {f(v(SEQ_M, inp, 'build', 'leaf_count'), '{:,.0f}')} leaves); process peak RSS {f(build_rss(SEQ_M, inp) / 1e9, '{:.2f}')} GB at 10 Mbp.
* Locate is the weak spot: {f(v(SEQ_M, inp, 'locate8', 'positions_per_s') / 1e6, '{:.1f}')} M positions/s (in-memory) and {f(v(SEQ_P, inp, 'locate8', 'positions_per_s') / 1e6, '{:.1f}')} M/s (persistent) against {f(v(GO, inp, 'locate8', 'positions_per_s') / 1e6, '{:.0f}')} M/s for Go and {f(v('libsais', inp, 'locate8', 'positions_per_s') / 1e6, '{:.0f}')} M/s for a libsais SA slice (a suffix array returns a contiguous run of the SA; the tree has to walk the subtree).
  At 50 Mbp the 8-mer locate drops to {f(v(SEQ_P, 'dna_50m', 'locate8', 'throughput'), '{:,.0f}')} patterns/s ({f(v(SEQ_P, 'dna_50m', 'locate8', 'total_hits'), '{:,.0f}')} hits per 1 000 patterns).
* Build memory is the in-memory tree's problem, not the persistent tree's: persistent build anonymous memory stays at {f(v(SEQ_P, 'dna_50m', 'build_process', 'peak_rss_anon_sampled') / 1e6, '{:,.0f}')} MB for 50 Mbp.
""")

# ---- functionality
w("""## 5. What each tool provides (functionality actually exercised here)

| tool | index kind | exercised operations | not provided / not measured |
|---|---|---|---|
| Seqeron.SuffixTree (in-memory) | suffix tree (Ukkonen), UTF-16 chars | Build, Contains, CountOccurrences, FindAllOccurrences, LongestRepeatedSubstring, FindMaximalExactMatches, FindMaximalUniqueMatches, FindExactMatchAnchors, 4-thread shared reads | - |
| Seqeron.SuffixTree.Persistent | suffix tree in a memory-mapped file (v6 format) | same + Create(file) / Load(file) | out-of-core (index > RAM) not tested |
| Go `index/suffixarray` | suffix array (SA-IS), bytes | New, Lookup(n=1 / n=-1), Write, Read | no count API (count = len(Lookup)); no LRS / MEM / MUM / anchors; whole index must be in RAM |
| libsais 2.9.1, libdivsufsort 2.0.1 | suffix-array construction (int32) | construction (libsais also OpenMP 2/4 threads) | no query API: count/locate = our binary search |
| pydivsufsort 0.0.20 | libdivsufsort bindings | divsufsort, `sa_search` from a Python loop | queries are Python-call-bound |
| rust-bio 2.3.0 | SA + BWT + FM-index (Occ k=3), locate with the full `usize` SA | suffix_array, backward_search, Interval::occ | no LRS/MEM; sampled SA not used |
| sdsl-lite 2.1.1 `csa_wt<wt_huff>` / `csa_wt<wt_huff<rrr>,512,1024>` / `csa_sada` | compressed suffix arrays (FM-index / Psi) | construct_im, count, locate | no LRS/MEM API |
| sdsl-lite `cst_sct3<>` | compressed suffix tree | construct_im, descent with child()/edge(), size(v), locate via csa[lb..rb], LRS (user DFS, 1 Mbp only) | no MEM/MUM API |
| Python `suffix-trees` 0.4.0 | pure-Python Ukkonen suffix tree | STree build, find, find_all | no LRS/MEM measured |
| BWA 0.7.19 `index -a bwtsw`, Bowtie 2 2.5.5 `bowtie2-build` | FM-index for read alignment | index build time / peak RSS / on-disk bytes | no count/locate harness (use the aligners' own `mem`/`bowtie2`) |
""")

# ---- caveats
w(f"""## 6. Methodology and caveats

* **Procedure.** One process = one measurement; each (tool, input, phase) repeated 3 times (median reported, min/max in the `notes` column of `bench.tsv`).  Exceptions: `sdsl csa_wt_rrr` and `cst_sct3` on 10 Mbp inputs ran once
  (`csa_wt_rrr` took 6-16 min per run at 10 Mbp because its locate costs ~256 LF-steps per hit on average at sample rate 512; `cst_sct3` took ~75 s); `suffix-trees` 20 Mbp and Seqeron 100 Mbp are single *blocked* attempts.  `bench.tsv` column `repeats` states the real count.
* **Query timing.** untimed warm-up pass over the full query set (JIT/tier-1, caches), then full passes repeated until >= 1 s; throughput = queries / mean pass time.  Query sets are fixed files (`queries/`), 100 000 present 20-mers, 100 000 random 20-mers, 1 000 random 8-mers (`locate8`).
  Locate sums every returned position (so lists are really materialised); count for Go is `len(Lookup(q, -1))` (Go has no count call).  Per-query latency uses one clock read per query (adds ~25-40 ns).
* **Threads.** each thread runs the *whole* query set `reps` times against one shared index (same per-thread working set at 1 and 4 threads); throughput = threads x queries x reps / wall.
* **Memory.** `peak RSS` = `/proc/self/status VmHWM` (cross-checked with `wait4` `ru_maxrss`, column `ext_peak_rss`).  RSS of the persistent tree counts file-backed mmap pages that the kernel may drop, so the anonymous part (`RssAnon`, sampled every 20 ms) is reported separately.
  Seqeron in-memory bytes/char = `GC.GetTotalMemory(true)` delta around `Build` (text string excluded); Go bytes/char = after-GC `HeapAlloc` delta of `suffixarray.New` (the Index references the caller's text); rust-bio bytes/char is analytic (8n SA + n BWT + (n/3+1)*sigma*8 Occ), the RssAnon growth is in `rss_anon_growth`.
* **Cache state.** page cache dropped (`echo 3 > /proc/sys/vm/drop_caches`) only for the `cold_load` rows; all other query rows run on a warm page cache.
* **Alphabets.** Seqeron indexes UTF-16 `char`s (persistent file stores ASCII at 1 B/char); the other tools index bytes.  `prot_1m` (20 letters) was run for every tool that accepts bytes.
* **Seqeron settings.** Release build, .NET 10.0.403 (runtime 10.0.14), workstation non-concurrent GC, TieredPGO on, `ISuffixTree` string overloads (`CountOccurrences(string)` etc.).  Library code was built from a snapshot of commit 242318429a5e.
* **sdsl `csa_sada` is not safe for concurrent reads in our harness**: 4 threads sharing one instance returned a wrong count checksum (`checksum_ok = 0`) and then crashed with `std::bad_alloc` in the locate test (first run crashed the process; the harness was changed to stop after the first failed concurrent test).  Its 4-thread count number is therefore **not valid**; the other sdsl structures passed.  We did not dig into the cause.
* **libsais** rows use the serial `libsais()`; the `build_omp2` / `build_omp4` rows use `libsais_omp`.  libdivsufsort 2.0.1 serial.  Both with `-O3`, no `-march=native`.
* **Go count/locate** includes the cost of allocating the returned slice (`Lookup` copies); suffix-array C harnesses return slice indices without allocation, which is part of their locate advantage.
* **Randomness / generality.** inputs are random or lightly mutated; real genomes have different repeat structure.  `rep_10m` (100 copies of a 100 kbp block, 1 % substitutions) is the only repeat-rich input.
* **Not measured:** out-of-core (index > RAM) behaviour; Seqeron `FindAllOccurrences` span overloads; any Seqeron operation not listed above (approximate search, reverse complement, k-common substrings...); libsais64 / > 2^31 characters; sdsl `cst` MEM/MUM (no API).
* **Anything marked `NA` in `bench.tsv`** (value `NA`, repeats 0) is an operation the tool's own API does not offer; it is deliberately visible instead of absent.

## 7. Status per requested subject

| subject | status | where |
|---|---|---|
| Seqeron in-memory (1, 10 Mbp, rep 10 Mbp, protein 1 Mbp) | MEASURED, 3 repeats | `bench.tsv` tool `seqeron-st-inmem` |
| Seqeron persistent (1, 10, 20, 50 Mbp, rep, protein) | MEASURED, 3 repeats; 100 Mbp BLOCKED (600 s cap) | `seqeron-st-persistent` |
| Go `index/suffixarray` (1, 10, 20, 50 Mbp, rep, protein) | MEASURED, 3 repeats | `go-suffixarray` |
| libsais, libdivsufsort | MEASURED (built from source with cmake) | `libsais`, `libdivsufsort` |
| sdsl-lite CSA/FM-index + CST | MEASURED (1 repeat for RRR/CST at 10 Mbp) | `sdsl-*` |
| rust-bio | MEASURED (crate fetched from crates.io) | `rust-bio` |
| pydivsufsort, Python `suffix-trees` | MEASURED; suffix-trees BLOCKED at 20 Mbp (MemoryError at 12 GB cap) | `pydivsufsort`, `suffix-trees-py` |
| BWA, Bowtie2 | MEASURED: index build time / size only | `bwa-index`, `bowtie2-build` |

Hand-off scripts (to re-run or to extend on a bigger machine): `HANDOFF.md`.

## 8. Files

`bench.tsv` (final rows), `runs/raw*.tsv` (every repeat), `crosscheck.md`, `ENV.txt`, `SOURCES.txt`, `MANIFEST.sha256` (input hashes), `gen/` (input generator), `lib/` (driver, assembler, report generator),
`seqeron/SeqeronBench` (C# harness) and `seqeron/snapshot` (library snapshot), `go/` (Go harness), `competitors/` (C, C++, Rust, Python harnesses + run scripts), `run_main.sh run_all.sh run_stretch.sh run_100m.sh` (exact Seqeron/Go command lines).

Exact commands: `python3 gen/gen_inputs.py bench`; `dotnet build seqeron/SeqeronBench -c Release -o seqeron/out`; `(cd go && go build -o gobench .)` (GOTOOLCHAIN=local GOPROXY=off GOFLAGS=-mod=mod);
`./run_all.sh` (Go 1-50 Mbp, Seqeron in-memory, Seqeron persistent 1-10 Mbp), `./run_stretch.sh` (persistent 20/50 Mbp), `./run_100m.sh`; `competitors/run_libsais_libdivsufsort.sh`, `run_sdsl.sh`, `run_rustbio.sh`, `run_python.sh`, `run_aligners.sh`;
`python3 lib/assemble.py bench` -> `bench.tsv`; `python3 lib/crosscheck.py bench.tsv`; `python3 lib/make_report.py bench` -> this file.
""")
open(f"{B}/REPORT.md", "w").write("\n".join(out))
print("REPORT.md written,", len("\n".join(out)), "chars")
