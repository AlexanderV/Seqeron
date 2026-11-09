# SuffixTree.Esa.Bench

Measurement harness of the ESA backend (ST-S1 work package 1a.8): SA + LCP build time, peak RSS and index bytes of the `SaisReference` /
`LcpReference` pipeline (from 1o.1 on: of the optimized kernels), plus a phase profile from stopwatch counters inside the reference kernels
(`SaisTrace`, `LcpTrace`; no external profiler). It is an application, not a test project: `dotnet test` never runs it.

Protocol and gates: [`PROGRESS-S1.md`](../../docs/refactoring/suffix-tree-upgrade/PROGRESS-S1.md), "Measurement protocol for the gate". One process = one
measurement (JIT warm-up included, peak RSS needs a fresh process); the input is loaded into memory first and is not part of any timer.

```
dotnet build apps/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.csproj -c Release
python3 docs/refactoring/suffix-tree-upgrade/analysis/bench/gen/gen_inputs.py /tmp/stbench        # outside the repository
docs/refactoring/suffix-tree-upgrade/analysis/bench/run_esa.sh build   /tmp/stbench dna_10m 5      # 5 fresh processes
docs/refactoring/suffix-tree-upgrade/analysis/bench/run_esa.sh steady  /tmp/stbench dna_10m 3      # 5 builds per process, last 3 reported
docs/refactoring/suffix-tree-upgrade/analysis/bench/run_esa.sh profile /tmp/stbench dna_10m 3
docs/refactoring/suffix-tree-upgrade/analysis/bench/run_esa.sh aggregate                          # -> analysis/bench/bench_esa.tsv
```

`build` prints `build_prepare_s`, `build_sais_s`, `build_sa_s` (prepare + SA-IS), `build_lcp_s`, `build_total_s` (the gate metric), `peak_rss_hwm_after_build`
(`/proc/self/status` VmHWM; the driver adds the independent `ext_peak_rss` from `wait4`), GC allocation and collections, `index_bytes_*` (SA 4 B x (n+1),
compact LCP payload, remapped text), and the SHA-256 of SA[1..n] and LCP[1..n]; with `--verify analysis/sa_vectors.tsv` it compares them with the libsais
vectors (`vector_match`, exit code 3 on a mismatch). `build_lcp_jit_s` / `build_lcp_jit_methods` (1o.3c) are the JIT compilation time and method count of the LCP step on the building thread: in a cold process they are part of `build_lcp_s` (3.7 ms of the 6.8 ms at 1 Mbp), and a kernel that adds code to the common path shows it there first. Rows go to `analysis/bench/bench_esa.tsv` in the format of `bench.tsv`; `bench.tsv` is never modified.

## Child table and query throughput (1b.1, 1b.2)

`query` builds SA + LCP once (SA-IS pointer kernel, direct LCP kernel) and then, per child-table layout, times the child-table scan on its own, reports the table bytes, optionally the
write trace of the scan (`--trace true`: writes at the cursor / behind it, backward runs, page misses of an LRU cache of 16 / 64 / 1,024 pages of 4,096 bytes; prototype P3, ADR-004
section 3.2) and the query throughput of the in-memory ESA with that layout: 20-mer count of `present20` and `absent20`, locate of the 8-mers (`FindAllOccurrences`, which allocates the
result, and the lazy `TryMatch` + `Locate` into a reused buffer). Timing as in `BENCHMARK.md`: one untimed warm-up pass over the query file, then full passes until `--min-seconds`
(default 1.0) have elapsed, throughput = queries / mean pass time. One process = one text; repeat in fresh processes for medians. Without `--layout` the engine default for the
text length is used (`EsaIndex.DefaultLayoutFor`).

```
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll query --input /tmp/stbench/data/dna_10m.txt --name dna_10m --queries /tmp/stbench --min-seconds 2
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll query --input /tmp/stbench/data/rep_10m.txt --name rep_10m --queries /tmp/stbench --layout all --trace true
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll query --synthetic fibonacci:1000000 --layout all --trace true    # a^n, abab, fibonacci, thue-morse, n-runs
```

`--queries` is the bench directory that holds `queries/<name>.present20.txt` etc. (`gen_inputs.py`); a missing query file skips that row. Do not time anything while a build or a test run is active.

## Native analysis operations (1b.3)

`analysis` builds the whole in-memory ESA once (`EsaIndex.Build`: prepare, SA-IS, LCP, child table, statistics) and then times each analysis operation on that index: the first call (cold JIT)
and the median of three further calls, plus the bytes allocated by the first call: `LongestRepeatedSubstring` (and the zero-copy `...Memory`), `FindAllLongestRepeatedSubstrings` (ties, positions,
output characters), `CountDistinctSubstrings`, `CountDistinctSubstringsByLength` for 100, 10,000 and n, `EnumerateSuffixes().Take(20)` (the cost is the copy of the strings: ns per character),
`GetAllSuffixes` of the first 2,000 and 10,000 characters (Theta(n^2) output by definition). `--generic true` adds the operations that run on `Traverse` through the generic algorithms of
`SuffixTree.Core`: `FindMaximalRepeatedPairs(--min-length)` and `FindLongestCommonSubstrings` of the two halves of the text with support 2. `--check true` computes the same results through the
`Traverse`-default path of the same index (`SuffixTreeAlgorithms.CountDistinctSubstrings`, `...ByLength`, `FindAllLongestRepeatedSubstrings`; it also reports their times) and compares them with the native
results; a difference sets exit code 3. One process = one text; nothing else may run while it is timed.

```
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll analysis --input /tmp/stbench/data/dna_10m.txt --name dna_10m --check true --generic true --min-length 16
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll analysis --input /tmp/stbench/data/rep_10m.txt --name rep_10m --check true --generic true --min-length 500
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll analysis --synthetic fibonacci:1000000 --check true
```

## Suffix-link prototypes P1 and P2 (1b.4a)

Prototype code, not the navigator of 1b.4: it lives in `Links/` (compiled into this app and, by a `Compile Include`, into `SuffixTree.Esa.Tests`, which holds its oracle tests in `Links/`).
`LinkTables` (ISA, Psi), `LcpRangeHelper` (pyramid of block minima over the compact LCP with range minimum + leftmost position, previous / next smaller value), `LinkKernels`
(four interval searches for the link `SL(v) = LCA(Psi(l), Psi(r))` of ADR-001 section 5: (a) outward LCP scan, (b) child-table descent along the rank path of `Psi(l)` -- the child table
has no parent pointers, so a parent climb is not available --, (c) range minimum + previous / next smaller of its position, (c2) previous / next smaller than the known value `depth - 1`;
each with ISA or Psi endpoints), `LinkNodeSet` (all internal nodes, node depths, seeded sample, largest keys), `MemDriver` (throwaway matching-statistics driver: extend, follow the link of
the last node, rescan; MEMs by an outward LCP scan from the longest match).

```
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll links --input /tmp/stbench/data/dna_10m.txt --name dna_10m [--sample 200000] [--worst 1000] [--fanout 16] [--kernels all|scan,descent,rmq,psv] [--sources all|isa,psi]
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll links --synthetic fibonacci:1000000
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll mem --input /tmp/stbench/data/dna_10m.txt --name dna_10m --queries /tmp/stbench [--kernel rmq] [--source psi] [--loci 2] [--ukkonen true]
docs/refactoring/suffix-tree-upgrade/analysis/s1b/run_p1.sh /tmp/stbench dna_10m 5      # 5 fresh processes -> analysis/s1b/p1_links_raw.tsv
docs/refactoring/suffix-tree-upgrade/analysis/s1b/run_p2.sh /tmp/stbench dna_10m 5      # -> analysis/s1b/p2_mem_raw.tsv
python3 docs/refactoring/suffix-tree-upgrade/analysis/s1b/aggregate.py p1 docs/refactoring/suffix-tree-upgrade/analysis/s1b/p1_links_raw.tsv
```

`links` (one process = one text): builds the index (tagged child table), then ISA, Psi and the pyramid with their build times and bytes per character, enumerates the internal nodes, draws a
seeded sample A (per-link timing) and B (batch timing) of internal nodes, and measures every (kernel, source) combination: the cold first link, a counting pass with a step budget (steps
median / p99 / max, links that exceed `--link-budget`, the sample stops at `--total-budget` steps so an unbounded kernel on a^n ends as "budget exceeded"), a round-robin warm-up of
`--warmup` seconds (the shared readers reach their final JIT tier), the per-link timing of sample A (raw and net of the calibrated `Stopwatch` pair cost) and the batch mean on sample B. The
results of all combinations are compared with the first one (links whose kernel gave up are skipped) and 300 of them with `EsaIndex.TryMatch` of the label without its first symbol (the
definition); a difference sets exit code 3. `--sweeps` first runs kernel (c2) once over every internal node (untimed) to get the target width and the exact outward-scan cost of each node,
then measures all combinations on the worst nodes by scan cost, target width, string depth and node depth. A final `repeat_*` row repeats the first combination as a drift check.

`mem` (one process = one text): builds the ESA, ISA, Psi and the pyramid, runs the benchmark's `<name>.memq100k` query (100 kbp, 1 % substitutions, minLength 20) `--calls` times (the first
call is the cold-JIT time, `warm_median_s` is over calls 2..N, `steady_median_s` over the last five), reports the counters of the driver (links followed, extension steps, child lookups, rescan
lookups) and, with `--time-links`, the share of the kernel calls; `--ukkonen true` then builds the in-memory Ukkonen tree of the same text in the same process, answers the same query (cold and
warm) and compares the two MEM sets (`mem_set_equal`). `--loci 2` follows the second chain of the generic algorithm (its suffix links are paid, its result is not used), `--loci 1` is the
single-chain matching statistics of an ESA.

## Production navigator behind the generic algorithms (1b.4)

`nav` measures the real thing: `EsaIndex.FindMaximalExactMatches` / `FindMaximalUniqueMatches` / `FindExactMatchAnchors` / the LCS family, i.e. the generic algorithms of `SuffixTree.Core`
(`SuffixTreeAlgorithms.*<EsaNode, EsaNavigator>`) on the production navigator, the navigation pack (ISA + LCP pyramid of fanout 32, built on the first matching call) and the suffix-link kernel (c2).
One process = one text and one operation; the in-memory Ukkonen tree of the same text answers the same query in the same process (reference time, answer comparison, exit code 3 on a difference).

```
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll nav --input /tmp/stbench/data/dna_10m.txt --name dna_10m --queries /tmp/stbench --op mem|mum|anchors [--min-length 20] [--calls 10] [--ukkonen true]
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll nav --input /tmp/stbench/data/dna_1m.txt --name dna_1m --op lcs [--other-kind mutated|random | --other <file>]
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll nav --input /tmp/stbench/data/dna_10m.txt --name dna_10m --op links [--sample 200000]
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll nav --input /tmp/stbench/data/dna_10m.txt --name dna_10m --queries /tmp/stbench --op split
docs/refactoring/suffix-tree-upgrade/analysis/s1b/run_nav.sh /tmp/stbench dna_10m mem 5      # 5 fresh processes -> analysis/s1b/nav_raw.tsv
```

`mem` / `mum` / `anchors`: the benchmark's `<name>.memq100k` query (100 kbp, 1 % substitutions), minLength 20. Rows: ESA build, the pack build on its own (`pack_build_s`, `pack_bytes_per_char`, a second warm pyramid
build), then the operation `--calls` times: `*_first_call_s` is the cold-JIT time (the pack already exists), `*_warm_median_s` the median over calls 2..N, `*_steady_median_s` over the last five, result count and
checksum; the same for the Ukkonen tree (`nav_ukkonen`) and `nav_check` (answers equal). `lcs`: `LongestCommonSubstringInfo`, `FindAllLongestCommonSubstrings`, `FindAllDistinctLongestCommonSubstrings` against a second
text (a deterministic copy of the text with 1 % substitutions and 0.1 % single-letter indels, `mutated`, or an independent random text of the same letters, `random`, or a file). `links`: the production
kernel and the prototype kernel (c2) of 1b.4a (ISA, fanout 32, WITH the extra `SA[Lb]` read that completes the production handle) on the same sampled internal nodes in paired rounds (two interleaved halves of the
sample, the halves swap and the order alternates), ns per link net of the timer floor. `split`: MEM through a counting navigator (exact call counts per member), through a timing navigator (each heavy member
timed one by one; the timers serialise the pipeline, so the shares of the timed run, not its absolute times, are meaningful) and through a replay (the recorded arguments of each member played back on their own
with a data dependency); the three views bracket the share of the link kernel and of the child lookups, which are not additive because the processor overlaps them.

## Paired interleaved A/B runner (1o.0)

`ab` is the acceptance instrument of the optimization phase (PROCESS.md, amendments of 2026-11-03). It spawns FRESH `build` processes of this app, so every
number has cold-JIT start-up and its own peak RSS, and runs one session: quiet-machine check (`/proc/loadavg` below `--max-load`), libsais reference runs,
an A/A calibration (the configuration A against itself, `--calibration-pairs`, default 7), the A/B comparison (`--pairs`, default 9, at least 7; use 15 or more for
a decision near the bar), libsais again. Inside a pair the order alternates (A B, B A, ...); the statistic is the median of the pair ratios B/A. A and B differ by
`--a-args` / `--b-args` (build options, for example `--sa-kernel pointer`), `--a-env` / `--b-env` (runtime switches such as `DOTNET_TieredPGO=0`) or `--a-dll` /
`--b-dll` (two builds of the app, for example from two git worktrees). A child that does not reproduce the libsais vectors (`--verify analysis/sa_vectors.tsv`)
aborts the session. The libsais driver is outside the repository (`--libsais <driver> ` is called as `<driver> <input file>` and must print `total=<seconds>`).

```
dotnet build apps/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.csproj -c Release
dotnet apps/SuffixTree.Esa.Bench/bin/Release/net10.0/SuffixTree.Esa.Bench.dll ab --input <dir>/data/dna_10m.txt --name dna_10m \
    --verify docs/refactoring/suffix-tree-upgrade/analysis/sa_vectors.tsv --a-args "--sa-kernel reference" --b-args "--sa-kernel pointer" \
    --pairs 15 --calibration-pairs 9 --libsais <libsais driver> --out raw.tsv
```

Statistics (`PairedStats.cs`, pinned by `Bench/PairedStatsTests.cs`): noise band `[1 - h, 1 + h]`, `h = |median(A/A ratios) - 1| + 2 se`, `se = 1.2533 * 1.4826 * MAD / sqrt(N)`.
Verdict ACCEPT only if the median B/A ratio is at most 0.97, below the lower edge of the band, and the 95 % interval of the median (order statistics) is below 1.0;
a median within one percentage point of 0.97 is flagged BORDERLINE and the session is repeated with more pairs. The first live A/A session (S1a reference both sides, 9 pairs)
showed 0.9726 between identical binaries, which is why the flag exists.

## Checked-pointer configuration

`docs/refactoring/suffix-tree-upgrade/analysis/run_esa_checked.sh [dotnet-test args]` runs the Esa tests with the symbol `ESA_BOUNDS` (`-p:EsaBounds=true`): every `Buf<T>` access is
bounds-checked. The plain Release run is the raw-pointer configuration. Both must pass.
