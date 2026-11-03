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
vectors (`vector_match`, exit code 3 on a mismatch). Rows go to `analysis/bench/bench_esa.tsv` in the format of `bench.tsv`; `bench.tsv` is never modified.

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
