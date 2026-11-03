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
