# SuffixTree.PerfHarness

Repeatable macro performance harness for the suffix-tree library (ST-S0 work package 0.7). It rebuilds the **10 Mbp random-DNA row of
[`docs/refactoring/suffix-tree-upgrade/BENCHMARK.md`](../../docs/refactoring/suffix-tree-upgrade/BENCHMARK.md)** against the *current* sources
(project references, not a snapshot) and flags every metric that leaves a +-15 % band, so later stages (S1 ... S8) can prove "faster" or
"not slower" with the same numbers.

It is an application (`apps/`), not a test project: `dotnet test` never runs it, and it is not wired into CI.
It complements `apps/SuffixTree.Benchmarks` (BenchmarkDotNet micro-benchmarks on 10-100 kbp texts); this one measures a 10 Mbp build,
memory and query throughput end to end, in fresh child processes, with the line protocol of the harness that produced the published numbers
(`analysis/bench/seqeron/SeqeronBench`).

## Run

1. Generate the inputs **outside the repository** (deterministic, about 4 s; never commit them):

   ```bash
   python3 docs/refactoring/suffix-tree-upgrade/analysis/bench/gen/gen_inputs.py /tmp/stbench
   (cd /tmp/stbench && sha256sum -c /path/to/docs/refactoring/suffix-tree-upgrade/analysis/bench/MANIFEST.sha256)   # every file must say OK
   ```

2. Close everything else (timing runs only on a quiet machine; the harness prints the load average) and run:

   ```bash
   dotnet run -c Release --project apps/SuffixTree.PerfHarness -- run --data /tmp/stbench
   # about 3.5 minutes: 3 repeats x (in-memory build+queries, persistent build, persistent queries, load)
   ```

   Options: `--repeats 3`, `--tolerance 0.15`, `--targets inmem,persistent,load`, `--out result.tsv`, `--work <dir for the tree file>`,
   `--keep-trees`, `--report-only` (always exit 0). Exit code: 0 = every metric inside the band or better and every oracle equal,
   1 = regression / wrong query result, 2 = usage error.

Each measurement runs in its own child process (memory peaks need a fresh process); the median of the repeats is compared with the baseline.

## What is checked (input `dna_10m`)

| Subject | Metrics compared with BENCHMARK.md (band +-15 %) |
|---|---|
| in-memory tree | build time, heap bytes per character, peak RSS after build, count throughput (20-mers present), count latency p50, count throughput on 4 threads, locate throughput |
| persistent tree (mmap) | build time, file bytes per character, peak RSS of the build process, count throughput, count latency p50, count throughput on 4 threads, locate throughput, load time (ceiling x5: it is a few milliseconds) |
| oracles (exact equality) | node and leaf counts, sums of the count and locate results, LRS length, number of MEMs / MUMs / anchors of the 100 kbp query and the MEM checksum |

Verdicts: `OK`, `IMPROVED` (better than the band: informational, update the baseline when the improvement is intended), `REGRESSION`
(worse than the band), `MISMATCH` (an oracle differs: the library returns different answers - stop and fix), `MISSING`.
Metrics without a baseline (strict-validation time, p99 latency, ...) are printed as `info`.
Other inputs (`dna_1m`, `rep_10m`, ...) run too but are only reported, because only the `dna_10m` baseline is recorded (`Baseline.cs`,
copied from `analysis/bench/bench.tsv`).

## Reference result (ST-S0, 2026-11-03)

Library at the end of S0 work package 0.6, inputs regenerated and verified against `MANIFEST.sha256` (36 of 36 files OK), 4 vCPU, 3 repeats,
quiet machine: 39 of 39 checked metrics inside the band or better, all 19 oracle values equal. Ratios measured/baseline: in-memory build 0.81
(faster), bytes per character 1.00 / 1.00, peak RSS 1.00 / 1.00, count throughput 0.97 / 1.06, count p50 0.90 / 0.89, 4-thread count 1.08 / 1.16,
locate throughput 0.96 / 1.11, persistent build 0.89, load 1.67 (6.6 ms against 3.9 ms; a fresh-process first load including JIT; the pre-S0
code measures the same 5.8 ms on this machine in an A/B run, so it is not a regression). Full table and the old-versus-new A/B in
`docs/refactoring/suffix-tree-upgrade/PROGRESS-S0.md`.

## Updating the baseline

`Baseline.cs` holds the numbers of BENCHMARK.md (median of 3, 4 vCPU Intel Xeon under KVM, .NET 10.0.403, library commit `242318429a5e`).
Absolute values depend on the machine: on different hardware first re-measure the reference commit, then compare. Change the baseline only
together with a recorded measurement (`PROGRESS-*.md`) and a reason.
