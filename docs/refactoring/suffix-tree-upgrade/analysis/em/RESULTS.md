# ADR-004 prototypes E0-E3: measured results

Prototype code: `apps/SuffixTree.EmProto` (managed C#, no native code), tests `tests/SuffixTree/SuffixTree.EmProto.Tests`. Raw rows: the TSV files in this directory
(`R <TAB> experiment <TAB> case <TAB> metric <TAB> value <TAB> unit <TAB> notes`). Analysis in the ADR: [ADR-004 section 5a](../../adr/ADR-004-external-memory-construction.md).
Every number below is **measured here** (4 vCPU Intel Xeon, 15 GB, kernel 6.18.44, .NET 10.0.14, Release, workstation non-concurrent GC, nothing else running; build servers shut down)
unless it is marked estimated or literature. Inputs: `/tmp/stbench/data` (outside the repository).

## E0 harness (2026-11-09)

Command: `dotnet apps/SuffixTree.EmProto/bin/Release/net10.0/SuffixTree.EmProto.dll e0 --dir <scratch> --out e0.tsv` (raw rows: `e0.tsv`).

| gate | measured | result |
|---|---|---|
| RssAnon sampler (`/proc/self/status`, 10 ms) overhead < 2 % | 0.64 % (median sampled/unsampled ratio 1.0064 of 15 interleaved pairs of `SaisKernel` on 10 M random DNA, 95 % interval 0.994-1.052); A/A median 0.9972, A/A noise band [0.9917, 1.0083], so the overhead is inside the noise of this session | PASS |
| sampler reads the right quantity | 200 MiB (209.7 MB) of touched native memory raised the sampled RssAnon maximum by 212.1 MB | PASS |
| explicit I/O counters within 1 % of `/proc/self/io` | writes: 301,989,888 bytes counted, `wchar` delta 301,989,888 (difference 0 %); reads: 335,544,320 counted, `rchar` delta 335,702,685 (0.047 %, the difference is the sampler's own reads of `/proc/self/status`) | PASS |
| `ArenaBudget` on `NativeArena` | cap enforced (a block that would exceed it is refused and reserves nothing), `ReleaseTo` and `Dispose` return the bytes, shared by arenas of several threads (unit tests) | PASS |

Notes. The counters cover explicit `pread`/`pwrite` calls (class `CountedFile`); accesses through memory mappings are not explicit calls and are reported separately in E3 (`RssFile`, block-layer
`read_bytes` of `/proc/self/io`, `Dirty` of `/proc/meminfo`). The `ArenaBudget` is a separate class that wraps a `NativeArena` (`BudgetedArena`): `NativeArena` itself is not modified.

## E1: in-RAM simulation of PH2-PH3 without the difference-cover sample (2026-11-09)

What is simulated: the bucket partition (stratified random sample of 64 x 16 positions, sorted by the comparator, exact quantile splitters; one scan of the bucket per level, 64-bit key first and the
comparator only on key ties; a bucket larger than K is split again, fan-out 64) and the bucket sort (LSD radix on the 64-bit key, ties ordered by the comparator, SA slice written at the final
offset). The comparator is the direct suffix comparison with a vectorised LCE (AVX2, 256 bits = 128 symbols of a 2-bit text per step; scalar 64-bit fallback). Text: b-bit packed (2 bits for DNA,
3 for ACGNT, 5 for protein, 1 for a^n and the binary words). 1 thread, K = n/1000, suffix array only. `sa_s` is the time of `EmBuilder.Run` (the text is generated, loaded and packed before the timer;
packing costs 0.003-0.012 s). Every suffix array was checked: equal to the committed libsais hash (`dna_10m`, `rep_10m`, `prot_1m`, `all-same`, `fibonacci`, `thue-morse` at 10^6) or equal to `SaisKernel`
on the same text (tandem, N run, genome-like, 10 M variants): every one of the 55 completed builds checked is equal (`sa_correct = yes`: 32 in Table 1, 23 in Table 2). Raw rows: `e1.tsv` (tie runs by introsort, the code of the first measurement) and `e1_merge.tsv` (merge sort,
the final code).

Table 1: E1 as first measured (introsort for the tie runs). `runs` = fresh processes, median [min-max]; ratio = median `sa_s` divided by the median of `dna_10m` (1.006 s).

| family | n | runs | sa_s (s) | ratio to dna_10m | share of positions in key ties | LCE work (symbols matched) | comparator calls | gate 10x |
|---|---|---|---|---|---|---|---|---|
| dna_10m (gate: <= 3 s) | 10 M | 5 | 1.006 [0.973-1.041] | 1 | 0 % | 1.4 x 10^6 | 3.3 x 10^5 | **1.0 s: PASS** |
| dna_10m, 4 threads | 10 M | 3 | 0.556 [0.553-0.635] | 0.55 | 0 % | same | same | - |
| rep_10m (100 x 100 kbp, 1 % substitutions) | 10 M | 3 | 2.578 [2.430-2.652] | 2.6 | 77.5 % | 3.4 x 10^9 | 4.6 x 10^7 | pass |
| prot_1m | 1 M | 3 | 0.258 [0.235-0.277] | 0.26 | 0 % | 5.7 x 10^5 | 3.2 x 10^5 | pass |
| Thue-Morse | 1 M | 3 | 4.141 [4.031-4.273] | 4.1 | 99.99 % | 3.0 x 10^11 | 1.4 x 10^7 | pass |
| Fibonacci | 1 M | 3 | 9.414 [9.276-9.605] | 9.4 | 99.99 % | 7.3 x 10^11 | 1.5 x 10^7 | pass (marginal) |
| genome-like (random + repeats 10-15 % + 171-bp satellites + N runs) | 10 M | 3 | 2.216 [2.109-2.326] | 2.2 | 8.6 % | 2.6 x 10^10 | 4.3 x 10^6 | pass |
| a^n | 1 M | 3 | 118.0 [117.9-119.1] | 117 | 100 % | 1.0 x 10^13 | 2.1 x 10^7 | **FAIL** |
| a^n | 10 M | 1 | > 600 s (aborted by the watchdog) | > 600 | - | - | - | **FAIL** (counts as failure) |
| 171-bp tandem x 10^4 | 1.71 M | 3 | 505 [502-511] | 502 | 100 % | 2.1 x 10^13 | 2.6 x 10^7 | **FAIL** |
| N run of 1 M inside random DNA | 10 M | 3 | 411.8 [410.1-415.0] | 410 | 10 % | 1.2 x 10^13 | 2.4 x 10^7 | **FAIL** |

Table 2: sensitivity, the final code (merge sort for the tie runs; `e1_merge.tsv`). The merge sort saves 11-15 % on rep_10m / Thue-Morse / Fibonacci, 22-26 % on a^n 1 M and the genome-like text, 38-39 % on tandem and N run (fewer comparator calls, n / 2 per merge on reverse-ordered runs), and cannot save a^n 10 M:

| family | runs | sa_s (s) | ratio to dna_10m (1.044 s, same session) |
|---|---|---|---|
| dna_10m | 3 | 1.044 [0.975-1.047] | 1 |
| rep_10m | 3 | 2.293 [2.279-2.364] | 2.2 |
| prot_1m | 3 | 0.221 [0.214-0.222] | 0.21 |
| Thue-Morse 1 M | 3 | 3.575 [3.535-3.765] | 3.4 |
| Fibonacci 1 M | 3 | 7.978 [7.927-8.067] | 7.6 |
| genome-like 10 M | 3 | 1.637 [1.605-1.706] | 1.6 |
| a^n 1 M | 3 | 91.7 [90.9-95.0] | 88 |
| 171-bp tandem x 10^4 | 1 | 314.6 | 301 |
| N run 1 M in 10 M | 1 | 251.1 | 240 |
| a^n 10 M | 1 | > 600 s (aborted) | > 575 |

Noise: five runs of `dna_10m` span 0.973-1.041 s (+-3.4 % around the median), the 3-run families span 1-8 % (min to max); the slow cases (118-505 s) span 1-2 %. For scale, `SaisKernel` (in-repository,
warm, 1 thread) builds the 10 M random-DNA suffix array in 0.537 s (E0 median of 15 runs), so the simulation costs 1.9x the in-RAM kernel on random DNA (the gate is 3 s absolute).

Result against the E1 rules of ADR-004 section 5:
* Gate "SA equal to the libsais hashes": met (no mismatch ever occurred; hashes where libsais has one, `SaisKernel` otherwise, as the rules allow).
* Gate "`dna_10m` <= 3 s on 1 thread": met, 1.006 s (median of 5; 0.55 s with 4 threads).
* Rule "if every family is within 10x of `dna_10m`, the difference-cover sample is not needed": **not met**. Four families are far outside (a^n 1 M 117x, a^n 10 M above 600 s, tandem 502x, N run 410x) and Fibonacci
  at 1 M is at 9.4x. **PH1 (the difference-cover sample) is needed; E2 follows.**
* The mechanism is visible in the counters: on every failing family the total LCE work is 10^13 symbols (about 10^6 symbols per text symbol) and almost all positions are key ties (100 % for a^n and tandem; 10 % for the
  N run, which is exactly the share of the text inside the run). The number of comparator calls is the same order for all families (1.4-2.6 x 10^7 at 1-1.7 M, 2.4 x 10^7 for the N run at 10 M): the cost is the length of each call, not the count.
* Share of positions in key ties on the benign families: random DNA 0 %, protein 0 %, `rep_10m` 77 % (100 copies at 1 % divergence: a 32-mer recurs with probability 0.72), genome-like 8.6 %.

Provenance. The first measurement session (run 0 of a^n 1 M, tandem and N run in `e1.tsv`, 118 / 502 / 410 s) used the build before the scratch-file counters were added to the position store; all other rows come
from the frozen build `bin_e1` of the same source with introsort ties. The E1 code path (direct comparator) was not changed afterwards except for the added `--tie-sort` option.
