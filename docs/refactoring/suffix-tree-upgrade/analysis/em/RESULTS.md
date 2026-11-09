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

## E2: the difference-cover sample (2026-11-09)

What is added to E1: (a) a difference cover D of Z_X, **constructed and proved by test** (below); (b) the sample ranks, computed in RAM (`DcsBuilder`): the sample positions (`p mod X` in D, p in [0, n]) are sorted by
their X-prefix and named, equal names get the same number; when all names are distinct they are the ranks (`dna_10m`, `rep_10m`, `prot_1m`), otherwise the names of each residue class form a reduced string whose suffix
array (SA-IS `SaisReference<int>` up to 16 M symbols, this builder recursively above) gives the ranks; (c) the comparator of Karkkainen (TCS 2007): for two suffixes i, j find l < X such that i+l and j+l are both sample
positions (table `Witness[(j-i) mod X]`), compare the first l symbols directly, then compare the two sample ranks. Nothing else changes. The 32-bit fast division (Lemire) replaces five 64-bit divisions per call; the merge sort replaces the
introsort for key ties (both are part of the final code; see "Ablation" below).

**Difference covers (the ADR said the table was unverified).** There is no table: `DifferenceCover.Create(X)` constructs a cover for any X. For `X = q^2 + q + 1` with q prime it builds Singer's perfect difference set (q + 1 residues, every non-zero
difference exactly once; the residues i with a zero constant coefficient of alpha^i in GF(q^3), alpha a primitive root found by search); for any other X it prunes the square-root construction `{0..m-1} + {m, 2m, .., m*m}` greedily (random order, fixed seed).
Proof by test (`DifferenceCoverTests`, default run, a few seconds): for **every X from 2 to 2300** an independent exhaustive loop checks that all X residue differences occur among the pairs of D, and that the witness table `Witness[d]` names an
a in D with (a + d) mod X in D for every d; for the 15 Singer moduli q = 2..47 the set is perfect and has q + 1 elements. Sizes: exact periods 512 / 1024 / 2048 / 4096 give |D| = 33 / 47 / 67 / 95 (6.5 / 4.6 / 3.3 / 2.3 % of the positions; the ADR guess
was 0.04-0.05 for X = 1024); Singer periods 553 / 993 / 1893 give |D| = 24 / 32 / 44 (4.3 / 3.2 / 2.3 %). `NearestSingerPeriod` maps the nominal 512 / 1024 / 2048 of the ADR to 553 / 993 / 1893.

**Round 1 (the design as written: nominal X = 1024 = Singer 993, K = n/1000, 1 thread).** `analysis/em/e2_r1.tsv`; all 38 suffix arrays equal to the libsais hash or to `SaisKernel`. `total` = rank table construction + bucket build; ratio = total / total of `dna_10m` (1.052 s).

| family | n | runs | rank table (s) | bucket build `sa_s` (s) | total (s) [min-max] | ratio to `dna_10m` | gate 10x |
|---|---|---|---|---|---|---|---|
| dna_10m | 10 M | 5 | 0.073 | 0.977 | 1.052 [1.048-1.084] | 1 | - |
| rep_10m | 10 M | 3 | 0.208 | 2.766 | 2.974 [2.951-3.103] | 2.8 | pass |
| prot_1m | 1 M | 3 | 0.016 | 0.233 | 0.249 | 0.24 | pass |
| Thue-Morse | 1 M | 3 | 0.221 | 0.878 | 1.094 [0.978-1.126] | 1.0 | pass |
| Fibonacci | 1 M | 3 | 0.232 | 0.754 | 0.987 [0.961-1.020] | 0.9 | pass |
| genome-like | 10 M | 3 | 0.234 | 1.291 | 1.529 [1.518-1.551] | 1.5 | pass |
| a^n | 1 M | 3 | 0.228 | 0.648 | 0.877 [0.875-0.886] | 0.8 | pass |
| 171-bp tandem x 10^4 | 1.71 M | 3 | 0.245 | 0.953 | 1.199 [1.196-1.228] | 1.1 | pass |
| N run 1 M in 10 M | 10 M | 3 | 0.260 | 1.647 | 1.907 [1.872-1.945] | 1.8 | pass |
| a^n | 10 M | 3 | 0.392 | 5.919 | 6.266 [6.119-6.532] | **6.0** | pass |
| Fibonacci | 10 M | 3 | 0.414 | 7.895 | 8.285 [8.275-8.960] | **7.9** | pass (largest) |
| Thue-Morse | 10 M | 3 | 0.388 | 6.641 | 7.038 [6.973-7.095] | 6.7 | pass |

(E1 for comparison: a^n 1 M 118 s, tandem 505 s, N run 412 s, a^n 10 M above 600 s.)

**Gates.** (1) Every adversarial family <= 10x `dna_10m`: **met** (largest: Fibonacci 10 M, 7.9x; a^n 10 M 6.0x; the three that failed E1 by 117-502x are 0.8-1.8x). (2) `dna_10m` slowdown caused by the DCS < 25 %: **met**:
paired interleaved A/B, 9 fresh-process pairs (`e2_ab_dcs.txt`), rank table construction included: median ratio 1.074 (7.4 %; pairs 1.035-1.155; the A/A of the same session is 1.006 [0.946-1.077]). No tuning round was needed (the ADR allows two before declaring "use F1").
`dna_10m` with 4 threads: 0.606 s total; a^n 10 M 3.45 s; Fibonacci 10 M 3.23 s.

**Sweep (informational; 3 runs for dna_10m / a^n / Fibonacci, 1 run for rep_10m, genome, tandem, N run; `e2_sweep.tsv`).** Total seconds (ratio to `dna_10m` of the same configuration):

| config | dna_10m | rep_10m | genome-like 10 M | a^n 10 M | Fibonacci 10 M | tandem 1.71 M | N run 10 M |
|---|---|---|---|---|---|---|---|
| X=512 exact (|D|=33) | 1.09 | 3.36 (3.1x) | 1.35 (1.2x) | 4.09 (3.8x) | 7.02 (6.5x) | 0.97 (0.9x) | 1.65 (1.5x) |
| X=553 Singer (24) | 1.06 | 3.43 (3.3x) | 1.44 (1.4x) | 5.17 (4.9x) | 7.48 (7.1x) | 1.09 (1.0x) | 1.89 (1.8x) |
| X=993 Singer (32), K=n/1000 (round 1) | 1.05 | 2.97 (2.8x) | 1.53 (1.5x) | 6.27 (6.0x) | 8.28 (7.9x) | 1.20 (1.1x) | 1.91 (1.8x) |
| X=1024 exact (47) | 1.09 | 2.91 (2.7x) | 1.29 (1.2x) | 4.49 (4.1x) | 8.32 (7.6x) | 1.02 (0.9x) | 1.74 (1.6x) |
| X=1893 Singer (44) | 1.04 | 2.81 (2.7x) | 1.36 (1.3x) | 7.16 (6.9x) | 8.95 (8.6x) | 1.42 (1.4x) | 2.11 (2.0x) |
| X=2048 exact (67) | 1.09 | 2.97 (2.7x) | 1.36 (1.2x) | 6.29 (5.7x) | 9.31 (8.5x) | 1.28 (1.2x) | 2.07 (1.9x) |
| X=993, K=n/300 | 0.93 | 3.07 (3.3x) | - | 5.80 (6.2x) | 7.59 (8.1x) | 1.11 (1.2x) | 1.76 (1.9x) |
| X=993, K=n/3000 | 1.24 | 3.24 (2.6x) | - | 6.32 (5.1x) | 8.29 (6.7x) | 1.32 (1.1x) | 2.22 (1.8x) |
| X=993, K=n/1000, 4 threads | 0.61 | 1.35 (2.2x) | - | 3.45 (5.7x) | 3.23 (5.3x) | - | - |

Reading: every configuration passes both gates (largest ratio 9.3 / 1.09 = 8.5x for Fibonacci 10 M with the exact period 2048). A smaller period is faster on the repetitive texts (a^n 10 M: 4.1 s at X = 512, 6.3 s at 2048) because
the direct part of the comparator is shorter (average l is about X/2), at the price of a denser sample (6.5 % of the positions at 512, 2.3 % at 1893: rank file 0.26 B/char against 0.09 B/char). On random DNA the choice of X is irrelevant (1.04-1.09 s).
Fibonacci is the worst family at every X (many ties, long LCP, pairs of near-neighbours in the order).

**Ablation of the implementation changes (single unrepeated runs of a^n 10 M at X = 993, made while another process was running on a different core, so only the order of magnitude counts):** first version (introsort, 64-bit `%` and `/`) 11.1 s;
with the 32-bit fast division 9.2 s; with the merge sort for tie runs 5.8 s. The E1 numbers above are the first version (introsort); `e1_merge.tsv` has the merge-sort variant.

**What E2 does not show.** The rank table is built in RAM by safe code (`DcsBuilder`, `SaisReference` for the reduced string): its construction time (0.07-0.4 s per 10 M) and memory (about 60 B per sample position transient) are not those of the
external PH1 of the design; n is at most 10 M, so the rank table is 1.3 MB-2.6 MB and fits in the CPU cache: the cost of rank lookups that miss the cache is what E3 measures; a text of 3 x 10^9 symbols was not tried.
