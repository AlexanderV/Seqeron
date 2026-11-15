# RESULTS-E7: scaling of the whole pipeline and the true 0.5x page-cache test (2026-11-15)

Prototype E7 of ADR-004 section 5 (the owner authorised it on 2026-11-15): the end-to-end pipeline PH0-PH5 of E6, unchanged, run at 50 M, 100 M, 200 M and 500 M symbols in two series, plus the page-cache-limited
runs of the whole pipeline (the "true 0.5x test" of section 6b). Gates S1-S4 and C were fixed by the orchestrator before the work (task statement; quoted in the tables). Machine A5 (4 vCPU, 15 GB RAM, virtio disk, ext4,
~15 GB free), cold page cache for every timed run (`sync` + `drop_caches`), fresh process, 4 threads unless a row says otherwise, `DOTNET_GCHeapHardLimit` = 0x6000000 (96 MiB) for the scaling runs. Binary:
`dotnet publish apps/SuffixTree.EmProto -c Release -o /tmp/stbench/bin_e7_run` from commit `ba1fbc0e` (the code change of this prototype, section 2; nothing in `apps` changed after it).
Raw rows (this folder): `e7_dna.tsv` (dna 50 / 100 / 200 M), `e7_dna500.tsv` (dna 500 M), `e7_genome.tsv` (genome-like 50 / 100 / 200 / 500 M), `e7_r500.tsv` (cross-check against the E4 files), `e7_setup.tsv` (references),
`e7_cache.tsv` (cache sweep), `e7_cache_calib.tsv` (calibration of the sweep, kept). Scripts: `run_e7_setup.sh`, `run_e7_protocol.sh`, `sample_e7.py`, `run_e7_cache.sh`, `run_e7_cache_sweep.sh`, `sample_e7_cache.py`,
`e7_unpack_packed.py`, `tables_e7.py` (every table below is its output).

**No restart of the container happened** during the work. Two defects of the measuring scripts (not of the product) were found and fixed (section 9).

## 1. Result in one table

| gate | threshold | measured | result |
|---|---|---|---|
| S1 exponent, dna series A (fixed parameters), 50-500 M | b <= 1.15 (abandon > 1.3) | **b = 1.153** [1.146-1.183]; PH1 1.66, PH2/PH3 1.06, PH4+PH5 1.16 | **not met** (by 0.003 on the medians, 0.033 on the minima; between 1.15 and 1.3: not an abandon) |
| S1 exponent, dna series B (ratio emulation), 50-500 M | b <= 1.15 | **b = 0.985** [0.954-0.977]; 207 / 244 / 238 / 205 ns per symbol | **met** |
| S1 exponent, genome-like series B, 50-200 M | b <= 1.15 | **b = 1.073** [1.049-1.062] (0.987 [0.974-0.983] with the 500 M run added) | **met** |
| S2 RssAnon <= 256 MB at every n, both series | <= 256 MB | maximum over all runs **214.5 MB** (dna A, 500 M), **124.4 MB** (dna B, 500 M), 132.5 MB (genome-like B, 200 M), 146.1 MB (genome-like A, 200 M) | **met** |
| S2 series A slope of the median peak against log2 n | <= 5 MB per doubling (abandon > 15) | **36.0 MB per doubling** (dna, four sizes), 33.0 (genome-like, three sizes) | **not met; the abandon rule "series A slope > 15 MB per doubling" fires** (cause: PH4's buffers per chunk with P = n / R, section 5) |
| S2 series B slope (report) | - | 9.5 MB per doubling (dna), 8.8 (genome-like); part follows from parameters that scale with n | reported |
| S3 identity at every n: five outputs and statistics equal to the references, between series A, B and 1 / 2 / 4 threads | byte-identical | 54 timed runs (430 comparisons) + 43 finished cache runs, all equal, 0 different; 500 M CHILD against the e5-oracle reference; extra cross-check with the E4 files (pk, SA, LCP) equal | **met** |
| S4 I/O (program counters, cross-checked with rchar + wchar) | <= 100 B/char | 65.9-75.9 B/char (rchar + wchar 59.6-71.0); scratch peak 12.0 B/char; outputs 9.25 (dna) / 9.50 (genome-like) B/char | **met** |
| C dna_100m, effective cache 0.5x (resident fraction of text + rank over PH3 = 0.49 [0.46-0.51]), whole pipeline, PH0-PH3 limited | slowdown <= 3x | **3.6x [3.59-3.63]** (98.3 s against 27.3 s, PH2/PH3 5.0x); with the limit held for PH4 + PH5 too the run does not finish below an effective 0.65 (15 min) | **not met**; a true 0.5x was reached (variant D) and is 3.6x, the 3x line is at about 0.55; the abandon criterion (needs > 2x) does not fire |
| C genome-like 100 M, effective >= 1x | slowdown <= 3x (the floor of ~1x) | 0.97x [0.96-1.03] (effective 0.97, 3 runs), 0.82x (0.96), 1.23x (0.99); below: 2.8x at 0.81, 5.2-5.3x at 0.65-0.67, **8.5x [8.3-8.6] at 0.51** | **floor holds** (needs about 0.8-0.85 for 3x) |

## 2. What was added to the code

One additive option, **`--rule-n <n>`** (`EndToEndOptions.RuleN`, `Experiments/E6.cs`, one usage line in `Cli.cs`): the scale-down rule of E6 is evaluated at `RuleN` instead of the text length (`ruleN = RuleN > 0 ? RuleN : n` in
`ScaleRules.Resolve`, used for the factor `f` and for K of PH3; ignored without `--scale-n`; negative refused). With `--scale-n 3000000000 --rule-n 100000000` every text, whatever its length, runs with the values the rule gives at
10^8 symbols. All the other values of series A (the explicit options exist for them) come from the same rule, so no further plumbing was needed. **Tests first:** 8 tests (`ScaleRulesTests` +6 cases, `EndToEndPlumbingTests` +1 test and 2 assertions added to the command-line test:
the values for 50 / 100 / 200 / 500 M equal those of 10^8, default 0 = E6 rule, unscaled ignores it, an explicit value still wins, negative refused, `--rule-n` on the command line, an end-to-end run on 30,000 symbols with
`RuleN` = 6 M whose parameters and output are checked): with a stub property that nothing reads **6 red** (the case at n = 10^8 and the default-0 test pass by construction), after the implementation 32 of 32 green in the three classes
(`ScaleRulesTests`, `EndToEndPlumbingTests`, the command-line test); the whole EmProto suite was run at the end (section 8). One conditional is the whole production logic: option plumbing, no mutation audit.
**Series A reproduces E6 at 10^8:** the `parameters` row of dna_100m series A is character for character that of E6
(`n=100000000 f=3.333E-002 scaled=True: PH1 K_s=34953 N_base=559241 R=139810 K_sort=34953 run-spill=2185; PH3 K=33333 F=64; PH4 R=266667 P=375 B2 workers=2; PH5 block=262144 window=1048576 slots=8`), the five outputs and the statistics
are identical to the E6 references (hash of the E6 files) in all 3 runs, and the time (24.65 s) and RssAnon (96.7 MB) are those of E6 within its spread (24.3-24.7 s, 88.7-96 MB).

## 3. Method

* **Inputs.** Random DNA `dna:N` (seed 20240001) at 50 / 100 / 200 / 500 M, genome-like `genome:N` (seed 20240002; ACGNT, repeats, satellites, N runs) at 50 / 100 / 200 / 500 M (the 500 M genome-like text was run: disk and time allowed it),
  as ASCII files of one byte per symbol (`e6-mktext`; the 100 M texts and references are E6's).
* **Series A** (fixed parameters): `--scale-n 3000000000 --rule-n 100000000`. **Series B** (ratio emulation): `--scale-n 3000000000` (the E6 rule, f = n / 3x10^9). Median of 3 runs per cell (dna A and B at four sizes,
  genome-like B at four sizes); genome-like A 1 run per size (50 / 100 / 200 M, not run at 500 M); 1 and 2 threads once per size and input in series B (dna also 500 M).
  Runs of one size were interleaved A, B, A, B, ... (drift of the host shows in both).
* **Per run:** `sample_e7.py` (independent, 20 ms: RssAnon of `/proc/<pid>/status`, system Dirty, scratch directory size every 100 ms, `/proc/<pid>/io`), the program's own per-phase rows (10 ms RssAnon sampler, byte counters),
  the exit code (0 = RssAnon <= 256 MB, scratch gone, clean), then, outside the timed region, SHA-256 of the five outputs against the hashes of the references, the statistics, the scratch directory, and the outputs are deleted.
* **References (S3), made once per input by the in-RAM kernels of the earlier stages, as `run_e6_setup.sh`:** `pack` (packed text), `oracle --out` (SaisKernel suffix array in RAM), `lcp-oracle` (LcpKernel), `e5-oracle` (the product's
  `ChildTable.Build` in RAM and the statistics). Their SHA-256 are in `<input>.refs.sha` (not committed, in the work directory `/tmp/stbench/e7`); at 200 M and 500 M the big files were deleted after hashing (the disk is
  tight at 500 M: the SA reference was deleted as soon as the LCP files existed, `E7_DROP_SA`). The 100 M references are the files of E6 (identical to the orchestrator's E4 / E5 files, RESULTS-E6.md). Every run is compared by hash
  with the same reference, so identity between series A, series B and 1 / 2 / 4 threads follows from equality to one reference.
* **Extra cross-check (not trusted blindly):** `/tmp/stbench/orch4/r500.*` is a `packrand` text (seed 20240009, `r500.pk.meta`: 500,000,000 symbols, 2 bits, letters ACGT, a different text from `dna:500000000`) with a suffix array
  built by the E3 builder (`r500.sa`, rows `orch.r500.sa`: K = n / 3000, 4 threads) and an LcpKernel LCP (`r500.oracle.lcp`). Its packed form is MSB first, codes 0..3 = A, C, G, T, so `e7_unpack_packed.py` turns it into
  the ASCII text `r500.txt`, and the pipeline's own packed text must equal `r500.pk`. Nothing under `orch4` was modified (only read and hashed).

## 4. S1: the exponent (gate: b <= 1.15 for dna A, dna B, genome-like B; abandon above 1.3)

Total wall time of PH0-PH5 (the sum of the phase times, without the final sync), median [min-max] of the runs per size, 4 threads, cold cache; `b` by least squares of `log t` on `log n` over the medians (and over the per-size
minima and maxima); the per-phase exponents are fitted the same way. `PH2/PH3` is the bucket sort (the SA), `PH4+PH5` the fused LCP and child scan. Fits over four sizes (50 / 100 / 200 / 500 M) except where a table says three.


#### dna (four sizes)

### dna series A (4 threads, cold cache, runs per size: 3, 3, 3, 3)

| size | n | total s | PH0 s | PH1 s | PH2/PH3 s | PH4+PH5 s | ns / symbol | local exponent b (previous size) |
|---|---|---|---|---|---|---|---|---|
| 50 M | 50000000 | 11.46 [10.55-11.49] | 0.2773 [0.2624-0.2994] | 0.6392 [0.637-0.7132] | 6.514 [5.923-6.686] | 3.83 [3.716-4.006] | 229 | - |
| 100 M | 100000000 | 24.65 [23.1-26.02] | 0.4869 [0.4602-0.6187] | 1.772 [1.716-1.777] | 14.95 [14.07-15.42] | 7.448 [6.857-8.202] | 247 | 1.105 |
| 200 M | 200000000 | 58.2 [56.49-58.58] | 1.007 [0.9932-1.02] | 5.426 [5.399-5.435] | 32.11 [31.05-32.48] | 19.64 [19.01-19.71] | 291 | 1.239 |
| 500 M | 500000000 | 160.1 [157.6-160.6] | 2.397 [2.358-2.42] | 28.77 [28.63-29.62] | 74.17 [73.75-75.54] | 52.81 [52.56-55.25] | 320 | 1.104 |

Exponents b of t = a n^b (least squares in log-log over 4 sizes):

| series | b (medians) | b (per-size min) | b (per-size max) | a (medians) | ESTIMATE t(3x10^9) from the medians, min fit, max fit |
|---|---|---|---|---|---|
| dna A | **1.153** | 1.183 | 1.146 | 1.52e-08 | 1278 s / 1334 s / 1274 s (ESTIMATE, not measured) |

Per-phase exponents (medians):

| phase | b | share of the total at the largest size |
|---|---|---|
| PH0 | 0.948 | 1.5 % |
| PH1 | 1.656 | 18.0 % |
| PH2/PH3 | 1.056 | 46.3 % |
| PH4+PH5 | 1.162 | 33.0 % |
| total | 1.153 | 100.0 % |
| PH1+PH2/PH3 | 1.157 | |

### dna series B (4 threads, cold cache, runs per size: 3, 3, 3, 3)

| size | n | total s | PH0 s | PH1 s | PH2/PH3 s | PH4+PH5 s | ns / symbol | local exponent b (previous size) |
|---|---|---|---|---|---|---|---|---|
| 50 M | 50000000 | 10.37 [10.37-11.03] | 0.2911 [0.262-0.3673] | 0.9887 [0.9249-1.012] | 5.528 [5.303-5.787] | 3.797 [3.625-3.889] | 207 | - |
| 100 M | 100000000 | 24.39 [23.6-26.72] | 0.5891 [0.4652-1.224] | 1.68 [1.659-1.71] | 14.63 [14.46-14.94] | 7.513 [6.97-8.873] | 244 | 1.234 |
| 200 M | 200000000 | 47.54 [46.79-49.4] | 1.016 [0.9068-1.073] | 2.968 [2.925-3.086] | 25.95 [25.69-26.75] | 17.61 [17.27-18.49] | 238 | 0.962 |
| 500 M | 500000000 | 102.3 [99.52-102.4] | 2.426 [2.291-2.492] | 6.571 [6.499-6.714] | 44.35 [43.48-45.58] | 47.57 [47.25-48.96] | 205 | 0.836 |

Exponents b of t = a n^b (least squares in log-log over 4 sizes):

| series | b (medians) | b (per-size min) | b (per-size max) | a (medians) | ESTIMATE t(3x10^9) from the medians, min fit, max fit |
|---|---|---|---|---|---|
| dna B | **0.985** | 0.977 | 0.954 | 2.93e-07 | 640 s / 612 s / 610 s (ESTIMATE, not measured) |

Per-phase exponents (medians):

| phase | b | share of the total at the largest size |
|---|---|---|
| PH0 | 0.909 | 2.4 % |
| PH1 | 0.824 | 6.4 % |
| PH2/PH3 | 0.886 | 43.3 % |
| PH4+PH5 | 1.110 | 46.5 % |
| total | 0.985 | 100.0 % |
| PH1+PH2/PH3 | 0.877 | |

#### genome-like, the sizes of the task (50 / 100 / 200 M; series A has 1 run per size)

### genome series A (4 threads, cold cache, runs per size: 1, 1, 1)

| size | n | total s | PH0 s | PH1 s | PH2/PH3 s | PH4+PH5 s | ns / symbol | local exponent b (previous size) |
|---|---|---|---|---|---|---|---|---|
| 50 M | 50000000 | 10.55 | 0.2953 | 0.9385 | 5.766 | 3.547 | 211 | - |
| 100 M | 100000000 | 24.74 | 0.4861 | 2.099 | 14.31 | 7.851 | 247 | 1.230 |
| 200 M | 200000000 | 57.4 | 0.9598 | 6.235 | 32.03 | 18.17 | 287 | 1.214 |

Exponents b of t = a n^b (least squares in log-log over 3 sizes):

| series | b (medians) | b (per-size min) | b (per-size max) | a (medians) | ESTIMATE t(3x10^9) from the medians, min fit, max fit |
|---|---|---|---|---|---|
| genome A | **1.222** | 1.222 | 1.222 | 4.12e-09 | 1574 s / 1574 s / 1574 s (ESTIMATE, not measured) |

Per-phase exponents (medians):

| phase | b | share of the total at the largest size |
|---|---|---|
| PH0 | 0.850 | 1.7 % |
| PH1 | 1.366 | 10.9 % |
| PH2/PH3 | 1.237 | 55.8 % |
| PH4+PH5 | 1.178 | 31.7 % |
| total | 1.222 | 100.0 % |
| PH1+PH2/PH3 | 1.257 | |

### genome series B (4 threads, cold cache, runs per size: 3, 3, 3)

| size | n | total s | PH0 s | PH1 s | PH2/PH3 s | PH4+PH5 s | ns / symbol | local exponent b (previous size) |
|---|---|---|---|---|---|---|---|---|
| 50 M | 50000000 | 10.96 [10.95-11.47] | 0.2867 [0.2806-0.3048] | 1.215 [1.208-1.23] | 5.667 [5.562-6.111] | 3.858 [3.754-3.897] | 219 | - |
| 100 M | 100000000 | 24.54 [24.06-25.39] | 0.5265 [0.4792-1.149] | 2.159 [2.155-2.279] | 14.55 [14.2-14.73] | 7.232 [7.219-7.311] | 245 | 1.163 |
| 200 M | 200000000 | 48.5 [47.74-49.1] | 0.9839 [0.9835-1.112] | 4.128 [4.034-4.184] | 25.56 [25.01-26.1] | 17.83 [17.58-17.83] | 243 | 0.983 |

Exponents b of t = a n^b (least squares in log-log over 3 sizes):

| series | b (medians) | b (per-size min) | b (per-size max) | a (medians) | ESTIMATE t(3x10^9) from the medians, min fit, max fit |
|---|---|---|---|---|---|
| genome B | **1.073** | 1.062 | 1.049 | 6.11e-08 | 906 s / 862 s / 860 s (ESTIMATE, not measured) |

Per-phase exponents (medians):

| phase | b | share of the total at the largest size |
|---|---|---|
| PH0 | 0.890 | 2.0 % |
| PH1 | 0.882 | 8.5 % |
| PH2/PH3 | 1.087 | 52.7 % |
| PH4+PH5 | 1.104 | 36.8 % |
| total | 1.073 | 100.0 % |
| PH1+PH2/PH3 | 1.054 | |

#### genome-like series B with the 500 M run added
### genome series B (4 threads, cold cache, runs per size: 3, 3, 3, 3)

| size | n | total s | PH0 s | PH1 s | PH2/PH3 s | PH4+PH5 s | ns / symbol | local exponent b (previous size) |
|---|---|---|---|---|---|---|---|---|
| 50 M | 50000000 | 10.96 [10.95-11.47] | 0.2867 [0.2806-0.3048] | 1.215 [1.208-1.23] | 5.667 [5.562-6.111] | 3.858 [3.754-3.897] | 219 | - |
| 100 M | 100000000 | 24.54 [24.06-25.39] | 0.5265 [0.4792-1.149] | 2.159 [2.155-2.279] | 14.55 [14.2-14.73] | 7.232 [7.219-7.311] | 245 | 1.163 |
| 200 M | 200000000 | 48.5 [47.74-49.1] | 0.9839 [0.9835-1.112] | 4.128 [4.034-4.184] | 25.56 [25.01-26.1] | 17.83 [17.58-17.83] | 243 | 0.983 |
| 500 M | 500000000 | 107.7 [106.3-109.6] | 2.446 [2.379-2.604] | 9.782 [9.576-9.815] | 47.68 [46.45-48.65] | 47.86 [47.82-48.5] | 215 | 0.871 |

Exponents b of t = a n^b (least squares in log-log over 4 sizes):

| series | b (medians) | b (per-size min) | b (per-size max) | a (medians) | ESTIMATE t(3x10^9) from the medians, min fit, max fit |
|---|---|---|---|---|---|
| genome B | **0.987** | 0.983 | 0.974 | 2.91e-07 | 665 s / 649 s / 656 s (ESTIMATE, not measured) |

Per-phase exponents (medians):

| phase | b | share of the total at the largest size |
|---|---|---|
| PH0 | 0.931 | 2.3 % |
| PH1 | 0.910 | 9.1 % |
| PH2/PH3 | 0.906 | 44.3 % |
| PH4+PH5 | 1.113 | 44.4 % |
| total | 0.987 | 100.0 % |
| PH1+PH2/PH3 | 0.906 | |


**Reading.**

* **dna series B (ratio emulation): b = 0.985 [0.954-0.977], met.** The local exponent falls with n (1.23, 0.96, 0.84): the per-symbol time is 207 / 244 / 238 / 205 ns, i.e. flat, and at 500 M lower than at 100 M. Nothing in the
  pipeline grows faster than linear in this series (PH0 0.91, PH1 0.82, PH2/PH3 0.89, PH4+PH5 1.11; PH4+PH5 is the largest at 500 M: 46.5 % of the total).
* **genome-like series B: b = 1.073 [1.049-1.062] over the sizes the task fixed (50 / 100 / 200 M), met;** with the 500 M run added (it was run, three runs, 107.7 s) b = 0.987 [0.974-0.983]. The local exponents
  are 1.16, 0.98 and 0.87.
* **dna series A (fixed parameters): b = 1.153 [1.146-1.183], gate (<= 1.15) not met, by 0.003 on the medians and by 0.033 on the per-size minima.** It lies between 1.15 and 1.3, so by the rule of the task it is "gate not met,
  not an abandon". The cause is in the per-phase exponents: **PH1 b = 1.66** (5.6 % of the time at 50 M, 18.0 % at 500 M: with K_s fixed at 34,953 the number of sample buckets, 460 at 500 M, grows with n and every bucket is one
  scan of the sample, so PH1 is quadratic in n), **PH4+PH5 b = 1.16** (P = n / R = 1,875 chunk files at 500 M, 750 at 200 M: more files, more buffers, more seeks) and PH2/PH3 1.06. With three sizes (50-200 M) b is 1.172
  [1.175-1.210] (the fit of the first E7 batch; the fourth size lowers it). Series A is not the regime of the real 3x10^9 run (there K_s, R and K are the reference values, which is what series B emulates); it shows
  what fixed buffers cost, and that PH1 is the first stage that must scale its bucket size with n.
* **genome-like series A** (one run per size, not gated): b = 1.222 over 50 / 100 / 200 M (PH1 1.37, PH2/PH3 1.24, PH4+PH5 1.18), the same picture; **not run at 500 M.**
* Threads (one run each, series B): 1 / 2 / 4 threads at 500 M dna 184.8 / 119.6 / 102.3 s (4 threads 1.81x faster than 1), at 100 M 39.7 / 27.2 / 24.4 s (E6: 39.3 / 28.6 / 24.7); table in section 6.

**ESTIMATES for 3x10^9 (extrapolated 6x beyond the largest measured size; not measured).** Series B, a n^b with the medians / per-size minima / maxima: dna **640 s [610-640] = 10.7 min**, genome-like (50-200 M fit) 906 s
[860-906], with the 500 M run 665 s [649-665]; the linear extrapolation of E6 was 12.4 min (12-15). Two reasons not to take the lower numbers at face value: E4 measured the per-symbol time of PH4 at real block sizes up to 1.5x the
scale-down one, which on the 46 % share of PH4+PH5 gives up to +23 % (about 13 min), and the real run has the real block sizes. **A reasonable reading is 11-15 minutes; the real
run (E8) decides.** For series A the fit gives 1278 s [1274-1334]; that number has no meaning for a real build (it extrapolates fixed buffers).

## 5. S2: RssAnon (gate: <= 256 MB at every n in both series; series A slope <= 5 MB per doubling, abandon above 15; series B slope reported)

Peak RssAnon of the whole pipeline (program sampler at 10 ms; the independent sampler of `/proc/<pid>/status` at 20 ms read the same value in every cell), median [min-max], and the per-phase peaks:


### dna series A: RssAnon (MB, median [min-max])

| size | program sampler (10 ms) | independent sampler (20 ms) | PH0 | PH1 | PH2/PH3 | PH4+PH5 | arena peak MB |
|---|---|---|---|---|---|---|---|
| 50 M | 74.26 [72.22-75.55] | 74.26 [72.22-75.56] | 14.1 [13.7-14.2] | 23.4 [23-23.5] | 55.6 [54.1-57.1] | 74.3 [72.2-75.6] | 18.7 [18.7-18.7] |
| 100 M | 96.69 [86.09-100.4] | 96.69 [86.09-100.4] | 15 [13.7-15.1] | 30.2 [29-30.5] | 82.3 [57.1-89.1] | 96.7 [86.1-100] | 30.9 [30.9-30.9] |
| 200 M | 132.4 [125.1-134.3] | 132.4 [125.1-134.3] | 13.9 [13.9-15.1] | 34 [34-35.3] | 85.9 [82.8-89.2] | 132 [125-134] | 55.5 [55.5-55.5] |
| 500 M | 192.6 [187.5-214.5] | 192.6 [187.5-214.5] | 13.9 [13.9-14.4] | 35.2 [35.2-35.8] | 85.5 [79.4-110] | 193 [187-214] | 129 [129-129] |

Slope of the median peak (program sampler) against log2(n): **35.99 MB per doubling** (least squares over 4 sizes; the last doubling alone: 60.26 MB)
Slope with the independent sampler: 35.99 MB per doubling
Doublings step by step: 50 M -> 100 M: +22.4 MB, 100 M -> 200 M: +35.7 MB, 200 M -> 500 M: +45.6 MB

### dna series B: RssAnon (MB, median [min-max])

| size | program sampler (10 ms) | independent sampler (20 ms) | PH0 | PH1 | PH2/PH3 | PH4+PH5 | arena peak MB |
|---|---|---|---|---|---|---|---|
| 50 M | 93.94 [92.3-95.35] | 93.93 [92.3-95.35] | 13.5 [13.5-14.5] | 29 [29-29.7] | 65 [64.5-66.1] | 93.9 [92.3-95.3] | 28.8 [28.8-28.8] |
| 100 M | 93.45 [92.84-105.7] | 93.45 [92.84-105.7] | 13.8 [13.7-15.5] | 29.1 [29-30.7] | 67.6 [64.9-94.7] | 93.4 [92.8-106] | 30.9 [30.9-30.9] |
| 200 M | 104.5 [97.36-110.1] | 104.5 [97.36-110.1] | 13.9 [13.9-15.5] | 29.7 [29.7-31.3] | 75.6 [74-106] | 104 [97.4-110] | 35.2 [35.2-35.2] |
| 500 M | 124.1 [108-124.4] | 124.1 [108-124.4] | 13.9 [13.8-13.9] | 31.5 [31.5-31.5] | 89.2 [86.7-124] | 114 [108-124] | 48 [48-48] |

Slope of the median peak (program sampler) against log2(n): **9.48 MB per doubling** (least squares over 4 sizes; the last doubling alone: 19.67 MB)
Slope with the independent sampler: 9.48 MB per doubling
Doublings step by step: 50 M -> 100 M: -0.5 MB, 100 M -> 200 M: +11.0 MB, 200 M -> 500 M: +14.9 MB

### genome series A: RssAnon (MB, median [min-max])

| size | program sampler (10 ms) | independent sampler (20 ms) | PH0 | PH1 | PH2/PH3 | PH4+PH5 | arena peak MB |
|---|---|---|---|---|---|---|---|
| 50 M | 80.17 | 80.17 | 13.9 | 29.5 | 63.9 | 80.2 | 18.7 |
| 100 M | 92.39 | 92.39 | 13.9 | 34.1 | 84.5 | 92.4 | 30.9 |
| 200 M | 146.1 | 146.1 | 14.1 | 43.3 | 123 | 146 | 55.5 |

Slope of the median peak (program sampler) against log2(n): **32.98 MB per doubling** (least squares over 3 sizes; the last doubling alone: 53.75 MB)
Slope with the independent sampler: 32.98 MB per doubling
Doublings step by step: 50 M -> 100 M: +12.2 MB, 100 M -> 200 M: +53.8 MB

### genome series B: RssAnon (MB, median [min-max])

| size | program sampler (10 ms) | independent sampler (20 ms) | PH0 | PH1 | PH2/PH3 | PH4+PH5 | arena peak MB |
|---|---|---|---|---|---|---|---|
| 50 M | 98.13 [94.85-98.68] | 98.13 [94.85-98.68] | 13.9 [13.9-14.4] | 31.2 [31.2-31.6] | 68.1 [67.1-71.2] | 98.1 [94.8-98.7] | 28.8 [28.8-28.8] |
| 100 M | 101.4 [100.1-104] | 101.4 [100.1-103.9] | 13.9 [13.8-14.8] | 34.1 [34-37.6] | 79.7 [78.8-87.3] | 101 [100-104] | 30.9 [30.9-30.9] |
| 200 M | 109.1 [105.8-132.5] | 109.1 [105.8-132.5] | 14.2 [14-14.5] | 39.4 [38.9-43.6] | 97.1 [79.4-106] | 109 [106-133] | 35.2 [35.2-35.2] |
| 500 M | 127.2 [109.7-129] | 127.2 [109.7-129] | 14 [14-14] | 46.8 [46.1-48.9] | 101 [99.5-117] | 127 [110-129] | 48 [48-48] |

Slope of the median peak (program sampler) against log2(n): **8.84 MB per doubling** (least squares over 4 sizes; the last doubling alone: 18.17 MB)
Slope with the independent sampler: 8.84 MB per doubling
Doublings step by step: 50 M -> 100 M: +3.2 MB, 100 M -> 200 M: +7.7 MB, 200 M -> 500 M: +13.7 MB


* **Peak over all 4-thread runs: series A 214.5 MB** (dna 500 M run 3; median 192.6), **series B 124.4 MB** (dna 500 M) and 132.5 MB (genome-like 200 M, one run; medians 127.2 MB at 500 M): **<= 256 MB at every n in both
  series: met.** The peak is the fused PH4 + PH5 phase in the medians; in a few runs the PH2/PH3 phase has a spike above its usual 65-90 MB (one run each at 200 M and 500 M, up to 124 MB; not investigated).
* **Series A slope: 36.0 MB per doubling (dna, four sizes; 22.4, 35.7 and 45.6 MB per doubling between consecutive sizes) and 33.0 (genome-like, three sizes): gate (<= 5) not met, and the abandon rule "series A
  slope > 15 MB per doubling" fires.** Cause, read from the code and the arena rows: PH4's writers (`SpoolWriter`) hold one buffer of `--tuple-kb` (32 KB) **per chunk**; in B2 each of the 2 workers has one, so the native buffers are about 2 x 32 KB x P
  with P = n / R chunks. At fixed R = 266,667 the arena peak is 18.7 / 30.9 / 55.5 / **129 MB** at P = 188 / 375 / 750 / 1,875 (69 KB per chunk); its steps (+12, +25, +74 MB) are of the size of the steps of the peak RssAnon (+22, +36, +61 MB), i.e. it explains most of the growth; the cap of 160 MiB would stop the build (`ArenaBudget` throws) above P of about 2,400, i.e. a text of about 640 M symbols at this R. In series B P stays 375 (R grows with n), the arena peak is
  28.8 / 30.9 / 35.2 / 48 MB (the arrays of R entries per B2 worker grow with R).
* **Series B slope: 9.5 MB per doubling (dna) and 8.8 (genome-like), reported.** Plainly: **part of this growth follows from parameters that scale with n** (R = 8,000,000 f grows with n, the Phi and rank arrays of the 2 B2
  workers with it; PH3's K = n / 3000 too), not from anything that grows at the real 3x10^9, where R is 8,000,000 and the arena is capped; the steps are -0.5 (50 -> 100 M), +11.0 (100 -> 200 M) and +14.9 MB per doubling
  (200 -> 500 M), so the slope is not constant and is largest where f approaches its real value. **ESTIMATE** for 3x10^9 by extending the last step: 124 + 2.6 doublings x 15 = about 165 MB (not measured; E4's PH4 alone
  measured 154 MB at a real 500 Mbp, so the real 3x10^9 value, E8, must be measured). The abandon rule does not fire for series B (it names series A).
* What the 40 MB "residue" of E6 does: the PH4 + PH5 peak minus the PH2/PH3 peak is 25-29 MB at the four sizes of series B; not investigated (the slope test does not point at it: the growth is the arena).

## 6. S3: identity, threads, series

* **54 timed runs of the scaling groups** (dna A and B, 50 / 100 / 200 / 500 M, 3 runs each = 24; dna B with 1 and 2 threads at the four sizes = 8; genome-like B at four sizes = 12, genome-like A 3, genome-like B with 1 and 2
  threads at 50 / 100 / 200 M = 6; the r500 cross-check 1), **430 identity comparisons (pk, SA, LCP base, LCP overflow, CHILD, the four statistics, exit code, scratch gone): all equal, none different.** Because every run is compared with
  the hash of one reference per input, series A, series B and 1 / 2 / 4 threads produce identical files at every n. The 43 runs of the cache sweep and its calibration that finished (26 + 17, section 8) were compared the same way: all equal; 7 runs were killed (watchdog or OOM) and have no outputs.
* **At 500 M** (dna): the references were the in-RAM kernels (SaisKernel SA 2.0 GB, LcpKernel LCP, the product's `ChildTable` and statistics: NODE_COUNT 810,933,159, internal 310,933,157, LRS_DEPTH 29, LRS_RANK 461,412,445);
  the disk allowed all of them, so **CHILD was compared with the e5-oracle reference** (not with a second run): 3 + 3 runs at 4 threads (series B and A) and 1 + 1 runs at 1 and 2 threads, each against the same hashes.
  The genome-like 500 M text (NODE_COUNT 829,508,909, LRS_DEPTH 2,499,999) the same: 3 runs equal to all five references.
* **Cross-check against the E4 files (`orch4/r500.*`), run once:** the packed text of the pipeline on `r500.txt` equals `r500.pk`, its suffix array equals `r500.sa` (the E3 builder's SA of the E4 verification) and its LCP
  and overflow files equal `r500.oracle.lcp` / `.ovf`: 4 of 4 byte-identical (there is no CHILD reference in `orch4`, shown `n/a`). How they relate: `r500.pk` is a `packrand` text (seed 20240009), its symbols are unpacked MSB-first
  as ACGT (`e7_unpack_packed.py`), so the pipeline's own packing reproduces the file, which was checked by that equality, not assumed. It is the same in-repo code family (E3 builder, LcpKernel), so independent of PH1-PH5 of
  the pipeline but not an independent implementation.

**Threads and series (one table, `tables_e7.py threads`):** total time, per-phase times and RssAnon of every group (median [min-max] where 3 runs):

| input | series | threads | n | total s | PH2/PH3 s | PH4+PH5 s | RssAnon MB (program) | runs |
|---|---|---|---|---|---|---|---|---|
| dna_50m | A | 4 | 50000000 | 11.46 [10.55-11.49] | 6.514 [5.923-6.686] | 3.83 [3.716-4.006] | 74.26 [72.22-75.55] | 3 |
| dna_100m | A | 4 | 100000000 | 24.65 [23.1-26.02] | 14.95 [14.07-15.42] | 7.448 [6.857-8.202] | 96.69 [86.09-100.4] | 3 |
| dna_200m | A | 4 | 200000000 | 58.2 [56.49-58.58] | 32.11 [31.05-32.48] | 19.64 [19.01-19.71] | 132.4 [125.1-134.3] | 3 |
| dna_500m | A | 4 | 500000000 | 160.1 [157.6-160.6] | 74.17 [73.75-75.54] | 52.81 [52.56-55.25] | 192.6 [187.5-214.5] | 3 |
| dna_50m | B | 4 | 50000000 | 10.37 [10.37-11.03] | 5.528 [5.303-5.787] | 3.797 [3.625-3.889] | 93.94 [92.3-95.35] | 3 |
| dna_50m | B | 2 | 50000000 | 11.5 | 6.342 | 3.424 | 67.06 | 1 |
| dna_50m | B | 1 | 50000000 | 16.96 | 9.005 | 5.081 | 45.47 | 1 |
| dna_100m | B | 4 | 100000000 | 24.39 [23.6-26.72] | 14.63 [14.46-14.94] | 7.513 [6.97-8.873] | 93.45 [92.84-105.7] | 3 |
| dna_100m | B | 2 | 100000000 | 27.24 | 17.12 | 6.947 | 68.93 | 1 |
| dna_100m | B | 1 | 100000000 | 39.68 | 24.5 | 9.946 | 54.08 | 1 |
| dna_200m | B | 4 | 200000000 | 47.54 [46.79-49.4] | 25.95 [25.69-26.75] | 17.61 [17.27-18.49] | 104.5 [97.36-110.1] | 3 |
| dna_200m | B | 2 | 200000000 | 53.55 | 29.49 | 18.07 | 72.7 | 1 |
| dna_200m | B | 1 | 200000000 | 77.89 | 44.15 | 23.59 | 54.94 | 1 |
| dna_500m | B | 4 | 500000000 | 102.3 [99.52-102.4] | 44.35 [43.48-45.58] | 47.57 [47.25-48.96] | 124.1 [108-124.4] | 3 |
| dna_500m | B | 2 | 500000000 | 119.6 | 58.39 | 46.84 | 84.64 | 1 |
| dna_500m | B | 1 | 500000000 | 184.8 | 94.12 | 65.25 | 63.68 | 1 |
| genome_50m | A | 4 | 50000000 | 10.55 | 5.766 | 3.547 | 80.17 | 1 |
| genome_100m | A | 4 | 100000000 | 24.74 | 14.31 | 7.851 | 92.39 | 1 |
| genome_200m | A | 4 | 200000000 | 57.4 | 32.03 | 18.17 | 146.1 | 1 |
| genome_50m | B | 4 | 50000000 | 10.96 [10.95-11.47] | 5.667 [5.562-6.111] | 3.858 [3.754-3.897] | 98.13 [94.85-98.68] | 3 |
| genome_50m | B | 2 | 50000000 | 13.15 | 7.271 | 3.864 | 72.27 | 1 |
| genome_50m | B | 1 | 50000000 | 19.9 | 11.63 | 4.958 | 46.75 | 1 |
| genome_100m | B | 4 | 100000000 | 24.54 [24.06-25.39] | 14.55 [14.2-14.73] | 7.232 [7.219-7.311] | 101.4 [100.1-104] | 3 |
| genome_100m | B | 2 | 100000000 | 28.84 | 18.01 | 6.995 | 78.4 | 1 |
| genome_100m | B | 1 | 100000000 | 43.35 | 27.01 | 10.29 | 53.21 | 1 |
| genome_200m | B | 4 | 200000000 | 48.5 [47.74-49.1] | 25.56 [25.01-26.1] | 17.83 [17.58-17.83] | 109.1 [105.8-132.5] | 3 |
| genome_200m | B | 2 | 200000000 | 56.44 | 30.92 | 17.89 | 86.29 | 1 |
| genome_200m | B | 1 | 200000000 | 88.54 | 51.89 | 24.51 | 56.51 | 1 |
| genome_500m | B | 4 | 500000000 | 107.7 [106.3-109.6] | 47.68 [46.45-48.65] | 47.86 [47.82-48.5] | 127.2 [109.7-129] | 3 |


## 7. S4: I/O, scratch, footprint (gate: I/O <= 100 B/char)

| case | n | I/O B/char (program counters) | rchar + wchar B/char (program / independent whole process) | mapped write B/char | scratch peak B/char (tracked / apparent / allocated) | output footprint B/char | block-layer read B/char | dirty peak MB (system / mapped outputs) | disk free at start GB |
|---|---|---|---|---|---|---|---|---|---|
| dna_50m.e7.B | 50000000 | 65.87 [65.87-65.87] | 59.58 [59.55-59.76] / 59.58 [59.55-59.76] | 9.38 [9.38-9.38] | 12 [12-12] / 11.65 [11.58-11.69] / 11.66 [11.6-11.69] | 9.25 [9.25-9.25] | 1.04 [1.03-1.05] | 602 [601-603] / 9.31 [9.18-9.44] | 12 [12-12] |
| dna_50m.e7.A | 50000000 | 68.64 [68.64-68.64] | 62.64 [62.38-62.67] / 62.64 [62.38-62.67] | 9.38 [9.38-9.38] | 12 [12-12] / 11.74 [11.66-12] / 11.76 [11.67-12.01] | 9.25 [9.25-9.25] | 1.04 [1.03-1.04] | 601 [599-602] / 9.31 [9.18-9.44] | 12 [12-12] |
| dna_100m.e7.B | 100000000 | 71.59 [71.59-71.59] | 65.77 [65.65-66.09] / 65.78 [65.65-66.09] | 9.38 [9.38-9.38] | 12 [12-12] / 11.81 [11.71-11.92] / 11.82 [11.72-11.92] | 9.25 [9.25-9.25] | 1.02 [1.02-1.02] | 1.2e+03 [1.2e+03-1.2e+03] / 9.44 [9.31-9.44] | 12 [12-12] |
| dna_100m.e7.A | 100000000 | 71.59 [71.59-71.59] | 65.8 [65.59-66.01] / 65.8 [65.59-66.01] | 9.38 [9.38-9.38] | 12 [12-12] / 11.81 [11.8-11.84] / 11.82 [11.81-11.85] | 9.25 [9.25-9.25] | 1.02 [1.02-1.02] | 1.2e+03 [1.2e+03-1.2e+03] / 9.31 [9.25-9.44] | 12 [12-12] |
| dna_200m.e7.B | 200000000 | 73.77 [73.77-73.77] | 67.9 [67.78-68.02] / 67.9 [67.78-68.02] | 9.38 [9.38-9.38] | 12 [12-12] / 11.98 [11.87-12] / 11.98 [11.88-12] | 9.25 [9.25-9.25] | 1.01 [1.01-1.01] | 1.57e+03 [1.57e+03-1.6e+03] / 9.44 [9.44-9.44] | 12 [12-12] |
| dna_200m.e7.A | 200000000 | 74.25 [74.25-74.25] | 69.07 [68.98-69.09] / 69.07 [68.98-69.09] | 9.38 [9.38-9.38] | 12 [12-12] / 11.92 [11.86-11.95] / 11.93 [11.87-11.96] | 9.25 [9.25-9.25] | 1.01 [1.01-1.01] | 1.61e+03 [1.57e+03-1.61e+03] / 9.44 [9.44-9.44] | 12 [12-12] |
| dna_50m.e7.B.t1 | 50000000 | 65.87 | 61.19 / 61.19 | 9.38 | 12 / 11.98 / 11.99 | 9.25 | 1.04 | 601 / 9.18 | 12 |
| dna_50m.e7.B.t2 | 50000000 | 65.87 | 59.78 / 59.78 | 9.38 | 12 / 11.52 / 11.53 | 9.25 | 1.05 | 601 / 9.44 | 12 |
| dna_100m.e7.B.t1 | 100000000 | 71.59 | 67.65 / 67.65 | 9.38 | 12 / 11.81 / 11.82 | 9.25 | 1.02 | 1.2e+03 / 9.31 | 12 |
| dna_100m.e7.B.t2 | 100000000 | 71.59 | 66.03 / 66.03 | 9.38 | 12 / 11.81 / 11.82 | 9.25 | 1.02 | 1.2e+03 / 9.44 | 12 |
| dna_200m.e7.B.t1 | 200000000 | 73.77 | 69.71 / 69.71 | 9.38 | 12 / 11.99 / 12 | 9.25 | 1.06 | 1.61e+03 / 9.44 | 12 |
| dna_200m.e7.B.t2 | 200000000 | 73.77 | 68.14 / 68.14 | 9.38 | 12 / 12 / 12 | 9.25 | 1.06 | 1.62e+03 / 9.44 | 12 |
| dna_500m.e7.B | 500000000 | 74.08 [74.08-74.08] | 67.64 [67.55-67.67] / 67.64 [67.55-67.67] | 9.38 [9.38-9.38] | 12 [12-12] / 11.97 [11.94-12] / 11.97 [11.94-12] | 9.25 [9.25-9.25] | 1.02 [1.02-1.02] | 1.61e+03 [1.61e+03-1.62e+03] / 9.44 [9.44-9.44] | 16 [16-16] |
| dna_500m.e7.A | 500000000 | 75.77 [75.77-75.77] | 70.97 [70.94-71.05] / 70.97 [70.94-71.05] | 9.38 [9.38-9.38] | 12 [12-12] / 11.98 [11.97-12] / 11.99 [11.99-12.01] | 9.25 [9.25-9.25] | 1.02 [1.02-1.02] | 1.61e+03 [1.61e+03-1.62e+03] / 9.44 [9.44-9.44] | 16 [16-16] |
| dna_500m.e7.B.t1 | 500000000 | 74.08 | 69.74 / 69.74 | 9.38 | 12 / 11.97 / 11.97 | 9.25 | 1.02 | 1.62e+03 / 9.44 | 16 |
| dna_500m.e7.B.t2 | 500000000 | 74.08 | 68.01 / 68.01 | 9.38 | 12 / 12 / 12 | 9.25 | 1.02 | 1.61e+03 / 9.44 | 16 |
| genome_50m.e7.B | 50000000 | 67.43 [67.43-67.43] | 60.94 [60.92-61.08] / 60.94 [60.92-61.08] | 9.76 [9.76-9.76] | 12 [12-12] / 11.55 [11.42-11.87] / 11.55 [11.43-11.88] | 9.502 [9.502-9.502] | 1.24 [1.23-1.24] | 603 [603-603] / 10.1 [10.1-10.2] | 15 [15-15] |
| genome_100m.e7.B | 100000000 | 73.06 [73.06-73.06] | 66.89 [66.83-67] / 66.89 [66.83-67] | 9.75 [9.75-9.75] | 12 [12-12] / 11.81 [11.79-11.84] / 11.82 [11.8-11.85] | 9.486 [9.486-9.486] | 1.12 [1.12-1.12] | 1.2e+03 [1.2e+03-1.2e+03] / 10.1 [10-10.1] | 15 [15-15] |
| genome_200m.e7.B | 200000000 | 75.45 [75.45-75.45] | 69.25 [69.21-69.28] / 69.25 [69.21-69.28] | 9.75 [9.75-9.75] | 12 [12-12] / 11.94 [11.89-11.97] / 11.94 [11.89-11.97] | 9.491 [9.491-9.491] | 1.06 [1.06-1.06] | 1.61e+03 [1.6e+03-1.61e+03] / 10.5 [10.5-10.5] | 15 [15-15] |
| genome_50m.e7.A | 50000000 | 70.04 | 63.41 / 63.41 | 9.76 | 12 / 11.8 / 11.81 | 9.502 | 1.24 | 602 / 9.83 | 15 |
| genome_100m.e7.A | 100000000 | 73.06 | 66.93 / 66.94 | 9.75 | 12 / 11.8 / 11.81 | 9.486 | 1.12 | 1.2e+03 / 10.5 | 15 |
| genome_200m.e7.A | 200000000 | 75.93 | 70.36 / 70.36 | 9.75 | 12 / 11.97 / 11.98 | 9.491 | 1.06 | 1.6e+03 / 10.5 | 15 |
| genome_50m.e7.B.t1 | 50000000 | 67.43 | 63.21 / 63.21 | 9.76 | 12 / 11.93 / 11.94 | 9.502 | 1.23 | 601 / 10.2 | 15 |
| genome_50m.e7.B.t2 | 50000000 | 67.43 | 61.42 / 61.42 | 9.76 | 12 / 11.68 / 11.69 | 9.502 | 1.24 | 602 / 9.83 | 15 |
| genome_100m.e7.B.t1 | 100000000 | 73.06 | 69.28 / 69.28 | 9.75 | 12 / 11.87 / 11.88 | 9.486 | 1.12 | 1.2e+03 / 10.1 | 15 |
| genome_100m.e7.B.t2 | 100000000 | 73.06 | 67.4 / 67.4 | 9.75 | 12 / 11.67 / 11.67 | 9.486 | 1.12 | 1.2e+03 / 10.5 | 15 |
| genome_200m.e7.B.t1 | 200000000 | 75.45 | 71.83 / 71.83 | 9.75 | 12 / 11.94 / 11.94 | 9.491 | 1.06 | 1.61e+03 / 10.4 | 15 |
| genome_200m.e7.B.t2 | 200000000 | 75.45 | 69.7 / 69.7 | 9.75 | 12 / 11.94 / 11.94 | 9.491 | 1.06 | 1.61e+03 / 10.4 | 15 |
| genome_500m.e7.B | 500000000 | 75.57 [75.57-75.57] | 68.94 [68.9-68.99] / 68.94 [68.9-68.99] | 9.76 [9.76-9.76] | 12 [12-12] / 11.98 [11.94-11.98] / 11.98 [11.94-11.99] | 9.499 [9.499-9.499] | 1.02 [1.02-1.02] | 1.61e+03 [1.61e+03-1.61e+03] / 10.4 [10.4-10.5] | 15 [15-15] |


**I/O (program counters, explicit I/O + writes through the mapped windows) is 65.9-75.9 B/char in every cell (gate 100: met);** `rchar + wchar` (syscalls only) 59.6-71.0, the same figure from the independent sampler of the whole
process (the two agree to the digit shown; the program's number is 6-7 B/char higher: it adds the writes through the mapped windows, 9.4-9.8 B/char, and its explicit counters are about 5 % below rchar + wchar, as in E6). It grows slowly with n in series B (65.9 -> 74.1 over a factor of 10 in n; not attributed) and a little faster in series A (68.6 -> 75.8). **Scratch peak 12.0 B/char** (tracked; the independent sampler of the directory
11.5-12.0 apparent, the same allocated) at every n: not above E6's figure, and linear (6.0 GB at 500 M). **Output footprint 9.25 B/char (dna), 9.49-9.50 (genome-like)**: 4.6 GB at 500 M. The block layer read
1.0-1.24 B/char in every run: **at 500 M the page cache of the 15 GB machine was not exceeded** (text 0.5 GB + outputs 4.6 GB + scratch up to 6 GB do not exceed it at one time, and the reads of the block layer are only
the cold text, 1.02 B/char as at 50 M); the 500 M runs are therefore not a cache-pressure test, and the task's expectation that "the page cache is really exceeded" did not materialise on this machine. System-wide dirty
pages peak at 0.6-1.6 GB (PH4's scratch; the kernel limit of this machine), the mapped dirty pages of the outputs at 9.2-10.5 MB (E6's clause: <= 64 MB, met).
The disk is a virtio disk, probably host-cached: the scaling runs at 500 M counted 37 GB of I/O in 100-107 s (about 350 MB/s) without any visible limit by the device, which is what a host cache does; real disks are slower
(ESTIMATE: the 8-10 s of 25 s of scratch I/O at 100 M in E6 would become more).

## 8. C: the page-cache test of the whole pipeline (the "true 0.5x test" of section 6b)

**What is limited and how the effective cache is measured.** A cgroup-v1 sub-group (`/sys/fs/cgroup/memory/emproto_e7_*`, created and removed by `run_e7_cache.sh`) holds the whole process; its `memory.limit_in_bytes` is
the TOTAL of anonymous memory, page cache, mapped pages and charged kernel memory, so the limit is not the cache. The **effective cache** is therefore measured, not set: `sample_e7_cache.py` takes the page residency of the packed
text and of the rank file with `mincore` (fresh read-only mapping, nothing touched) every 200 ms while the rank file exists (PH1 end to PH3 end), and the effective cache is the **time-weighted resident fraction of
(text + rank file) over PH3** (it cannot exceed 1: cells at "1x" and above are told apart by their limit). The sampler also records the cgroup's `memory.stat`, `memory.failcnt`, RssAnon (20 ms) and `/proc/<pid>/io`.

**Settings that were needed to get below ~0.85 at all ("lean").** With the E6 settings the process holds 105 MB of anonymous memory and the runtime image's 30-60 MB of file pages, so a limit of 100 MB leaves almost nothing for the
cache (E3 found the same: the nominal 0.5x was 1.3-1.6x effective, the process died below ~98 MB). The test runs therefore use **`DOTNET_GCHeapHardLimit` = 48 MiB, `--spill-kb 16 --tuple-kb 8 --budget-mb 64
--ph5-budget-mb 24 --b2-workers 2`** (4 threads, `--scale-n 3000000000`, series B rule; the output does not depend on any of them), **which cut the peak RssAnon from 105 MB to 43 MB (dna_100m; 54-68 MB genome-like) at a cost of
+12 % time unlimited (27.3 s against 24.4 s)**, and the .NET runtime image (CoreLib, coreclr, the JIT, the app's assemblies) is **read into the page cache from outside the sub-group before the run** so that its pages are charged
to the parent and the limit applies to what the design budgets (anonymous memory and the data files). Every slowdown below is against the **unlimited cell of the same configuration** (a sub-group without a limit, same sampler,
cold cache, 3 runs); the E6 configuration unlimited is 24.4 s (dna_100m) and 24.5 s (genome-like) without a group. Two variants:

* **C (constant limit):** the limit holds for the whole run, PH4 + PH5 included. This is the literal "whole pipeline in a sub-group".
* **D (limited PH0-PH3):** the same limit until the rank file is deleted (PH3 done), then the limit is raised to 1 GiB for PH4 + PH5. It exists because of what C showed (below): PH4's 12 B/char of scratch is a
  separate cache consumer from text + rank, and C cannot reach a true 0.5x without it collapsing.

The cell name is the variant and the limit in MB (D53 = limit 53 MB until PH3 ends). Time series of four runs (every 200 ms: phase flags, residency of text and rank file, the cgroup's memory counters, RssAnon): `e7_series_dna_U.tsv`, `e7_series_dna_D53.tsv`, `e7_series_dna_C65.tsv` (the run that did not finish), `e7_series_genome_D66.tsv`. 1 run per cell, 3 runs at the unlimited, 1x and 0.5x cells (and for the genome-like text at the 0.5x cell); watchdog 15 min.

**Sweep (`e7_cache.tsv`)**, slowdown against the unlimited cell of the same input, "PH2/PH3 slowdown" against its PH2/PH3 time; block-layer read = `/proc/<pid>/io` read_bytes of the whole process; killed = watchdog (15 min) or OOM:

| input | cell | limit MB (raised to MB after PH3) | runs finished / started | total s | slowdown (median) | PH1 s | PH2/PH3 s | PH2/PH3 slowdown | PH4+PH5 s | effective cache (PH3 resident fraction of text + rank) | block-layer read MB (whole process) | memory.failcnt | killed | RssAnon peak MB | identity |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dna_100m | U | unlimited | 3 / 3 | 27.33 [26.99-27.6] | 1 [0.988-1.01] | 1.69 [1.65-1.73] | 17.6 [17-17.6] | 1 [0.965-1] | 7.76 [7.58-7.78] | 1 [1-1] | 105.5 [105.3-105.5] | 0 [0-0] | no | 43.2 [42.3-52.3] | 3 of 3 runs equal to the references |
| dna_100m | C400 | 400 | 1 / 1 | 37.93 | 1.39 | 1.66 | 17.1 | 0.971 | 18.7 | 0.982 | 2653 | 3.487e+05 | no | 48.5 | 1 of 1 runs equal to the references |
| dna_100m | C250 | 250 | 3 / 3 | 37.88 [35.36-41.71] | 1.39 [1.29-1.53] | 1.7 [1.67-1.76] | 18.1 [17.4-18.5] | 1.03 [0.986-1.05] | 17.1 [15.1-22.2] | 0.967 [0.965-0.968] | 2845 [2818-2859] | 3.506e+05 [3.334e+05-3.556e+05] | no | 47.3 [45-47.7] | 3 of 3 runs equal to the references |
| dna_100m | C80 | 80 | 1 / 1 | 70.52 | 2.58 | 1.76 | 23 | 1.31 | 45.2 | 0.777 | 4381 | 5.344e+05 | no | 43.2 | 1 of 1 runs equal to the references |
| dna_100m | D70 | 70 (1024) | 1 / 1 | 33.72 | 1.23 | 1.67 | 23.4 | 1.33 | 8.21 | 0.725 | 3204 | 3.778e+05 | no | 52.1 | 1 of 1 runs equal to the references |
| dna_100m | D53 | 53 (1024) | 3 / 3 | 98.3 [98.03-99.16] | 3.6 [3.59-3.63] | 1.72 [1.69-1.79] | 87.5 [86-88.6] | 4.97 [4.89-5.03] | 8.58 [8.35-9.76] | 0.493 [0.463-0.51] | 4.463e+04 [4.087e+04-4.741e+04] | 4.711e+06 [4.371e+06-4.992e+06] | no | 49 [48.4-51.4] | 3 of 3 runs equal to the references |
| dna_100m | C65 | 65 | 0 / 1 | > watchdog | - | - | - | - | - | 0.652 | 1.136e+06 | 1.789e+08 | watchdog | 48.4 | no finished run |
| genome_100m | U | unlimited | 3 / 3 | 41.9 [38.85-48.46] | 1 [0.927-1.16] | 3.67 [3.52-3.72] | 20.4 [20.3-20.5] | 1 [0.994-1] | 17.4 [14.1-23.9] | 1 [1-1] | 105.2 [103.4-105.5] | 0 [0-0] | no | 53.8 [53.6-57.9] | 3 of 3 runs equal to the references |
| genome_100m | C400 | 400 | 1 / 1 | 51.44 | 1.23 | 3.83 | 20.2 | 0.99 | 26.8 | 0.99 | 2739 | 3.116e+05 | no | 64.4 | 1 of 1 runs equal to the references |
| genome_100m | C250 | 250 | 3 / 3 | 40.5 [40.28-43.33] | 0.967 [0.961-1.03] | 2.39 [2.12-3.8] | 18.4 [17.8-19.9] | 0.899 [0.873-0.973] | 19.1 [19-20] | 0.974 [0.97-0.974] | 2847 [2832-2859] | 3.727e+05 [3.577e+05-3.777e+05] | no | 67.4 [57.8-67.8] | 3 of 3 runs equal to the references |
| genome_100m | D90 | 90 (1024) | 1 / 1 | 157.7 [95.83-219.7] | 3.76 [2.29-5.24] | 2.21 [2.14-2.28] | 146 [83.7-208] | 7.12 [4.09-10.2] | 9.32 [9.17-9.47] | 0.753 [0.667-0.839] | 1.341e+05 [4.206e+04-2.26e+05] | 1.919e+07 [6.428e+06-3.194e+07] | no | 63.2 [53.6-72.7] | 1 of 1 runs equal to the references |
| genome_100m | D83 | 83 (1024) | 1 / 1 | 222.9 | 5.32 | 2.23 | 211 | 10.3 | 9.25 | 0.652 | 2.431e+05 | 3.38e+07 | no | 68.1 | 1 of 1 runs equal to the references |
| genome_100m | D66 | 66 (1024) | 3 / 3 | 356 [349-358.6] | 8.5 [8.33-8.56] | 2.53 [2.36-2.71] | 344 [337-347] | 16.8 [16.5-17] | 8.96 [8.53-9.04] | 0.512 [0.508-0.544] | 4.508e+05 [2.974e+05-4.532e+05] | 5.874e+07 [3.98e+07-5.956e+07] | no | 59.9 [37.5-64.2] | 3 of 3 runs equal to the references |
| genome_100m | D100 | 100 (1024) | 1 / 1 | 117.4 | 2.8 | 2.28 | 106 | 5.17 | 8.95 | 0.814 | 6.507e+04 | 9.983e+06 | no | 72.5 | 1 of 1 runs equal to the references |
| genome_100m | D120 | 120 (1024) | 1 / 1 | 34.28 | 0.818 | 2.35 | 22.5 | 1.1 | 8.84 | 0.964 | 2472 | 2.518e+05 | no | 59.5 | 1 of 1 runs equal to the references |
| genome_100m | D75 | 75 (1024) | 0 / 1 | > watchdog | - | - | - | - | - | 0.628 | 2.677e+05 | 3.707e+07 | watchdog | 58.5 | no finished run |
| genome_100m | D58 | 58 (1024) | 0 / 1 | > watchdog | - | - | - | - | - | 0.409 | 1.631e+05 | 2.402e+07 | watchdog | 36.4 | no finished run |


**Reading, dna_100m.**

* **A true effective 0.5x was reached, in variant D:** D53 gives an effective cache of **0.49 [0.46-0.51]** and a whole-pipeline slowdown of **3.60x [3.59-3.63]** (98.3 s against 27.3 s; PH2/PH3 alone 5.0x: 87.5 s against
  17.6 s; PH1 is unaffected, 1.7 s, and PH4 + PH5 are 8.6 s because their limit was raised). **Gate (slowdown <= 3x at an effective 0.5x for dna_100m): not met, 3.6x.** Block-layer reads of the process 44.6 GB (the data is
  100 MB: every miss is a random 4 KB read and the rank file and text are read again and again), `memory.failcnt` 4.7 million, RssAnon 49 MB, no run killed, outputs identical in all 3 runs.
* **The knee:** the slowdown against the effective cache (calibration and sweep, variant D, against the unlimited 27.3 s, each 1 run unless noted): 0.72 -> 1.23x, 0.64 -> 2.0x, 0.59 -> 2.2x, **0.55 -> 3.1x**, 0.49 -> 3.6x (3 runs), 0.41 -> 3.5x (the plateau is the
  device: the 100 MB of data are read in 4 KB pieces). **The 3x line is crossed at an effective cache of about 0.55 x (text + rank file)**, which is far from the abandon criterion of section 5 ("the build needs an effective
  page cache above 2x (text + rank file) to stay within 3x on dna"): **that rule does not fire.** It is also a smaller cache than E3's "documented floor", but E3's number was a budget (limit minus anonymous memory) of which
  the data got a part.
* **The whole pipeline in one constant limit (variant C) does not get there:** C250 (effective 0.97) 1.39x [1.29-1.53], C400 1.39x (the unlimited cell has no reclaim; PH4 + PH5 take 15-22 s against 7.8 s because the 1.2 GB of
  scratch must be written back inside the group), **C80 (effective 0.78) 2.6x**, and **below an effective 0.7 PH4 + PH5 collapse: C65 (effective 0.65) did not finish in 15 minutes** (PH3 was done after 36 s, PH4 + PH5 had
  not finished after 865 s, `memory.failcnt` 179 million: the file budget of about 18 MB cannot hold PH4's 375 scratch files being written and read), the calibration run at 65 MB was stopped at 150 s and the one at 55 MB
  was OOM-killed (calibration table). **So a true 0.5x for the literal whole pipeline in one constant limit cannot be reached with this design and this machine: the smallest effective ratio at which the whole pipeline in one
  constant limit finished is 0.76-0.78 (2.4-2.6x); the smallest overall with a finished run is 0.41 (variant D, PH4 + PH5 outside the limit).** The cause is not text + rank but PH4's write-heavy scratch, which needs
  its own cache (or disk bandwidth: E4 saw 2.8-4.6x on PH4 alone with the cache nearly removed).

**Reading, genome-like 100 M** (text 37.5 MB + rank file 12.9 MB = 50.4 MB; the unlimited cell is noisy, 38.9-48.5 s, PH4 + PH5 14-24 s against 7.2 s in the scaling runs: its median 41.9 s makes the slowdowns
**understated by up to 1.2-1.3x**, which is why the PH2/PH3 slowdown, which has no such noise (20.4 s), is given next to it):

* **At an effective cache of 0.96-0.99 the whole pipeline runs at the unlimited speed**: C250 (0.97) 0.97x [0.96-1.03], D120 (0.96) 0.82x, C400 (0.99) 1.23x; PH2/PH3 1.1x. **The documented floor of about 1x (text + rank) holds for the
  whole pipeline including PH1** (PH1 does not slow down at any limit: 2.1-2.5 s in the limited cells against 3.7 s in the unlimited one).
* **Below it the cliff of E3 is there, softer:** effective 0.81 (D100) 2.8x (PH2/PH3 5.2x), 0.65-0.67 (D83, D90) **5.2-5.3x** (PH2/PH3 10x), **0.51 [0.508-0.514] (D66, 3 runs): 8.5x [8.3-8.6]** (356 s [349-359]; PH2/PH3 16.9x; 453 GB of
  block-layer reads for 50 MB of data; `memory.failcnt` 59 million), none killed within 15 minutes, outputs identical. The slowdown is not a function of the limit: **the same limit (90 MB) gave an effective 0.84 and 2.3x in the calibration and 0.67 and
  5.2x in the sweep** (once the rank file starts to be evicted the lookups of tied suffixes keep re-reading it); the effective cache, not the limit, predicts the time. **Genome-like text needs an effective cache of
  about 0.8-0.85 or more to stay at 3x (2.8x at 0.81, 3.4x against the 34.3 s of the non-binding D120 cell; 2.3x at 0.84 in the calibration), against 0.55 for random DNA.** For the abandon rule of 4.6 / 5 ("above 2x text + rank to stay within 3x on dna") the genome-like text is not the reference input and the
  need (0.8-1x) is far below 2x.

**What this says about the design.** (1) The cache the construction needs for text + rank is between 0.55x (random DNA) and about 0.85x (tie-rich DNA) of their size; the "1x" floor of 6b is a safe upper bound, and above
2x nothing is gained. At 3x10^9 text + rank are 1.1 GB (random DNA, 2 bits) to 1.5 GB (3 bits), times 0.55-1 = 0.6-1.5 GB (ESTIMATE by size; not measured). (2) **The part that was not in the E3 picture is PH4's scratch (12 B/char):** under a
tight total limit it, not the rank file, is what stops the build, so a deployment that really limits the page cache needs a few hundred MB beyond text + rank, or the bounded mode with its own 4 B/char PLCP file (E4). (3) The
anonymous memory of the design (43-68 MB in the lean settings against 105 MB) is not a fixed cost of the runtime: buffers per chunk and per bucket decide it and can be turned down for 12 % of the time.

### Calibration runs (kept; watchdog 150-300 s, stock settings C0, lean settings C1-C3 with 1 B2 worker, the others lean with 2 workers)

The limits of the sweep were found with these runs (`e7_cache_calib.tsv`; L = limit without the runtime pre-read, P = with it, D = limited until PH3 ends; `calib C0` = E6 settings and 96 MiB GC limit, unlimited; `C1`-`C3` = lean settings with 1 B2 worker, unlimited; slowdowns against the sweep's unlimited cell):

| input | cell | limit MB (raised to MB after PH3) | runs finished / started | total s | slowdown (median) | PH1 s | PH2/PH3 s | PH2/PH3 slowdown | PH4+PH5 s | effective cache (PH3 resident fraction of text + rank) | block-layer read MB (whole process) | memory.failcnt | killed | RssAnon peak MB | identity |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dna_100m | calib C0 | unlimited | 1 / 1 | 24.09 | 0.881 | 1.71 | 14.4 | 0.817 | 7.44 | 1 | 137.4 | 0 | no | 105 | 1 of 1 runs equal to the references |
| dna_100m | calib C1 | unlimited | 1 / 1 | 30.37 | 1.11 | 1.68 | 16.9 | 0.962 | 11.3 | 1 | 137.5 | 0 | no | 38.9 | 1 of 1 runs equal to the references |
| dna_100m | calib C2 | unlimited | 1 / 1 | 30.81 | 1.13 | 1.69 | 17 | 0.964 | 11.6 | 1 | 137.6 | 0 | no | 39.2 | 1 of 1 runs equal to the references |
| dna_100m | calib C3 | unlimited | 1 / 1 | 30.84 | 1.13 | 1.8 | 17 | 0.968 | 11.5 | 1 | 137.6 | 0 | no | 43.7 | 1 of 1 runs equal to the references |
| dna_100m | calib L160 | 160 | 1 / 1 | 37.64 | 1.38 | 1.68 | 17.5 | 0.996 | 17.9 | 0.917 | 3084 | 3.738e+05 | no | 49.7 | 1 of 1 runs equal to the references |
| dna_100m | calib L120 | 120 | 1 / 1 | 53.94 | 1.97 | 1.66 | 19.6 | 1.11 | 32.2 | 0.869 | 3855 | 4.363e+05 | no | 54 | 1 of 1 runs equal to the references |
| dna_100m | calib L100 | 100 | 1 / 1 | 49.93 | 1.83 | 1.66 | 21.1 | 1.2 | 26.6 | 0.849 | 4216 | 5.032e+05 | no | 40.7 | 1 of 1 runs equal to the references |
| dna_100m | calib L80 | 80 | 0 / 1 | > watchdog | - | - | - | - | - | 0.711 | 2.997e+05 | 4.953e+07 | watchdog | 56.2 | no finished run |
| dna_100m | calib P100 | 100 | 1 / 1 | 52.6 | 1.92 | 1.64 | 20.2 | 1.15 | 30.3 | 0.855 | 3834 | 4.285e+05 | no | 48.3 | 1 of 1 runs equal to the references |
| dna_100m | calib P80 | 80 | 1 / 1 | 65.15 | 2.38 | 1.64 | 24 | 1.36 | 39 | 0.76 | 4397 | 5.759e+05 | no | 42.4 | 1 of 1 runs equal to the references |
| dna_100m | calib P65 | 65 | 0 / 1 | > watchdog | - | - | - | - | - | 0.603 | 1.286e+05 | 1.997e+07 | watchdog | 47.8 | no finished run |
| dna_100m | calib P55 | 55 | 0 / 1 | > watchdog | - | - | - | - | - | 0.513 | 7.861e+04 | 1.053e+07 | oom-or-signal | 50.3 | no finished run |
| dna_100m | calib D70 | 70 (1024) | 1 / 1 | 36.41 | 1.33 | 1.94 | 25 | 1.42 | 8.9 | 0.7 | 3218 | 3.808e+05 | no | 50.7 | 1 of 1 runs equal to the references |
| dna_100m | calib D65 | 65 (1024) | 1 / 1 | 55.35 | 2.02 | 1.81 | 44.1 | 2.51 | 8.93 | 0.636 | 9096 | 1.162e+06 | no | 44.4 | 1 of 1 runs equal to the references |
| dna_100m | calib D60 | 60 (1024) | 1 / 1 | 60.33 | 2.21 | 1.82 | 48.7 | 2.77 | 9.3 | 0.593 | 1.364e+04 | 1.621e+06 | no | 43.3 | 1 of 1 runs equal to the references |
| dna_100m | calib P250 | 250 | 1 / 1 | 38.34 | 1.4 | 1.71 | 18.1 | 1.03 | 18 | 0.966 | 2854 | 3.142e+05 | no | 51.5 | 1 of 1 runs equal to the references |
| dna_100m | calib D56 | 56 (1024) | 1 / 1 | 84.14 | 3.08 | 1.75 | 73.2 | 4.16 | 8.64 | 0.546 | 2.678e+04 | 3.039e+06 | no | 56.4 | 1 of 1 runs equal to the references |
| dna_100m | calib D53 | 53 (1024) | 1 / 1 | 94.56 | 3.46 | 1.73 | 83.5 | 4.74 | 8.81 | 0.503 | 4.225e+04 | 4.519e+06 | no | 51.1 | 1 of 1 runs equal to the references |
| dna_100m | calib D50 | 50 (1024) | 1 / 1 | 94.55 | 3.46 | 1.7 | 83.4 | 4.74 | 8.96 | 0.413 | 5.293e+04 | 5.459e+06 | no | 52.6 | 1 of 1 runs equal to the references |

The genome-like calibration cells (D90, D75, D66, D58 with a 240 s watchdog) are labelled `C.<cell>` in the same file:

| input | cell | limit MB (raised to MB after PH3) | runs finished / started | total s | slowdown (median) | PH1 s | PH2/PH3 s | PH2/PH3 slowdown | PH4+PH5 s | effective cache (PH3 resident fraction of text + rank) | block-layer read MB (whole process) | memory.failcnt | killed | RssAnon peak MB | identity |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| genome_100m | D90 | 90 (1024) | 1 / 1 | 95.83 | - | 2.14 | 83.7 | - | 9.47 | 0.839 | 4.206e+04 | 6.428e+06 | no | 53.6 | 1 of 1 runs equal to the references |
| genome_100m | D75 | 75 (1024) | 0 / 1 | > watchdog | - | - | - | - | - | 0.628 | 2.677e+05 | 3.707e+07 | watchdog | 58.5 | no finished run |
| genome_100m | D66 | 66 (1024) | 0 / 1 | > watchdog | - | - | - | - | - | 0.544 | 2.974e+05 | 3.98e+07 | watchdog | 37.5 | no finished run |
| genome_100m | D58 | 58 (1024) | 0 / 1 | > watchdog | - | - | - | - | - | 0.409 | 1.631e+05 | 2.402e+07 | watchdog | 36.4 | no finished run |


## 9. Defects found

* **No defect in the product code.** The one change of this prototype is the option (section 2).
* **Measuring script (found by the first cache sweep, fixed, no data lost):** the watchdog of `run_e7_cache.sh` was a subshell around `sleep 900` that kept the stdout pipe open after the run ended, so the sweep driver waited 15 minutes
  for the pipe after the first cell; the rows of the unfinished second run were purged (the protocol script purges the rows of an interrupted run) and the sweep restarted; the first unlimited run, already complete, was kept. The
  watchdog now polls once a second and ends with the run; the same change was made in `run_e7_protocol.sh` (the scaling runs, made before the change, used the `sleep` version with a 45-minute limit, which never fired). The first `run_e7_setup.sh` made all hashes at the end; it was rewritten with a hash per step
  and step markers when the 500 M disk plan needed the SA reference deleted early.
* **Method note:** the E3-style "budget" (limit minus anonymous memory over text + rank, `eff_budget_ph3` in the rows) is a poor predictor of the effective cache (3.7 at an effective 0.92, 1.4 at 0.76): the other files of the group
  take their share; the measured residency is the figure used.

## 10. Not measured, not verified, open points

* **Not measured:** anything above 500 M (3x10^9 figures are ESTIMATES); genome-like series A at 500 M and genome-like 1 / 2 threads at 500 M; a^n or any chain at E7 (E6 did a^n 50 M; E7 did not ask); a real disk (virtio, probably
  host-cached); the page cache really exceeded at a large n (it was not: 15 GB of RAM); resume, 40-bit positions, PH6, PH7; the cache test at a size other than 100 M (the scale-down parameters of 3x10^9, not real block sizes), with
  other kernels or with cgroup v2; more than one run in most cache cells (3 at the unlimited, 1x and 0.5x cells); the effective cache of PH1 (only PH3 is weighted; PH1 did not slow down at any limit).
* **Not verified:** the cause of the series A growth is read from the code and the arena rows (buffers per chunk), the cap-fires-at-P-2,400 statement is arithmetic, not a run; the 3x10^9 estimates; why PH4 + PH5 of the genome-like
  unlimited cell are 14-24 s while the scaling runs of the same input take 7.2 s (an unlimited sub-group on a 15 GB machine: dirty write-back of PH3 competing with PH4; not investigated); the reference code is the in-repo
  in-RAM code (not an independent implementation; the 500 M `r500` check adds the E3 builder's SA).
* **Tests:** EmProto default build **341 of 341** (334 + 7 new cases; 10 m 29 s); the checked build was not run (the change is option plumbing in the driver, no code that the checked build compiles differently); no mutation audit (no new logic beyond
  `ruleN = RuleN > 0 ? RuleN : n`, covered by 8 tests, 6 of them red against the stub).
* **Open points for 1b.7 / E8 / E9:** (a) PH4's buffers are `P x 32 KB` per writer: with R fixed, memory grows with n (series A), so R (or the buffer) must follow n or a second distribution level is needed above about 600 M at
  R = 266,667, and the real R = 8,000,000 gives P = 375 only up to 3x10^9; (b) PH1's number of sample buckets must stay fixed (K_s proportional to n, as the rule does), or PH1 is quadratic (b = 1.66); (c) PH4 + PH5 are
  the most cache-hungry part of a limited deployment (12 B/char scratch): decide between the base and the bounded mode with the cache in mind; (d) the lean settings (small buffers, 48 MiB GC limit) are available as options and
  cost 12 %; (e) E8 (1 Gbp, 3 Gbp) should be run with the page cache really limited, because here it never was at the scaling sizes; (f) the 3x10^9 build needs roughly 0.6-1.5 GB of cache for text + rank (ESTIMATE) plus the scratch's.

## 11. Commits, files, how to repeat

Commits (all `ST-S1: ADR-004 E7 ...`): `ba1fbc0e` (the option, tests first), `8bdd6c27`, `b086a1a8`, `6644c6f9`, `16aa71be` (scripts and raw rows) and the commit of this document, ADR-004 section 5a and PROGRESS-S1 after them. Work directory `/tmp/stbench/e7` (texts, `*.refs.sha`, the unfinished scratch is deleted by the scripts; the outputs of every run are deleted after their hashes), binary
`/tmp/stbench/bin_e7_run` (`dotnet SuffixTree.EmProto.dll e6 ...`), references of the 500 M inputs deleted after hashing (rebuild with `run_e7_setup.sh`). Repeat: `run_e7_setup.sh <dir> <rows> dna_50m ...`, `run_e7_protocol.sh <dir> <rows> AB dna_50m dna_100m dna_200m`
(and `B` / `A` / `threads` for the other groups), `run_e7_cache_sweep.sh <dir> <rows> dna_100m U:max:-:900:3 C250:250:-:900:3 D53:53:1024:900:3 ...`, then `python3 -I tables_e7.py fit|rss|io|threads|cmp|cache <rows...>`.
