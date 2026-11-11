# E3b: PH1 (the difference-cover rank file) in external memory - measured results

Work package of the owner decision of 2026-11-11 (ADR-004 section 6c: the "external-memory PH1" prototype that 6b called E4). Code: `apps/SuffixTree.EmProto` (`Dcs/ExternalDcs.cs`, `SampleNamer.cs`, `SampleGeometry.cs`, `Runs.cs`, `NameSinks.cs`, `SaInverter.cs`, `MappedWindowWriter.cs`; command `e3b` and `packrand` in `Experiments/E3b.cs`), tests `tests/SuffixTree/SuffixTree.EmProto.Tests/Ph1`.
Raw rows: `e3b_*.tsv` in this directory (`R <TAB> E3b <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes`), scripts `run_e3b_protocol.sh`, `summarize_e3b.py`, `tables_e3b.py`. Everything below is **measured here** (4 vCPU Intel Xeon, 15 GB, kernel 6.18.44, .NET 10.0.14, Release,
workstation non-concurrent GC, nothing else running, build servers shut down, `DOTNET_GCHeapHardLimit` = 96 MiB, every timed run in a fresh process after `sync; echo 3 > /proc/sys/vm/drop_caches`, i.e. a **cold page cache**) unless it is marked ESTIMATE.
The oracle of every comparison is the in-RAM `DcsBuilder` (`mkranks`), the table of E2/E3; the output of the file builder is compared with it byte for byte (`--verify`).

## 1. What PH1 does here

For a text T_l (level 0: the packed input, read through a read-only mapping) with the period X (Singer 993, |D| = 32 for the nominal 1024), sample positions `p mod X in D`, `Samples = 3.2 %` of n:

1. **Naming (ADR 4.2, "buckets of the sample").** The sample is sorted by its X-prefix in sample buckets; one task per bucket scans all sample positions in text order (a 64-bit key per position, read through the mapping), keeps the records of its range `[lo, hi)`, sorts them in RAM (LSD radix on the key, merge sort with the X-prefix comparator on key ties), names equal X-prefixes (names are 1, 2, ... in sorted order) and writes a **run** of `(slot, local name)` records sorted by slot (a scratch file, or memory when tiny). The bucket order is `(64-bit key, X-prefix, position)`, a strict total order, so equal prefixes (a^n, satellites, N runs) are spread over buckets like everything else and no bucket ever exceeds K_s. A bucket larger than K_s (the root, or a rare overflow) is split by one more scan that counts its records and keeps a uniform reservoir; the quantiles of the reservoir are the splitters of the children. Groups of equal prefixes that meet at a bucket boundary share one name (a boundary check after all buckets are known).
2. **Merge back into position order.** The runs are scattered, window by window of R slots, into an array of R entries (`rank = global name`); each run is read through a bounded buffer (all buffers share one budget). No heap, no comparison; slots are a permutation of `0..Samples-1`.
3. **All names distinct:** the names are the ranks; the window goes to the rank file through bounded mapped windows (1 MiB views, no explicit flush).
4. **Otherwise** the names form the reduced string (class by class, a 0 after every class; `Samples + |D|` symbols over a 32-bit alphabet): written as a big-endian 32-bit file by |D| sequential cursors (a transposition without an intermediate file). At most N_base symbols: SA-IS in RAM, inverse, rank file. Longer: the file is mapped as a text of 32 bits per symbol, **the same builder runs one level down** to give that text its own rank file, the suffix array of the reduced string is built by the bucket sorter of E3 (`EmBuilder`, unchanged) with that rank file, and the inverse permutation goes through a distribution (tuples `(slot, rank)` into P = Slots/R chunk files, then chunk by chunk into the rank file of this level).
5. The rank of a sample suffix is its index in the suffix array of the reduced string, separators included, as in `DcsBuilder`; the unused slots of the last block are zero. Output independent of thread count, K_s, N_base, R, buffer sizes and seed (tested).

## 2. Design choices the ADR does not fix (reason, alternative)

| choice | reason | alternative not taken |
|---|---|---|
| Bucket order `(key, X-prefix, position)`, buckets are ranges between two sample records | strict total order: bucket size bounded by K_s for every text without special cases (a^n would put all 3x10^9 x 3.2 % records into one equality bucket otherwise); the boundary check for shared names is 20 lines | samplesort "equality buckets" (a separate bucket per splitter value, constant name); sorting only by key and recursing on 64-bit pieces (X/21 rounds on a^n) |
| One task per bucket, **scan of all sample positions per bucket** (ADR 4.2 literally), `2 x Samples / K_s` = 185 scans at 3x10^9 | no position files, no partition passes (PH1 I/O 0.65 B/char as in the ADR table); parallel by bucket; 185 x 96.7 M record visits = 36 s on 4 threads | partition into position files first (E3 style, +0.26 B/char written and read twice, fewer text scans) |
| Splitters from a uniform reservoir of the bucket taken in a counting scan, fan-out `ceil(2 count / K_s)` (target bucket K_s/2) | exact count known before the split; an overflow of a bucket (about nine standard deviations with 16,384 reservoir records for 185 buckets; never seen in the measured runs: scans = leaves + 1 everywhere) costs one more scan and nothing else; the same code handles the root and the overflow | stratified random positions for the root (a second code path), an exact counting pass |
| Names back into position order by **scatter windows** of R slots | O(Samples) with no heap or compare; the buckets' runs are already sorted by slot, so each run is read once, sequentially | k-way merge with a heap (log 185 comparisons per record), DIST of `(slot, name)` (one more write and read of 8 B per sample) |
| Reduced string written class-major by one cursor buffer per class (32 of them) | the transposition costs no extra file and no pass | block-major file plus a transposing pass |
| Recursion: same X at every level; the reduced string is a **32 bits per symbol big-endian file = `PackedText` with b = 32** | level l+1 reuses the text, key, LCE and bucket sorter of E3 unchanged; 2 symbols per 64-bit key (names are rich); X_l = X is what `DcsBuilder` does | own X_l per level (shorter at deep levels), b = bits of the name count (a bit-packing writer) |
| Level >= 1 suffix array by `EmBuilder` (E3), inverse by one distribution straight into the slot-ordered rank file | no second implementation of PH2-PH3; the DIST of ADR 3.3 | a separate ISA file and a transposing pass |
| Base case by `SaisReference<int>` (safe SA-IS), not `SaisKernel` | `SaisKernel` has only byte and 16-bit entries (its int recursion is private); `DcsBuilder` uses the same routine; measured below | an int32 entry in `src/**` (not allowed in this package) |
| N_base reference 16 Mi (the ADR value) | measured full-size: 150.8 MB arena, **RssAnon 180.7 MB** [180.70-180.90] at 520 Mbp a^n; level structure at 3x10^9 is identical for any N_base from 3.2 M to 16 M (L1 = 96.7 M is external, its reduced string 3.1 M is the base case) | 4 Mi: 64 MB at the base case (measured at 130 Mbp a^n), but a further external level above about 4x10^9 symbols (both strands) |
| K_s reference 1 Mi records: 24 B per record per worker, 100.7 MB for 4 workers | fits the 160 MiB arena cap with the merge and the base case phases separate | 512 Ki (halves the arena, doubles the scans: 370 instead of 185) |
| One shared budget (16 MiB) for the merge's run buffers | the only term that grew with n / K_s (one 64 KB buffer per run); capped | 64 KB per run (11.8 MB at 3x10^9, 24 MB at 6x10^9) |
| `--touch-buffers` | RssAnon of a scale-down run counts only touched pages; with this flag every arena block is made resident at its allocation, the worst case used for the slope of gate G3 | none |

## 3. Verification

| item | result |
|---|---|
| EmProto tests before the work | 58, all green |
| tests first, against the stub (`26109a36`) | **17 of 75 red** (`NotImplementedException`), the 58 older green; the explicit fuzz class is not in the count |
| after the implementation | 75 green; with the component tests written after the code, the option and boundary tests, **106 green** (default run, 48 new, about 50 s) and **106 green in the checked variant** (`-p:EsaBounds=true`, `ESA_EXPECT_BOUNDS=true`; the project has no other configuration) |
| exhaustive (oracle class, in the default run) | every binary string of length 0-11 with X = 3, 0-12 with X = 7 (external levels at N_base = 3), 0-10 with X = 9 (2 bits per symbol), every ternary string of length 0-6 with X = 5: **15,426 strings**, rank file equal to `DcsBuilder` every time |
| differential, adversarial families at n = 1,000 and 5,000 (30 families, 2 random parameter sets each: K = 3..50, X = 7..64, b in {min,2,3,4,8,16,32}, threads 1/2/4, N_base 4..200, F 2..6, reservoir 8..208, R 4..104, spill, buffer, window, merge budget, K_sort 3..50) | >100 cases equal; asserted: most have several buckets and splits, >20 have a recursion level, >5 reach level 2, >20 end in the base case, >20 use run files; the level-0 distinct-name and shared-name counts equal a naive dictionary count in every case |
| `Fuzz` class (explicit, `EMPROTO_FUZZ_SECONDS=120`, seed 1; texts 10^3..10^5, same parameter ranges) | **632 cases, 430 with at least one recursion level, all equal**, 120 s (an exploratory run of 60 s with seed 77 before the measurements: 219 cases) |
| written after the code (labelled in the test file) | `SampleGeometry`, `MappedWindowWriter`, `BucketRun`/`RunReader`, `WorkPool`, the three name sinks, `SaInverter` against a naive inverse, `PrefixOrder` against a naive comparison, `packrand`, the touch mode, the `e3b` pipeline and `CompareFiles`, failure cleanup; later added after the first measurements: the merge-buffer budget, the base-case boundaries, buffer sizing, the parameters row |
| one expectation of the first test version was wrong and corrected after the first run | `RecursionDepth_OfAllSame...` expected a base case at the last level of a^20000 (X = 31, N_base = 50); the chain actually ends one level earlier with unique names (levels 0, 1, 2 of 20000, 3878, 758 symbols); the coverage thresholds of the exhaustive configurations were set from a first run (X = 7 strings of length <= 12 have only 8 non-unique cases) and are lower bounds, not tuned to pass |

## 4. Gate G1: byte-identical rank file

Every run of the protocol compares the output with the rank file of the in-RAM `DcsBuilder` (`--verify`, byte for byte, 1 MB pieces): **188 of 188 runs equal, none different** (`rank_equals_oracle` in the rows), cold cache, fresh process each. Coverage:
* **threads 1, 2 and 4** on every input of the protocol (dna_100m, genome-like 100 M, dna_10m, rep_10m, a^n 10 M, Fibonacci 10 M, N run, genome-like 10 M, 171-bp tandem x 10^4, prot_1m, dna_1m and the committed small adversarial families a^n, Fibonacci, Thue-Morse at 10^5 and 10^6, a 171-bp tandem of 10^5 and an N run of 10^6);
* **three bucket capacities** (scaled K_s / 2, K_s, 2 K_s = 17,476 / 34,953 / 69,905 at 100 M, giving 369 / 185 / 93 leaf buckets) on dna_100m, genome-like 100 M and all ten families;
* sizes 10 M, 50 M, 100 M, 200 M (random DNA; genome-like to 100 M) and the unscaled real-size runs of section 7 (3x10^9 random DNA, 2x10^9 genome-like, 10^9 a^n, 520 M a^n): equal.

**G1: PASS.**

## 5. Gate G2: RssAnon <= 256 MB at 100 Mbp, including the runtime and the thread stacks, 4 threads

Scale-down emulation (section 6), `DOTNET_GCHeapHardLimit` = 96 MiB, sampler = `/proc/self/status` every 10 ms (13 threads at the peak):

| input | RssAnon peak, 3 runs of the baseline: median [min-max] | all runs of the slope series at 100 M (8 runs) | cap |
|---|---|---|---|
| dna_100m (random DNA) | **21.1 MB** [20.6-21.2] | 21.2 [20.5-22.6] | 256 |
| genome-like 100 M (0.84 % shared names, level 1) | **25.2 MB** [24.8-28.3] | 29.5 [27.7-32.7] | 256 |

The worst RssAnon of any of the 176 scale-down runs (all groups) is 32.7 MB. Mapped file pages (RssFile 54 and 79 MB) are not counted (ADR-004 section 6a). The same builder at the **real, unscaled parameters** (buckets of 1 Mi records, N_base 16 Mi): 3x10^9 random DNA **78.9 MB** [78.8-79.5], 2x10^9 genome-like 149.8 MB [149.5-149.9], 10^9 a^n 99.3 MB [98.6-101.7], and the worst phase, an in-RAM base case of 16.76 M reduced symbols (520 M a^n), **180.7 MB** [180.7-180.9]; all under 256 MB (section 7).

**G2: PASS** (21-33 MB at 100 Mbp; at most 181 MB at any measured real size).

## 6. Gate G3: scale-down emulation and the slope of RssAnon against n

**The 3x10^9 reference and how it is scaled** (`--scale-n 3000000000`; `E3b.KsRef` and friends): X = 993 (Singer, |D| = 32), so a text of 3x10^9 symbols has 96,676,740 sample positions. Every quantity that decides the number of buckets, levels, windows or runs is multiplied by n / 3x10^9:

| quantity | at 3x10^9 | derivation | 10 M | 100 M | 200 M |
|---|---|---|---|---|---|
| K_s, records per sample bucket | 1,048,576 | 4 workers x 24 B x K_s = 100.7 MB of the 160 MiB arena cap | 3,495 | 34,953 | 69,905 |
| sample buckets (leaves) | **185** = 2 x 96,676,740 / K_s | root split into `ceil(2 count / K_s)` buckets; measured 185 leaves at 3x10^9 unscaled | 185 | 185 | 185 |
| N_base | 16,777,216 | the ADR value | 55,924 | 559,241 | 1,118,481 |
| recursion | level 1 = 96,676,772 symbols > N_base: external; its reduced string (about 3.1 M) < N_base: base case. a^n: levels 0 -> 1 -> base | Samples + classes per level | the same structure (levels = 1 for every non-unique input, 0 for random DNA) | same | n/a |
| R, slots per merge window / chunk | 4,194,304 | `ceil(96.7 M / R)` = 24 windows | 13,981 | 139,810 | 279,620 |
| K_sort, leaf of the external suffix sort | 1,048,576 | as K_s | 3,495 | 34,953 | 69,905 |
| run spill threshold | 65,536 B | so that every bucket run is a scratch file as at 3x10^9 (4 MB runs) | 218 | 2,185 | 4,369 |
| fixed buffers (not scaled) | 64 KiB I/O buffers, 1 MiB windows, 16 MiB merge budget | they do not depend on n | same | same | same |

**Slope of RssAnon against n.** Least-squares fit of the RssAnon peak over log2 n (the median of each size; the fit through every single run is in brackets), scaled emulation, 4 threads, fresh process and cold cache each, 8 runs per size in the plain series and 3 in the touched series:


| series | RssAnon peak per n (median [min-max], runs per n in the text) | least-squares slope MB per doubling | step slopes MB per doubling |
|---|---|---|---|
| random DNA, plain buffers, 10-200 M | 10 M: 10.7 MB [7.9-11.5] (8 runs) ; 50 M: 21.1 MB [20.1-21.3] (8 runs) ; 100 M: 21.2 MB [20.5-22.6] (8 runs) ; 200 M: 22.6 MB [22.1-23.2] (8 runs) | 2.79 (all runs: 2.88) | 4.45 / 0.15 / 1.35 |
| random DNA, plain buffers, 10-100 M | 10 M: 10.7 MB [7.9-11.5] (8 runs) ; 50 M: 21.1 MB [20.1-21.3] (8 runs) ; 100 M: 21.2 MB [20.5-22.6] (8 runs) | 3.39 (all runs: 3.50) | 4.45 / 0.15 |
| genome-like, plain buffers, 10-100 M | 10 M: 19.0 MB [18.0-19.8] (8 runs) ; 50 M: 24.3 MB [23.3-26.5] (8 runs) ; 100 M: 29.5 MB [27.7-32.7] (8 runs) | 3.00 (all runs: 3.18) | 2.27 / 5.19 |
| random DNA, all buffers touched, 10-200 M | 10 M: 20.1 MB [20.0-20.2] (3 runs) ; 50 M: 21.1 MB [20.3-22.8] (3 runs) ; 100 M: 21.1 MB [20.9-21.1] (3 runs) ; 200 M: 22.2 MB [22.1-23.1] (3 runs) | 0.44 (all runs: 0.48) | 0.42 / 0.05 / 1.08 |
| random DNA, all buffers touched, 10-100 M | 10 M: 20.1 MB [20.0-20.2] (3 runs) ; 50 M: 21.1 MB [20.3-22.8] (3 runs) ; 100 M: 21.1 MB [20.9-21.1] (3 runs) | 0.33 (all runs: 0.33) | 0.42 / 0.05 |
| genome-like, all buffers touched, 10-100 M | 10 M: 22.3 MB [21.8-27.0] (3 runs) ; 50 M: 24.6 MB [24.1-27.1] (3 runs) ; 100 M: 25.2 MB [24.5-25.8] (3 runs) | 0.89 (all runs: 0.49) | 0.98 / 0.62 |

* "plain": as run. The step 10 M -> 50 M of random DNA is 4.45 MB per doubling and the step 50 M -> 100 M of the genome-like text 5.19 MB per doubling: the first because the 64 KiB run buffers of the merge are only partly filled when a run is shorter than a buffer (a run is 14 KB at 10 M, 4 MB at 3x10^9), the second is within the run-to-run spread of the genome-like text (27.7-32.7 MB at 100 M, the garbage of the managed arrays of the external level). The fit over the sizes, which is what the gate asks for, is 3.4 / 3.0 MB per doubling (10-100 M), 2.8 with 200 M.
* "touched" (`--touch-buffers`: every arena block resident at its allocation, so the buffers are fully counted at every size): 0.33 / 0.89 MB per doubling (10-100 M), 0.44 with 200 M.
* What grows with n in the emulation is the scaled part: the worker arrays (4 x 24 B x K_s: 3.4 MB at 100 M, 6.7 MB at 200 M); they are not the maximum over the phases until n is about 400 M, because the merge phase (run buffers 11.8 MB + window) is larger. Terms that grow with the number of buckets are capped (merge budget 16 MiB, tuple buffers P x 64 KiB with P = Slots / R, 24 at 3x10^9).

**G3 (emulation: bucket and level counts of 3x10^9, slope <= 5 MB per doubling): PASS** (slope 3.4 / 3.0 plain, 0.33 / 0.89 touched; two single steps of the plain series exceed 5 MB per doubling for the reasons above and are reported, the gate is the fit).

**ESTIMATE of RssAnon at 3x10^9, labelled so.** The estimate does not rest on the slope of the emulation but on a phase model checked against real-size runs: the maximum over the phases of (runtime and threads, about 8 MB idle) + (the arena of the phase) + (managed arrays), where the phases are: naming, 4 x 24 MiB = 100.7 MB if completely touched (measured 78.9 MB with the buckets half full: 3x10^9 random DNA, unscaled, **measured, not estimated**); merge (about 29 MB); the external suffix sort of a reduced string (EmBuilder, K_sort 1 Mi: **149.8 MB** measured on the 2x10^9 genome-like text, whose reduced string of 64.5 M symbols is two thirds of the 96.7 M of a 3x10^9 text); the in-RAM base case (**180.7 MB** measured at the full N_base of 16.8 M symbols, a case the genome-like text does not reach: its base case at 3x10^9 is 3.1 M symbols). So: **3x10^9 random DNA 79 MB (measured); 3x10^9 genome-like or repetitive text ESTIMATE 150-185 MB, below the 256 MB cap by at least 70 MB.** It rests on: the arena cap (160 MiB, enforced from inside; peak 100.7 / 150.8 MB), the n-independence of every other buffer, and the three real-size runs; it does not cover a text on which the external level has more than one level below it (more than about 4x10^9 symbols at N_base 16 Mi) or 40-bit positions.

## 7. G4 report (no gate)

Scale-down runs, 4 threads, cold cache (explicit I/O = `pread`/`pwrite` of runs, level files and tuple files counted by `CountedFile`; mapped = writes through the rank-file windows; block-layer read = `/proc/self/io read_bytes`, i.e. the cold text read):


| input | n | runs | PH1 s (4 threads) | ns/symbol | RssAnon MB | levels | shared names at level 0 | explicit I/O B/char (read+write) | mapped write B/char | block-layer read B/char (cold) | peak scratch B/char |
|---|---|---|---|---|---|---|---|---|---|---|---|
| dna_100m (random) | 100,000,000 | 3 | 1.59 [1.57-1.61] | 15.9 [15.7-16.1] | 21.1 [20.6-21.2] | 0 | 0.000 % | 0.516 (0.258+0.258) | 0.129 | 0.274 | 0.258 |
| genome-like 100 M | 100,000,000 | 3 | 2.14 [2.10-2.19] | 21.4 [21.0-21.9] | 25.2 [24.8-28.3] | 1 | 0.844 % | 1.657 (0.828+0.828) | 0.133 | 0.402 | 0.392 |
| dna_10m | 10,000,000 | 3 | 0.29 [0.28-0.31] | 28.5 [28.4-30.9] | 11.4 [10.6-12.0] | 0 | 0.000 % | 0.516 (0.258+0.258) | 0.129 | 0.441 | 0.258 |
| rep_10m | 10,000,000 | 3 | 0.38 [0.32-0.39] | 37.7 [31.8-39.0] | 11.5 [9.2-11.8] | 0 | 0.000 % | 0.516 (0.258+0.258) | 0.129 | 0.499 | 0.258 |
| a^n 10 M | 10,000,000 | 3 | 1.80 [1.80-1.82] | 180.2 [179.6-181.9] | 17.7 [17.7-17.8] | 1 | 99.990 % | 1.092 (0.546+0.546) | 0.133 | 0.389 | 0.387 |
| Fibonacci 10 M | 10,000,000 | 3 | 0.99 [0.85-1.07] | 99.4 [85.0-106.6] | 17.8 [17.8-17.8] | 1 | 99.990 % | 1.092 (0.546+0.546) | 0.133 | 0.375 | 0.387 |
| N run 1 M in 10 M | 10,000,000 | 3 | 0.85 [0.68-1.23] | 85.0 [67.6-122.6] | 17.9 [17.6-18.7] | 1 | 9.990 % | 1.092 (0.546+0.546) | 0.133 | 0.632 | 0.387 |
| genome-like 10 M | 10,000,000 | 3 | 0.37 [0.36-0.45] | 37.0 [36.3-44.6] | 19.6 [17.2-19.9] | 1 | 0.791 % | 1.092 (0.546+0.546) | 0.133 | 0.612 | 0.387 |
| 171-bp tandem x 10^4 | 1,710,000 | 3 | 0.38 [0.38-0.41] | 224.4 [223.3-242.7] | 14.4 [14.3-14.4] | 1 | 99.942 % | 0.790 (0.395+0.395) | 0.133 | 1.638 | 0.387 |

Observations. (a) **I/O** of PH1 alone: 0.52 B/char explicit plus 0.13 B/char mapped = 0.65 B/char for random DNA (the ADR table: about 0.6), 1.1-1.8 B/char when the names are not unique and a level below exists (reduced string, its rank file, suffix array tuples), plus the cold read of the text (0.25-0.4 B/char) and, for the external level, the E3 position buckets (counted in the explicit figure) and the SA of the reduced string (0.13 B/char through mappings, not in these numbers); proposed limit for the whole build: 100 B/char. (b) **Peak scratch** 0.26 B/char (random) to 0.39 B/char; proposed limit for the whole build 8 B/char; the output rank file (0.129 B/char) is not scratch. (c) **Recursion depth**: 0 for random DNA, rep_10m and the committed random inputs (all names distinct); 1 (one external level whose reduced string is sorted in RAM, or an external level with unique names) for every other input, at every size; the repeated-names share at level 0 is 0 % (random, rep_10m), 0.79-0.85 % (genome-like), 9.99 % (N run of 1 M in 10 M), 94-100 % (a^n, Fibonacci, Thue-Morse, tandem). (d) **Time**: 16-21 ns/symbol at 100 Mbp (4 threads), 49-55 ns with one thread; adversarial 10 M inputs 85-224 ns/symbol; the in-RAM oracle `DcsBuilder` needs 4.4-13 ns/symbol but 1.4 GB (100 Mbp) to 3.1 GB (3x10^9) of RAM.

Real-size runs, **unscaled parameters** (K_s = 1 Mi, N_base = 16 Mi, R = 4 Mi, 4 threads, cold cache, 3 runs each, `e3b_extra.tsv`); these are runs of PH1 only, not of E8:

| text | n | PH1 s | ns/symbol | RssAnon MB | arena peak MB | levels (shared names at level 0) | explicit I/O B/char (r+w) + mapped | peak scratch B/char | equal to oracle | oracle (`mkranks`, in RAM) |
|---|---|---|---|---|---|---|---|---|---|---|
| random DNA (`packrand`) | 3,000,000,000 | **36.7** [36.1-37.2] | 12.2 | **78.9** [78.8-79.5] | 100.7 | 0 (0 %) | 0.516 + 0.129 | 0.258 | 3 of 3 | 14.9 s, 3.1 GB |
| genome-like | 2,000,000,000 | **28.1** [27.5-28.3] | 14.0 | **149.8** [149.5-149.9] | 100.7 | 1 (0.85 %; level 1 = 64.5 M symbols, external, names 0.55 % shared, base case below) | 0.851 + 0.851 + 0.133 | 0.396 | 3 of 3 | 27.4 s, 3.3 GB |
| a^n | 1,000,000,000 | **93.4** [90.3-96.6] | 93.4 | 99.3 [98.6-101.7] | 100.7 | 1 (100 %; level 1 = 32.2 M symbols, external) | 0.782 + 0.782 + 0.133 | 0.391 | 3 of 3 | 75.6 s, 1.3 GB |
| a^n (base-case probe, N_base = 16 Mi) | 520,000,000 | 6.46 [6.42-6.46] | 12.4 | **180.7** [180.7-180.9] | 150.8 | 0 (100 %; reduced string 16.76 M symbols sorted in RAM) | 0.258+0.258+0.129 | 0.258 | 3 of 3 | 37.4 s, 0.7 GB |

(`a^n 130 M` with N_base = 4 Mi, run once before the final binary: RssAnon 63.6 MB, byte-identical; kept as the alternative N_base of section 2.) The 3x10^9 random run reads the 750 MB text 185 times through the mapping; the cold block-layer read is 0.25 B/char = one pass because the page cache then holds it.

## 8. What was not verified, and the limits of the evidence

* **No page-cache limit was applied.** Every run starts cold (dropped cache) and then has all of the 15 GB for the page cache. PH1 scans the whole text once per bucket (185 times at 3x10^9): if the text (750 MB at 2 bits, 1.1 GB at 3 bits, 0.25-0.375 B/char) is not resident, each scan is a pass of block-layer reads (up to 185 x 750 MB = 139 GB for random DNA). This is the page-cache floor of ADR-004 section 6b, now measured for PH1 only as "cold, then cached"; the 0.5x cache test is E9-class work and was not done here.
* **Disk**: the virtio disk of this VM (probably host-cached); real NVMe or network latency is not represented. Scratch files of PH1 are 0.26-0.40 B/char, written once and read once.
* **Sizes**: the gates are at 100 Mbp. The real-size, unscaled runs reach 3x10^9 only for random DNA (names unique, no recursion); the genome-like text was run to 2x10^9 (the in-memory generator cannot exceed 2^31 symbols), a^n to 10^9 (520 Mbp and 130 Mbp as base-case probes); nothing above 3x10^9, no both-strand 6.2x10^9 text, no 40-bit positions (the prototype's positions are 32-bit, as in E3).
* **The oracle is `DcsBuilder`** (E2), itself pinned to the suffix array by tests, not an independent implementation of the sample ranks; for 3x10^9 it needs 3.1 GB of RAM (one run, 15 s).
* **Threads**: 1, 2 and 4 only (4 vCPU). The merge, the inverse and the base case are single-threaded (1.9 s of 36 s at 3x10^9 random).
* **`SaisReference<int>` instead of `SaisKernel`** in the base case; its speed was not compared (0.4-4.2 s for 4-16 M symbols in the probes, not a decision driver).
* **No crash safety / resume** (the state.json of ADR 4.3 is not part of this prototype); a failed build deletes its tagged scratch files.
* **The time figures** of the scale-down runs are those of a build whose buckets are 30-300 times smaller than the real ones (the scans cost the same per record, the per-bucket overhead is higher); the real-size runs are the reference for time.

## 9. Open points for E4 (PH4) and for 1b.7

* **Interface.** `ExternalDcs.BuildRankFile(PackedText text, DifferenceCover cover, string rankPath, ExtDcsOptions options)` produces the rank file in the layout `DcsRankTable`/`EmBuilder` already read (slot order, `uint` little endian); the E3 sorter consumed such a file written by this code at level >= 1 in every recursive test, so PH3 can take it unchanged. 1b.7d needs a `PackedText` over the TEXT section (PH0) and the parameter record.
* **Page-cache floor of PH1.** 185 scans of the sample keys = 185 passes over the text through the mapping: with less cache than the text (0.75-1.1 GB at 3x10^9) each scan is a block-layer pass (up to 139 GB). The E9-class 0.5x test has to include PH1. Possible remedy if it hurts: several adjacent buckets per scan (q-fold fewer scans for q-fold arena), at the price of `q x 100 MB`.
* **Base-case memory.** 181 MB at N_base = 16 Mi is the largest phase; N_base = 4 Mi costs 64 MB with the same level structure up to about 4x10^9 symbols (130 M a^n probe). The decision belongs to 1b.7.
* **`SaisKernel` has no int32 entry** (its recursion on `int` is private to `src/**`); the base case uses `SaisReference<int>` (safe, slower). An int32 public entry in `SuffixTree.Esa` would be a change outside this package.
* **For PH4 (E4):** `SaInverter` is the distribution-and-scatter primitive (tuples into P chunk files by key / R, then chunk by chunk into a mapped output), `Runs`/`MappedWindowWriter` the bounded writers; PH4's Phi tuples need the same with 12 B tuples and a second distribution; they can be generalised from here. `ScratchTracker`/`MappedIo`/`IoCounters` already cover scratch, mapped writes and explicit I/O of PH1 completely.
* **Not covered by this prototype**: resume after a crash (state.json), 40-bit positions / both strands, a thread-parallel merge and inversion (1 thread; 1.9 of 36 s at 3x10^9), the cgroup-limited cache.

## 10. Mutation audit of the new code

(see below; filled in after the audit)
