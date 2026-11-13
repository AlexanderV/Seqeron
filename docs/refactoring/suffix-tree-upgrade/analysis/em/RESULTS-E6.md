# RESULTS-E6: the end-to-end run PH0 to PH5 (2026-11-13)

Prototype E6 of ADR-004 section 5: text in, packed text + suffix array + compact LCP + child table (layout A, codec Packed) + statistics out, **all stages in one process**, bounded RAM, every big output streamed
straight into a memory-mapped file, the scratch deleted at the end. Settings fixed by the owner on 2026-11-13 (ADR-004 section 6c): child layout A with the codec Packed, PH4 in base mode, PH5 fused with PH4 through the
B3 hook (`ILcpChunkConsumer`), no run-length encoding of the stack, 4 threads unless a run says otherwise. Gates G1-G4 were fixed by the orchestrator before the work (task statement; quoted in section 5).
Machine A5 (4 vCPU, 15 GB, virtio disk, ext4), cold page cache for every timed run (`sync` + `drop_caches`), fresh process, `DOTNET_GCHeapHardLimit` = 0x6000000 (96 MiB), real block sizes, scale-down emulation by
`--scale-n 3000000000`. Raw rows: `e6_main.tsv`, `e6_main_b1.tsv`, `e6_threads.tsv`, `e6_kvar.tsv`, `e6_dirty.tsv`, `e6_small.tsv`, `e6_small_b1.tsv`, `e6_setup.tsv` (this folder); scripts `run_e6_setup.sh`,
`run_e6_protocol.sh`, `sample_e6.py`, `tables_e6.py`; mutation audit `../mutation/e6_*`, `results_e6.tsv`.

**A restart happened.** The container was restarted once while the work was in progress (after the tests and the driver were written, before any measurement). The uncommitted tree survived; nothing was lost
except a test run and a build that were killed, both repeated. Steps that were repeated: the build and the WholePipeline test run (twice: the first run on a disk was too slow, see section 4). All measurements were
made after the restart, from a cold cache.

## 1. What was built

`apps/SuffixTree.EmProto`:

| file | content |
|---|---|
| `Pipeline/EndToEndOptions.cs` | `EndToEndOptions` (input, output prefix, scratch directory, threads, `ScaleN`, difference-cover period, codec, arena caps, hooks for the tests), `EndToEndTuning` (explicit values that win over the rule), `EndToEndParameters` (every size of one run), `ScaleRules.Resolve` (the scale-down rule of E3 / E3b / E4 applied to every stage, section 3) |
| `Pipeline/Prepare.cs` | PH0: alphabet scan, symbol counts, SHA-256 of the input, the b-bit packed text streamed into its file |
| `Pipeline/EndToEnd.cs` | the driver: PH0, PH1, PH2/PH3, PH4 fused with PH5, per-phase records, scratch bookkeeping, final-length check, clean-up |
| `Experiments/E6.cs` | commands `e6` (run + report + `--verify-*` against reference files) and `e6-mktext` (the text file of a family; setup, outside the pipeline) |

Reused unchanged: `ExternalDcs.BuildRankFile` (E3b), `EmBuilder` + `MappedSaSink` (E3), `ExternalLcp.Build` with its B3 hook (E4), `ChildScan` through `E5.LazyScan` (E5), `MappedWindowWriter`, `PackedText.PackInto`.
**Additive changes to existing code** (all existing tests stay green): `MappedSaSink.Write` adds its bytes to `MappedIo` (the SA writes were not counted before; E3b's mapped-write row now includes the SA writes of its external
levels, which its note said it did not); `DirtySampler` gets `BaselineKb` and `TakePhaseMax()` (per-phase maxima); `E5.LazyScan` is `internal` instead of `private`; the `oracle` command takes `--out` (writes the kernel's suffix
array as a file, the byte-level reference of G1); `Cli.cs` registers `e6` and `e6-mktext`. `src/**` is untouched.

Outputs of one run (prefix `P`): `P.pk` + `P.pk.meta` (packed text, n, bits, sigma, letters, counts, SHA-256), `P.sa` (4 B/char), `P.lcp` + `P.ovf` (compact LCP, ADR-003 form, as E4 writes it), `P.cld` (CHILD, layout A, Packed,
4 B/char), `P.stats` (NODE_COUNT, internal nodes, LRS_DEPTH, LRS_RANK, width, overflow count). The existing commands read these files unchanged (`e4`, `e5-oracle`, `lcp-oracle`, `check-lcp`, `hash-sa`).

## 2. Design choices

* **PH0 reads the input twice** (pass 1: alphabet, counts, SHA-256; pass 2: map to dense symbols in increasing byte order and pack) because the width b is known only after the scan: 2 B/char of reads instead of the 1.25 of the ADR
  table. Rejected: packing with a guessed width (2 bits) and an exception list (Q3 asked for b bits with a 2-bit fast path, no exception codec), or packing at 8 bits and re-packing (a third pass over a 1 B/char file). The pack
  works on blocks of a multiple of 8 symbols (so that every block ends on a byte boundary) with the existing `PackedText.PackInto` into a `MappedWindowWriter` (file created at its final length, windows of 1 MiB); the
  second pass recounts the symbols and fails if the input changed. The input length is taken through a handle: `FileInfo.Length` of a symbolic link is the length of the link (**a defect found by the protocol group `small`**,
  where `dna_10m`, `rep_10m`, `prot_1m` are links to the benchmark inputs; fixed, two tests).
* **One process, sequential stages, one scratch directory** (`<scratch>/ph1`, `ranks.bin`, `ph3`, `ph4`, `ph5`), created empty (a non-empty one is refused) and deleted with everything in it at the end, also when a stage throws.
  The rank file is the only scratch that crosses a stage boundary: it is tracked, and **deleted as soon as PH3 has used it** (its only consumer). The packed text and the SA are final outputs read by the later stages through
  read-only mappings and sequential reads.
* **PH4 base mode fused with PH5** through the existing hook: `E5.LazyScan` creates the scan at the first LCP chunk, so its stack blocks and windows do not exist during B1/B2; `Complete()` checks that every rank arrived.
  The LCP is consumed once, as it is written.
* **Arena budgets**: 160 MiB each for PH1, PH3 and PH4 (they run one after the other), 64 MiB for PH5 (next to PH4). The number of parallel B2 tasks is derived from the cap at the reference geometry (R = 8,000,000), as E4
  does (2 of 4).
* **Accounting**: per phase wall time, the program's own byte counters (explicit reads and writes through `CountedFile`, writes through counted mapped windows), `/proc/self/io` deltas, an RssAnon sampler per phase (10 ms),
  the dirty-page sampler of E5 (system-wide Dirty and the per-process mapped dirty pages of the output files), the tracked scratch bytes. CHILD's windows are not counted by `AlignedWindowWriter`, so its 4 (n + 1) bytes are
  added at the run level (the file is written once).
* **Sync**: the outputs are fsynced at the end as their own phase `sync` (0.006 s; the data was already on the device: see section 6), reported next to, not inside, the stage times.
* The final lengths of all five big outputs are checked at the end of the run (the packed text with its padding, SA and CHILD 4 (n + 1) bytes, base LCP (n + 1) x width, overflow pairs).
* Not built (out of scope): PH6, PH7 (seal, container), resume, positions >= 2^32, an empty text (refused with `ArgumentException`).

## 3. The scale-down rule (every stage, the same rule as E3 / E3b / E4)

With f = n / `ScaleN`: PH1 takes K_s, N_base, R, K_sort as their 3x10^9 reference values (1 Mi, 16 Mi, 4 Mi, 1 Mi) times f, `Math.Round`, floors 3 / 4 / 16 / 3, and the run-spill bytes 64 KiB x f (E3b); PH3 takes K = n x 10^6 / `ScaleN`
(= n / 3000, E3's `--k-div 3000`); PH4 takes R = ceil(8,000,000 f), floor 3, and 2 B2 tasks from the arena cap at the reference geometry (E4); PH5 has no size that depends on n. Values at n = 100,000,000 (printed by the run
and pinned by `ScaleRulesTests`, computed by hand):

| stage | parameter | value | at 3x10^9 |
|---|---|---|---|
| PH1 | K_s / N_base / R / K_sort / run spill | 34,953 / 559,241 / 139,810 / 34,953 / 2,185 B | 1 Mi / 16 Mi / 4 Mi / 1 Mi / 64 KiB |
| PH3 | K (fan-out 64) | 33,333 | 1,000,000 |
| PH4 | R, P = n / R, B2 tasks | 266,667, 375, 2 | 8,000,000, 375, 2 |
| PH5 | block, window x slots | 262,144 entries, 1 MiB x 8 | same |

Resulting shape (dna_100m / genome-like): DCS sample 3,222,557, 185 sample buckets, recursion depth 0 / 1 (0.84 % shared names); 5,381 / 5,374 leaf buckets of at most 33,333 positions, 711 / 713 partitions, depth 3; 375 LCP chunks;
U8 base LCP with 0 / 1,390,028 overflow pairs.

## 4. Tests

New: 79 tests in `tests/SuffixTree/SuffixTree.EmProto.Tests/E6` (namespace `WholePipeline`): `EndToEndTests` (small texts and their enumeration through CHILD and the Ukkonen logical hash; **exhaustive binary 1-12 and ternary 1-6, 9,282
strings, through the whole pipeline with tiny parameters**; 380 tiny texts with the enumeration through the output and the Ukkonen logical hash; the adversarial families at 1,000 / 5,000 / 100,000 symbols; a^n, Fibonacci,
Thue-Morse at 10^6; **thread (1, 2, 3, 4) and K identity** on five 10^5 texts, each also equal to the in-RAM structures; phase order and scratch release stage by stage; failure injection at every phase; refused inputs;
counters; the fused hook; alphabets of 5 bits and of all 256 bytes; the tagged codec; a symbolic link as input), `EndToEndPlumbingTests` (stage options against the resolved parameters, scale-down run, command line, `e6` against the
reference files of the other stages with corrupted references, the RssAnon cap in the exit code, the report rows, the harness additions), `PrepareTests`, `ScaleRulesTests`, `EndToEndFuzzTests` (a seeded 40-case sample in the
default run; the explicit differential fuzz with K = 3..50, F = 2..4, X = 7..64, tiny blocks, windows and chunks, texts of 10^3 to 10^5 symbols, 1-4 threads).
The tests run on `/dev/shm` when it exists (`EMPROTO_TEST_TMP` overrides): the tiny runs msync and fsync tens of small files each, which on the disk took 50-100 ms per run (the first complete run of the exhaustive test had not
finished after 7 minutes); on a memory file system the 79 tests take 7 minutes, of which 6 for the exhaustive enumeration.

**Tests first**: the 75 tests of the first version were run against a stub (`NotImplementedException` in `EndToEnd.Run`, `Prepare.Run`, `ScaleRules.Resolve` and the `e6` commands): **71 red** (every driver, PH0, scale-rule and command-line test;
the reason is the stub in each), 4 green (`TheConstants_AreThoseOfTheStageCommands` and the three tests of the additions to the harness, which were made before the stub). After the implementation the first run of 69 of them (without the exhaustive enumeration, the 10^6 families and the seeded sample) had 68 green; the one failure was an
assertion of mine (`> 12 B/char`, the peak is exactly 12.0 B/char: `>=`), the code was not changed for it; the others followed in the complete run. Existing EmProto tests: 255, unchanged and green.

| run | result |
|---|---|
| default build, whole project (tests commit `b067d9e2`) | **334 of 334** (255 old + 79 new), 10 m 53 s |
| checked build (`-p:EsaBounds=true`, `ESA_EXPECT_BOUNDS=true`), tests commit `59c44319` (before the two symbolic-link tests) | **332 of 332**, 11 m 30 s |
| explicit fuzz of the whole pipeline, seed 1, 120 s / seed 2, 60 s | **507 / 241 cases, all equal** (280 / 137 with a stack spill) |

## 5. Gates (fixed by the orchestrator before the work) and results

Inputs: `dna_100m` = 100,000,000 uniform random DNA (generator `dna:100000000`, seed 20240001; the text of E3, E4, E5) and the genome-like 100 M (`genome:100000000`, seed 20240002; ACGNT, repeats, satellites, six N runs, 12.1 % of the
positions in key ties). The input of the pipeline is the ASCII text file (`e6-mktext`, 100 MB, one byte per symbol). References, made by `run_e6_setup.sh` with the in-RAM code of the earlier stages: the packed file of `pack`,
the kernel suffix array `oracle --out` (SaisKernel in RAM; its SHA-256 equals the E3 setup rows `ed12eb83...` and `d1bfd73a...`), `lcp-oracle` (LcpKernel on that SA), `e5-oracle` (the in-RAM `ChildTable.Build`, packed layout, and the
`EsaStatistics`). For `dna_100m` the SA, the LCP base and overflow files and the packed child table are byte-identical to the orchestrator's files in `/tmp/stbench/orch4` and `/tmp/stbench/orch5` (cmp; the packed text differs
only in the 86 padding bytes that the E3 `pack` of 2026-11-09 wrote, the genome-like one is identical).

### G1 identity (hard)

| clause | result |
|---|---|
| 100 Mbp, 4 threads, dna_100m and genome-like: pk, SA, LCP base, overflow, CHILD byte-identical to the references (`cmp`), NODE_COUNT / internal / LRS_DEPTH / LRS_RANK equal | **met**: 6 of 6 runs per input (final binary), every comparison `yes` (32 of 32 rows per input; run 1 of each also compared by the program itself with `--verify-*`) |
| the same for 1 and 2 threads | **met**: 1 run per input and thread count, 9 of 9 comparisons each (identical to the same references, hence to the 4-thread output) |
| two values of K (K / 2 = 16,666 for PH3 with K_s = 17,476 for PH1, and 2 K = 66,666 with K_s = 69,905), 4 threads | **met**: 9 of 9 comparisons each, both inputs |
| the same with the scratch on tmpfs (3 + 3 runs) | met: 27 of 27 comparisons each |
| dna_10m, rep_10m, prot_1m: SA against the committed libsais hash (`analysis/sa_vectors.tsv`), everything else against the in-RAM references | **met**: `sa_equals_libsais` yes, 10 of 10 comparisons each (prot_1m has a 5-bit alphabet: supported) |
| a^n, Fibonacci, Thue-Morse at 10^5 and 10^6: SA against libsais, the rest against the in-RAM references; tandem (171 x 600, 171 x 6000), N run (10^5 and 10^6), genome-like 10^6: against the in-RAM references (no libsais hash exists for them) | **met**: 10 of 10 (a, fib, tm) and 9 of 9 (the others) |
| all adversarial families of `StringFamilies` (that fit in bytes) at 10^3 / 5 x 10^3 / 10^5 and a^n, Fibonacci, Thue-Morse at 10^6 against the in-RAM structures, tiny and moderate parameters | **met** (unit tests: more than 40 cases at 1,000 and 5,000, more than 25 at 10^5, 3 at 10^6) |
| exhaustive binary 1-12 and ternary 1-6 through the whole pipeline, tiny parameters | **met**: 9,282 of 9,282 equal (SA, LCP files, CHILD, statistics, metadata, packed text, clean-up) |
| differential fuzz of the whole pipeline, K = 3..50, F = 2..4, X = 7..64, tiny blocks / windows / chunks, texts 10^3-10^5, 1-4 threads | **met**: seed 1, 121 s: 507 cases (280 with a stack spill), seed 2, 60 s: 241 cases (137 with a spill), every case crossed a partition, a recursion level or a chunk boundary, all equal; plus the seeded 40-case sample of the default run |
| enumeration through the output equals the Ukkonen logical hash on tiny texts | **met**: 380 texts (binary to length 7, ternary to 4, six named), enumeration through the packed CHILD file equal to `EsaIndex.Traverse` and to the Ukkonen logical hash |

### G2 RssAnon <= 256 MB incl. the .NET runtime (peak over all stages), real block sizes, 96 MiB GC limit

| input | program sampler (10 ms), median [min-max] of 3 | independent sampler (20 ms, `sample_e6.py`, `/proc/<pid>/status`) | result |
|---|---|---|---|
| dna_100m, 4 threads | **89.5 MB** [87.2-94.8] | 89.5 [87.2-94.8] | **met** |
| genome-like 100 M, 4 threads | **95.4 MB** [87.0-97.4] | 95.4 [87.0-97.4] | **met** |

The two samplers read the same peak in every run (a plateau of seconds). Per phase (median, dna_100m / genome): PH0 13.8 / 13.9 MB, PH1 29.1 / 34.3, PH2/PH3 61.0 / 81.0, **PH4 + PH5 89.5 / 95.4 (the peak)**. PH4 + PH5 alone
measured 50.8 and 51.0 MB in E4 / E5, so about 40 MB of the fused phase of the pipeline is what the earlier stages leave behind in the process (freed native memory that the allocator keeps, the managed heap is 9.6 MB at the end);
this was not investigated further. Other runs: 1 thread 52.9 / 53.1 MB, 2 threads 67.7 / 69.7, K / 2 121.6 / 105.7 (the largest of the matrix), 2 K 92.8 / 104.9, scratch on tmpfs 107.7 [105.9-109.1] / 104.9 [99.8-115.9]. The peak of the
arenas charged to one stage budget is 30.9 MB (caps 160 / 64 MiB), live bytes at the end 0 in every run. Margin to the cap: 160 MB.

### G3 (row E6) at 4 threads, cold cache, median of 3 runs

| clause | dna_100m | genome-like 100 M | result |
|---|---|---|---|
| total wall (PH0-PH5, without the final sync) <= 90 s | **24.7 s** [24.0-25.6] | **25.0 s** [24.6-25.2] | **met** (3.6x margin) |
| SA stage (PH0-PH3) <= 30 s | **16.9 s** [16.5-17.5] | **17.2 s** [17.1-17.2] | **met** |
| total I/O, program counters <= 100 B/char | **71.6 B/char** | **73.1 B/char** | **met** |
| cross-check with `/proc/<pid>/io` rchar + wchar (syscalls only, the mapped writes are not in it) | 65.8 B/char (independent sampler, whole process, runs 2-3: 65.8) | 67.0 B/char | program's explicit counters 62.2 / 63.3 B/char, 5.5 % below the syscall total (the difference is not attributed: reads that are not counted by `CountedFile`); the mapped writes 9.4 / 9.7 B/char come on top of both |
| mapped dirty pages of the outputs <= 64 MB (per-process view of the output mappings, Private + Shared Dirty in smaps every 20 ms) | **9.4 MB** | **10.1 MB** | **met** |
| dirty pages of the outputs in the page cache, system-wide Dirty minus baseline, scratch on tmpfs (3 runs) | 11.3 MB [11.1-11.4] (independent 11.3) | 12.7 MB [12.6-12.9] | **met** (isolation run: no scratch on the disk) |
| system-wide Dirty with the scratch on the disk (reported separately, not part of the clause) | 1,204 MB (PH4's scratch, as E4 / E5 found: 1.2 GB; 396 MB in PH2/PH3 from the position files) | 1,201 MB | - |

Total I/O by phase and the explicit / mapped split are in section 6 (G4); the writes through the mapped windows (SA, LCP, CHILD, text, rank file: 9.4 B/char) stay below 11 MB of dirty pages at a time because the windows are bounded
(the samplers see no accumulation of the output pages; the mechanism, that the unmap path of the runtime flushes a writable view, was not verified).

### G4 (report only)

Per phase, median [min-max] of 3 runs, 4 threads, cold cache (`tables_e6.py`; B/char = per input symbol; "explicit" = program counters of `CountedFile`, "mapped" = writes through mapped windows, CHILD's 4 B/char added at the run level
and not in the per-phase column; rchar + wchar = `/proc/self/io` deltas of the phase):

**dna_100m**

| phase | wall s | explicit read | explicit write | mapped write | rchar + wchar | RssAnon peak MB | scratch peak so far B/char | dirty peak MB (system / mapped outputs) |
|---|---|---|---|---|---|---|---|---|
| PH0 prepare | 0.535 [0.526-0.536] | 2.00 | 0 | 0.25 | 2.07 | 13.8 | 0 | 1.1 / 1.05 |
| PH1 DCS | 1.71 [1.70-2.04] | 0.258 | 0.258 | 0.129 | 0.755 | 29.1 | 0.258 | 1.4 / 0 |
| PH2/PH3 SA | 14.7 [14.3-15.0] | 7.85 | 7.84 | 4.00 | 17.8 | 61.0 | 4.30 | 396 / 0.21 |
| PH4 + PH5 | 7.79 [7.49-8.07] | 24.0 | 20.0 | 1.00 (+4 CHILD) | 45.1 | 89.5 | 12.0 | 1,200 / 9.4 |
| sync | 0.006 | 0 | 0 | 0 | 0 | 65.0 | 12.0 | 1.7 / 0 |

**genome-like 100 M**

| phase | wall s | explicit read | explicit write | mapped write | rchar + wchar | RssAnon peak MB | scratch peak so far B/char | dirty peak MB (system / mapped outputs) |
|---|---|---|---|---|---|---|---|---|
| PH0 prepare | 0.532 [0.502-0.559] | 2.00 | 0 | 0.375 | 2.07 | 13.9 | 0 | 1.0 / 1.05 |
| PH1 DCS | 2.19 [2.15-2.44] | 0.828 | 0.828 | 0.262 | 1.96 | 34.3 | 0.391 | 24.7 / 0 |
| PH2/PH3 SA | 14.4 [14.3-14.5] | 7.83 | 7.82 | 4.00 | 17.8 | 81.0 | 4.25 | 396 / 0.26 |
| PH4 + PH5 | 7.80 [7.49-7.96] | 24.0 | 20.0 | 1.11 (+4 CHILD) | 45.1 | 95.4 | 12.0 | 1,200 / 10.1 |
| sync | 0.006 | 0 | 0 | 0 | 0.0007 | 80.2 | 12.0 | 0 / 0 |

Shares of the total (dna_100m): PH0 2.2 %, PH1 6.9 %, PH2/PH3 59.4 %, PH4/PH5 31.5 %. Inside the fused phase: PH4 B1 1.15 s, B2 3.41 s, B3 3.21 s, of which the PH5 scan 2.36 s (genome-like: 1.23 / 3.23 / 3.32, scan 2.47).
Per symbol: 247 ns (dna_100m) and 250 ns (genome-like), 4.05 and 4.00 Msym/s. Time against the earlier prototypes: SA stage 16.9 s = PH0 0.5 + PH1 1.7 (E3b: 16-21 ns/symbol at 100 Mbp) + PH3 14.7 (E3: 14.6 s alone); PH4 + PH5
7.8 s (E5: PH4 5.2 s + 2.4 s for the scan fused: 7.6 s).

| measure | dna_100m | genome-like 100 M |
|---|---|---|
| peak scratch, tracked by the program (apparent) | **12.0 B/char** (1.2 GB, the end of PH4-B1: all Phi tuples on the disk) | 12.0 B/char |
| peak scratch, independent sampler of the scratch directory: apparent / allocated | 11.8 / 11.8 B/char | 11.8 / 11.8 B/char |
| peak scratch of PH2/PH3 and PH1 (so far, before PH4) | 4.30 B/char (PH3 position files), 0.26 (rank file) | 4.25, 0.39 |
| disk footprint of the outputs | **9.25 B/char** (925 MB: text 0.25, SA 4.0, LCP 1.0, CHILD 4.0) | 9.49 B/char (text 0.375, LCP 1.11 with 1.39 M overflow pairs) |
| LCP width, overflow pairs | U8, 0 | U8, 1,390,028 |
| NODE_COUNT, LRS (depth, rank) | 162,149,501; 25, 14,863,022 | 166,050,648; 499,999, 75,208,039 |
| block layer: bytes read / pages dirtied over the run | 1.02 B/char / 40.7 B/char | 1.02 / 41.6 |

Scaling (4 / 2 / 1 threads, total s): dna_100m 24.7 / 28.6 / 39.3 (1.59x on 4 threads), genome-like 25.0 / 29.9 / 45.9 (1.83x); SA stage 16.9 / 20.9 / 29.8 and 17.2 / 22.6 / 34.9. K / 2 (more buckets): 27.1 / 27.6 s, 2 K: 20.7 / 21.2 s.
The scratch on tmpfs instead of the disk: 14.7 [14.6-16.3] s and 17.6 [17.0-17.7] s, i.e. **8-10 s of the 25 s are the disk I/O of the scratch** (virtio disk, probably host-cached).

**ESTIMATE (linear extrapolation, labelled; nothing above 100 Mbp was run in E6)**: 3x10^9 symbols = 30 x the 100 Mbp numbers: total about 740-750 s = 12.4 min (0.21 h), SA stage about 510 s, PH4/PH5 about 235 s; scratch 36 GB
(12 B/char; the base mode of PH4), outputs about 28 GB (9.25-9.5 B/char), I/O 72-73 B/char. The per-symbol time of PH4 at real size is up to 1.5x the scale-down one (E4: 13.0 against 19.2 Msym/s), PH1 at 3x10^9 was measured at
36.7 s (E3b) against 30 x 1.7 = 51 s linear, so the range is about 12-15 min; the abandon threshold of the row (4.5 min at 100 Mbp, i.e. 2.3 h at 3 Gbp) is far away. E7 measures the exponent.

## 6. Not verified, caveats

* Anything above 100 Mbp, 40-bit positions, a text of 2^32 - 2 symbols and more, an empty text (refused), a text with more than 256 distinct byte values (the input is a byte file).
* A chain at E6 level: a^n and the N run were run through the whole pipeline only up to 10^6 symbols (with tiny stack blocks, so that the spill path ran); the 50 M and 200 M runs of E5 were PH5 alone, and the E6 protocol did not repeat them.
* Page-cache-limited runs of the whole pipeline (E7), a real disk (the virtio disk is probably host-cached: the 8-10 s of scratch I/O and the 0.006 s of `fsync` are its), other kernels.
* The 100 Mbp runs of the thread, K and tmpfs groups (`e6_threads.tsv`, `e6_kvar.tsv`, `e6_dirty.tsv`, `e6_main_b1.tsv`) were made with the binary of commit `59c44319`; the group `main` and the small inputs were repeated with the
  final binary (commit `b067d9e2`), which differs only in how the input length is read (a symbolic link as input). `e6_small_b1.tsv` holds the first small-input run, in which three cases (`dna_10m`, `rep_10m`, `prot_1m`: links)
  ran on the 29 bytes of the link and failed every comparison: that is how the defect was found; the rows are kept.
* The independent sampler's whole-process `rchar` of run 1 includes the verification reads of the program itself (`--verify-*`); runs 2 and 3 do not verify, the medians use them.
* Reads through the mappings (text, rank file in PH3) are not in the I/O counters; the block layer read only 1.02 B/char (cold text and rank pages; the rest was written by the run itself and stays in the cache).
* The oracles are the in-repo in-RAM structures and kernels (SaisKernel, LcpKernel, `ChildTable`, `EsaStatistics`), independent of the external-memory code, but not an independent implementation; the libsais hashes cover SA only
  for the inputs that have one, the Ukkonen logical hash only for tiny texts.
* The fused phase's RssAnon is about 40 MB above PH4 + PH5 alone (section 5, G2): not investigated.

## 7. Mutation audit of the new driver code

Hand-made single-site mutants of the new code only (the stage code was audited with E3b, E4 and E5): `Pipeline/EndToEnd.cs`, `Pipeline/Prepare.cs`, the scale rules in `Pipeline/EndToEndOptions.cs`, `Experiments/E6.cs`. Written by
`make_e6_mutants.py` into `analysis/mutation/e6_<file>_listed.json` (**99 candidates**, 38 + 17 + 20 + 24) and the selection that was run into `e6_<file>_mutants.json`; run by `run_e6_mutants.py` (every build and every test run with a
watchdog that kills the whole process group, a hang counts as killed; stage 1 = the WholePipeline tests without the exhaustive enumeration, the enumeration through CHILD, the seeded fuzz sample and the thread / K identity,
about 20 s; a survivor would go to stage 2 = all but the exhaustive enumeration) in two separate git worktrees (`/tmp/stbench/wt_e6a`, `wt_e6b`, removed afterwards), two at a time on the 4 vCPUs. Only 34 of the 99 candidates were
run: the build of the test project costs 1-2 minutes per mutant and the economy limit of the task did not allow more; the selection covers every category of the task statement.

| file | run | compiled | killed | not compiled (replaced) | survivors |
|---|---|---|---|---|---|
| `EndToEnd.cs` | 14 + 2 | 12 + 2 | 14 | P13, P15 (S2737 catch that only rethrows, S1481 unused local) -> P13b, P15b | 0 |
| `Prepare.cs` | 5 + 2 | 3 + 2 | 5 | R01, R05 (S108 empty block, S1125 redundant boolean) -> R01b, R05b | 0 |
| scale rules (`EndToEndOptions.cs`) | 5 | 5 | 5 | - | 0 |
| `E6.cs` | 6 | 6 | 6 | - | 0 |
| **total** | **34** | **30** | **30** | 4 | **0** (0 equivalent) |

The mutants run: PH1 and PH3 and PH4 given another stage's parameter (P02, P06; S01 truncation instead of rounding, S03 rounding of R down, S04 the K reference 2^20, S08 an explicit R ignored, S11 the unscaled run on the scale rule),
the thread count ignored by PH3 and PH4 (P07, P09), PH4 in the bounded mode (P28), the statistics written after the sync (P01, phase order), the scratch directory not deleted after a run and after a failed run (P12, P13b), the rank file not
deleted after PH3 and the tracker not told (P14, P15b), the SA file of n entries (P17), PH5 not attached to the hook (P23), **the hook consumer not called for the last chunk (P24) and the chunk boundary handed to PH5 off by one (P25)**,
PH0 without counts / with a block that is not a multiple of 8 / without the comparison of the passes / with the last partial byte truncated / without padding (R01b, R05b, R07, R11, R12), the command line (default threads, sync flag,
the packed text and LRS_RANK and the exit code not compared, the text file holding the dense codes: C01, C02, C09, C12, C14, C18). The tests that kill them are by design the ones that look at the glue: the stage options against the resolved parameters, the phase order and scratch release by stage, the failure injection, the final-length and counter checks, the PH0 unit tests, the command-line and `e6` verification tests (the result file has the number of failing tests per mutant, not their names). One test was added during the work because a mutant class would otherwise have survived by design (the scratch tracker reset: `TheScratchPeak_IsThatOfTheRunItself_...`; its mutant
P28 is killed by 10 tests).

## 8. Commits

`37caa758` tests first (71 of 75 red), `4b2e448b` driver, `59c44319` tests on tmpfs / report tests / scripts / mutant lists, `b067d9e2` input-length fix, `fe99854e` raw rows, then the audit, results and documentation commits.
