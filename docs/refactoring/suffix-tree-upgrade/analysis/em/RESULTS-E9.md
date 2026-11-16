# RESULTS-E9: PH6, Psi by sigma-way distribution of the BWT, against ISA by DIST, and the suffix-link smoke test (2026-11-16)

Prototype E9 of ADR-004 section 5 (the owner authorised it on 2026-11-16, section 6c: "E4-E9, then 1b.7"). Gates G1-G4 were fixed by the orchestrator before the work (task statement; quoted in the table). Machine A5 (4 vCPU,
15 GB RAM, virtio disk, ext4, ~13 GB free), cold page cache for every timed run (`sync` + `drop_caches`), fresh process, `DOTNET_GCHeapHardLimit` = 0x6000000 (96 MiB) for every run of the protocol, 3 runs per cell, medians. Binary:
`dotnet publish apps/SuffixTree.EmProto -c Release -o /tmp/stbench/bin_e9_run`: the timed rows (G2, G4: `e9_protocol.tsv`) are from the binary of commit `8d79fc47`, the identity and smoke rows (G1, G3: `e9_identity.tsv`) from the binary of commit `a7c3b52a`; `apps` differs between the two only in `Ph6/LinkSmoke.cs` (two more counters in the smoke cells, the mismatch expression decomposed, no change to what is compared), which no timed command runs. Raw rows (this folder):
`e9_setup.tsv`, `e9_identity.tsv` (G1, G3), `e9_protocol.tsv` (G2, G4: 93 timed runs, all exit code 0, 90 of 90 output comparisons equal, scratch empty after every run), `e9_bwt_check.txt` (independent BWT check). Scripts: `run_e9_setup.sh`, `run_e9_identity.sh`, `run_e9_protocol.sh`,
`sample_e9.py` (independent sampler), `check_bwt_e9.py` (independent check of the BWT), `tables_e9.py` (the tables below are its output); mutation audit: `../mutation/*e9*`.

**No restart of the container happened.** One defect of the product code (an exception of the smoke driver on a text without internal nodes) was found by a new test and fixed; one measuring-script defect (a heap limit on the smoke test) was found
and fixed; one finding about the cost of small mapped windows changes the recommended window size (section 5).

## 1. Result in one table

| gate | threshold | measured | result |
|---|---|---|---|
| G1 Psi file byte-identical to Psi of the in-RAM oracle (P1 `LinkTables`), `PSI[0] = 0` | identical on dna_10m, prot_1m, rep_10m, a^n, Fibonacci, Thue-Morse (10^6); property tests, sigma = 1, 2, 4, 5, 20, 64, 65, 255, 256, window sizes down to a few bytes, the row of the suffix 0 anywhere | **identical in all 8 inputs x 5 comparisons (40 of 40 `cmp`)**: the 6 required inputs plus a random 256-symbol text (10^7 bytes, 4 passes at the default) and dna_100m; each input with the default parameters, with 1 KB windows and 2 symbols per pass (several passes on any alphabet), cyclic variant (tail equal, entry 0 = the row of the suffix 0), ISA by DIST equal to the in-RAM ISA, Psi from the DIST ISA equal; 51 tests in `EmProto.Tests.Ph6`, among them property tests over all 9 alphabets (windows 1-65536 bytes, 1-256 symbols per pass, both conventions), exhaustive texts (binary to length 9, ternary to 6, 4-ary to 4), the adversarial families at 9 lengths, a^n / Fibonacci / Thue-Morse at 30,000, the sentinel row at row 1, an inner row and row n | **met** |
| G1 ISA from DIST byte-identical to the in-RAM ISA | same inputs | **8 of 8 equal** (chunk counts P = 375 at the scale-down of 3x10^9; tests with R = 1 .. 2^30 and buffers of one tuple) | **met** |
| G2 I/O of the Psi pass (BWT 1 B/char read + Psi 4 B/char written) | <= 6 B/char, three counters | **program counters 5.000 B/char** (dna_100m, dna_10m, prot_1m; 1.000 read + 4.000 written through mapped windows); **rchar + wchar 1.000 B/char** of the pass alone (`--no-sampler`; the mapped writes are not syscalls: + 4.000 = 5.000); **block layer (net of the runtime start-up): 5.19 B/char with 1 MiB windows** (dna_100m: 1.02 read + 4.17 written), **7.16 B/char with the 64 KB windows of the task statement** (1.04 + 6.12) | **met with 1 MiB windows (all three counters); the block-layer counter is NOT met with 64 KB windows (7.16, dna_100m; 6.97 dna_10m)**; program counters and rchar + wchar are met at every window size. Cause: section 5 |
| G2 ISA by DIST next to it | about 24 B/char | **23.73 B/char counted** at dna_100m (reads 11.86 + writes 7.86 + mapped 4.00; P = 375, tuples through 375 chunk files, 7.86 B/char of scratch), 23.67 at dna_10m; + Psi from the ISA 8.00 (SA read 4.00 + Psi 4.00; the ISA is read at random through its mapping, not counted): 31.7 B/char against 5.0 | reported |
| G2 where the BWT comes from | stated | route (ii), a separate pass over the SA file and the packed text: **5.000 B/char counted** (SA read 4.000, BWT written 1.000; the text 0.25 B/char is read at random through its mapping, in the block layer's 4.30 read), 2.43 s [2.36-2.97] at dna_100m. **Psi pass + BWT pass = 10.0 B/char** (10.5 B/char at the block layer with 1 MiB windows: 5.30 + 5.19). The "for free in PH3" route (i) is **not built, not measured** | reported |
| G3 smoke test: every kernel x source of P1 with the Psi of the stream (and of the ISA route) and the ISA of DIST gives the same link, success flag and step count as with the in-RAM tables, on the P1 samples | 200,000 internal nodes (145,561 prot_1m), seed 1; equal to the definition on the small adversarial texts | **0 mismatches in 7 inputs x 2 Psi routes x 8 cells (112 cells)**: dna_10m (6,229,745 internal nodes, sample 200,000), rep_10m (8,568,247; 200,000), prot_1m (295,123; **145,561**), rnd256_10m (2,093,383; 200,000), and a^n / Fibonacci / Thue-Morse 10^6 (samples of 5,000 with a 20,000-step budget; on a^n 4,968 of 5,000 descent links exceed the budget, identically in both runs); the definition on 1,725 exhaustive texts and the adversarial families at 9 lengths (`LinkSmokeTests`, 3 Psi variants x 2 fan-outs each) | **met** |
| G4 RssAnon <= 256 MB incl. runtime; arena independent of n | <= 256 MB; sigma x window <= 160 MiB | independent 20 ms sampler, whole process: **10.9 MB** (64 KB windows) / **10.0 MB** (1 MiB) at dna_100m with the in-process samplers, **4.4 MB** (1 MiB) / 5.3 MB (64 KB) without them; ISA by DIST 41.3 MB; BWT pass 10.2 MB; arena peak of the pass **1.31 MB (64 KB windows) / 5.24 MB (1 MiB) at both dna_10m and dna_100m** (the 1 MiB read buffer plus the 4 windows; prot_1m, 20 symbols: 2.31 MB = its n + 1 read bytes plus 20 windows), bound sigma x window = 64 x 64 KB = 4 MiB (16 MiB for 256 symbols) | **met** |
| G4 time of the Psi pass at dna_100m, cold cache | reported; ESTIMATE for 3x10^9; abandon if projected > 5 min | **2.66 s [2.61-2.68] with 64 KB windows, 1.10 s [1.08-1.13] with 1 MiB** (process wall 2.84 / 1.30 s; 10.4 ns per entry with 1 MiB, 2.26 s [2.19-2.28] / 1.04 s without in-process samplers). **ESTIMATE for 3x10^9 (linear, x30, not measured): 80 s (64 KB) / 33 s (1 MiB); BWT pass 73 s; ISA by DIST 58 s + Psi gather 134 s** | **met** (the abandon rule does not fire: 80 s against 300 s) |

Reading. (1) **The streaming Psi is exact**: 40 of 40 file comparisons, and the P1 kernels cannot tell its tables from the in-RAM ones on up to 1.6 million links per input (all kernels, both sources) and on the definition on 1,725 exhaustive small texts and the adversarial families. (2) **The 6 B/char
gate holds on the program counters (5.000) and on rchar + wchar (1.000 + the 4.000 of mapped writes), and on the block layer only with large windows**: with 64 KB windows the device saw 6.12 B/char of writes for 4.00 of data (and the pass
took 2.4x as long), because the .NET `MemoryMappedViewAccessor` flushes (`msync`) every window when it is disposed (strace: 613 `msync` calls for 612 windows) each flush being a synchronous write (the extra device bytes are attributed to journal commits and to boundary pages that are written twice; the split was not measured); 1 MiB windows bring it to 4.17 B/char (4 MiB: 4.04), 16 KB windows to 12.9. No redesign: the window size is a parameter, and 1 MiB is the value of the ADR's PH6 row ("sigma windows x 1 MB"). (3) The ISA route costs 31.7 B/char counted against
5.0 for the same table (6.3x), and 1.9 s + 4.5 s = 6.4 s against 1.1 s (2.7 s with 64 KB windows) at dna_100m; the ISA is read from the device at random through its mapping (4.0 B/char of block-layer reads, cold cache), which at 3x10^9 means random reads over a 12 GB file that is larger than the page cache. (4) The BWT is the larger half of the bill: 5.0 B/char as a separate pass, 1.0 if PH3 wrote it,
which is the one thing this prototype does not show.

## 1a. Gates, one by one (as the task statement states them)

* **G1 identity: met.** Details in section 4. Sentinel row: the BWT symbol of the row of the suffix 0 is stored out of band (a dummy 0 byte in the file, the row in `<bwt>.meta`), because a 256-symbol alphabet leaves no byte for it;
  `PSI[0]` = 0 in the file (`PsiOptions.CyclicEntryZero` = false); P1's `psi[0] = -1`: `PsiFiles.LoadForP1` maps entry 0 to -1 on load and `WriteOracle` writes P1's -1 as 0; the kernels never read entry 0 (a test sets it to 0, to the row of the suffix 0, to -1
  and to 12345: no link changes). The cyclic variant (entry 0 = the row of the suffix 0) is the true successor permutation, one cycle of n + 1 rows (tested; at dna_100m the entry 0 of the cyclic file is 87,373,614, the primary row of the BWT meta).
* **G2 I/O: met on the program counters and on rchar + wchar, met on the block layer with 1 MiB windows and not met with 64 KB windows** (section 3).
* **G3 smoke test: met** (section 4): 7 inputs at scale (dna_100m has no smoke: an in-RAM index of 10^8 symbols is not part of the gate).
* **G4 memory and time: met** (section 5), the 3x10^9 numbers are ESTIMATES.

## 2. What was built (code: `apps/SuffixTree.EmProto/Ph6`, commands in `Experiments/E9.cs`)

Nothing that exists was changed: new files only, plus the dispatch and usage lines of `Cli.cs` and, in `SuffixTree.EmProto.csproj`, six `Compile Include` lines that link the P1 sources (`LinkTables`, `LinkContext`, `LinkKernels`, `LinkNodes`,
`LcpRangeHelper`, `StepCounters` of `apps/SuffixTree.Esa.Bench/Links`, unmodified; `MemDriver.cs`, which needs the Ukkonen tree, is not linked) and, in the test project, two `Compile Include` lines for the definitional suffix array and
the definitional lcp-interval tree of the ESA tests. `apps/SuffixTree.Esa.Bench` and `src/**` are untouched.

* **`PsiStream`** (PH6): one scan of the BWT in row order with one output cursor per symbol. `C[c]` = 1 + the number of text symbols below c (the empty suffix, row 0, comes first); the cursor of c owns the Psi entries
  `[C[c], C[c] + count(c))`, writes them strictly sequentially, and holds one **mapped window** of its own (`--window-kb`, default 64 KB; one view of the Psi file per cursor, disposed and re-created as the cursor advances; no explicit
  flush, as every other writer of the prototype). At most `--symbols-per-pass` cursors (default **64**, the ADR's number) are open in one pass; an alphabet with more symbols takes ceil(sigma / 64) passes, each reading the whole BWT.
  The row of the suffix 0 has no symbol (its BWT symbol is the sentinel) and is skipped; `PSI[0]` is 0 (the convention of ADR-003 section 4) or, with `--cyclic`, the row of the suffix 0 (the true cyclic Psi: one cycle of length n + 1).
  The pass refuses a BWT that does not match its counts (a cursor that overflows or ends short, a symbol outside the alphabet, a file of the wrong length) and charges (cursors x window) to an `ArenaBudget` cap (160 MiB default) before it opens
  anything, so the working set is bounded by sigma and the window, not by n.
* **`BwtPass`** (route (ii), below): `BWT[j] = T[SA[j] - 1]` from the SA file (read once, sequentially) and the packed text (random reads through its mapping); one byte per row, n + 1 rows, the row of the suffix 0 holds a dummy 0
  and its number is recorded, with n and the symbol counts, in `<bwt>.meta`. The sentinel cannot be a byte value for a 256-symbol alphabet, which is why it is recorded out of band (the standard "primary index" of a BWT).
* **`IsaByDist`** (the comparison route): ISA by one-level radix distribution exactly as ADR section 3.3 describes it: pass 1 reads the SA and writes the tuples `(SA[i], i)` (8 bytes) to P = (n + 1) / R chunk files by `SA[i] / R`
  through P write buffers (`--io-kb`; a chunk that never outgrows its buffer stays in memory), pass 2 scatters each chunk into an array of R entries and appends it to the ISA file through a mapped window. **`PsiFromIsa`** then reads the
  SA sequentially and gathers `ISA[SA[i] + 1]` at random through the mapping of the ISA file (PSI[0] = 0).
* **`LinkSmoke`**: every kernel (a scan, b descent, c RMQ + PSV/NSV, c2 PSV/NSV) with every endpoint source (ISA, Psi) of P1 is run on the same sample of internal nodes (P1's own sampling: sample A of 2 x 200,000 + 4000 draws, seed 1) with
  P1's in-RAM tables and with the tables under test; the success flag, the target and the step count of every link must be identical (a kernel that throws on wrong tables counts as a mismatch), and, in the tests, equal to the definition.
* **How PSI[0] meets P1.** The file stores 0 (or the cyclic value); P1 keeps -1 in RAM (`LinkTables.BuildPsi`). `PsiFiles.LoadForP1` maps entry 0 to -1 when it loads a file. No link kernel ever reads Psi[0]: the endpoints read are
  `Psi[lb]` and `Psi[rb]` of an internal node, whose lcp-interval value is at least 1, so `lb >= 1` (the empty suffix, row 0, shares no symbol with any suffix). A test sets entry 0 to 0, to the cyclic value, to -1 and to 12345 and the smoke
  test is unchanged.

**Route decision for the BWT: (ii), a separate pass over the SA file and the packed text; not (i), an opt-in BWT output of the SA sink.** Reason: (ii) touches no existing code, so the e6 pipeline cannot change (no flag, nothing to prove "flag off
unchanged"); the cost is a second read of the SA (4 B/char) and the random text reads that (i) would also pay, and the price of choosing it is that the ADR's "emitted for free when the SA slice is produced" is **not measured here**: what is
measured is the separate pass (section 5). What (i) would add to PH3 is, per suffix, one random read of the packed text (the same read this pass makes) and one byte written at the final offset of the BWT file; the leaf sorter holds the sorted
positions at that moment (`LeafSorter.Sort`, before `sink.Write`), so the change is local, but it would put 1 B/char of writes and a text read into the hottest loop of the pipeline, and the measurement of that belongs with 1b.7d, where PH3 is rebuilt.

## 3. G2: I/O volume, three counters

Per-character values are bytes per input character (n = 10^8, 10^7, 10^6). "Counted" = the program's own counters: explicit reads + explicit writes + writes through mapped windows (reads through the mapping of the text or of the ISA are not
in it). "rchar + wchar" = `/proc/<pid>/io` of the pass: the column **phase** is the pass inside the process (it includes the in-process RssAnon and smaps readers, which read `/proc/self/*`: +0.3 B/char at 10^8), the column
**whole process** is the independent sampler's last reading (it includes the .NET start-up: 273 KB); the rows `psiclean*` run with `--no-sampler` and give the pass alone. "Block layer" = `/sys/block/vda/stat` sectors x 512 between the start of the sampler and the exit of the process
(reads) and after `sync` (writes), **minus the same numbers of the `e9.base` runs** (the 27.0 MB of reads of the runtime image and 0.1 MB of writes of a 1,001-row run). Medians of 3 cold-cache runs; the full rows are in `e9_protocol.tsv`.

**The Psi pass, windows of 64 KB (the task statement's value, the default of `e9-psi`) and 1 MiB (the value of the ADR's PH6 row):**

| input | window | counted B/char | rchar + wchar, pass alone | + mapped writes | block layer: read + written = total B/char | phase s | process wall s |
|---|---|---|---|---|---|---|---|
| dna_100m | 64 KB | 5.000 | 1.000 | 5.000 | 1.04 + 6.12 = **7.16** | 2.26 [2.19-2.28] (2.66 with the samplers) | 2.33 [2.28-2.34] |
| dna_100m | 1 MiB | 5.000 | 1.000 | 5.000 | 1.02 + 4.17 = **5.19** | 1.04 [0.99-1.04] (1.10 with the samplers) | 1.10 [1.07-1.10] |
| dna_10m | 64 KB | 5.000 | 1.000 | 5.000 | 1.16 + 5.81 = **6.97** | 0.295 [0.238-0.299] | 0.364 |
| dna_10m | 1 MiB | 5.000 | 1.000 | 5.000 | 1.07 + 4.12 = **5.19** | 0.127 [0.126-0.133] | 0.198 |
| prot_1m | 64 KB | 5.000 | 1.000 | 5.000 | 1.92 + 4.79 = 6.71 (see note) | 0.050 [0.047-0.052] | 0.119 |
| prot_1m | 1 MiB | 5.000 | 1.000 | 5.000 | 1.92 + 4.09 = 6.01 (see note) | 0.036 [0.032-0.040] | 0.117 |

Note on prot_1m: the block-layer reads are 1.92 B/char for a BWT of 1.00 B/char: about 0.9 MB of other pages (code paths of the pass that the 1,001-row baseline did not run) is a large share of a 1 MB input; at 10^7 the excess is 0.1-0.2 B/char and at 10^8 0.03. The
program counters and rchar + wchar are exact at every size (the pass reads exactly n + 1 bytes: `psi_bwt_bytes_read` = 100,000,001 at dna_100m).

**The window size at dna_100m (3 runs each, cold cache; samplers in the process):**

| window | phase s | block layer written, B/char | RssAnon, independent MB | windows opened |
|---|---|---|---|---|
| 16 KB | 6.15 [6.08-6.31] | 12.90 | 13.7 | 24,415 |
| 64 KB | 2.69 [2.36-2.72] | 6.12 | 10.9 | 6,105 |
| 256 KB | 1.41 [1.29-1.44] | 4.65 | 10.2 | 1,528 |
| 1 MiB | 1.10 [1.08-1.13] | 4.17 | 10.0 | 384 |
| 4 MiB | 1.07 [0.98-1.13] | 4.04 | 10.0 | 96 |

**The comparison routes at dna_100m (cold cache, windows of the writers 1 MiB, tuple buffers 64 KB):**

| route | counted B/char (explicit read + explicit write + mapped write) | rchar + wchar (+ mapped) | block layer read + written | phase s | RssAnon MB |
|---|---|---|---|---|---|
| BWT pass (route (ii)) | 5.000 = 4.000 + 0 + 1.000 | 4.27 (+1.00) | 4.30 + 1.00 | 2.43 [2.36-2.97] | 10.2 |
| **Psi pass** (1 MiB windows) | **5.000** = 1.000 + 0 + 4.000 | 1.000 (+4.000) | 1.02 + 4.17 | 1.04-1.10 | 4.4-10.0 |
| ISA by DIST | 23.73 = 11.86 + 7.86 + 4.00 (P = 375; 7.86 B/char of tuples written and read back; scratch peak 7.86) | 19.9 (+4.00) | 4.04 + 4.02 (the tuple files were deleted before write-back: the page cache absorbed 15.7 B/char) | 1.94 [1.79-1.94] | 41.3 |
| Psi from the ISA (gather) | 8.000 = 4.000 + 0 + 4.000 (the random reads of the ISA through its mapping are not counted: block layer +4.0) | 4.48 (+4.00) | 8.05 + 4.02 | 4.48 [4.26-4.48] | 11.3 |
| **ISA route in total** | **31.7** | | | **6.4** | |

The ISA by DIST measures what the ADR said (about 24 B/char: 23.73 counted, the 1.3 % below 24 are the 375 write buffers that never reach their file). The sum of the two programs' RssAnon is not a peak of one process.
The 16 B/char of tuples never reached the device here, because 786 MB of scratch files that are written and read within two seconds fit the page cache of a 15 GB machine; at 3x10^9 the 24 GB of tuples would not (E7's cache findings apply).

## 4. G1 and G3 in detail

* **Files.** For every input of `e9_setup.tsv` the references are made by the in-RAM code of P1: `oracle` (SaisKernel suffix array in RAM, for dna_100m the SA of the E6 pipeline, whose SHA-256 and that of the packed text equal the
  hashes of the E4-E7 references in `/tmp/stbench/e7/dna_100m.refs.sha`: row `sa_pk_equal_e7_reference_hashes`), then `e9-oracle`: `LinkTables.BuildIsa` / `BuildPsi` in RAM written as files (psi[0] = -1 stored as 0). The BWT of the pass is made by `e9-bwt` and **checked independently**
  by `check_bwt_e9.py` (Python, no code shared: dense symbol of `T[SA[j] - 1]` for every row, the dummy byte, the meta file): 8 of 8 `bwt_check yes`, 0 mismatches (`e9_bwt_check.txt`). `cmp` of the Psi file with the reference: dna_10m, prot_1m, rep_10m, a_1m, fib_1m, tm_1m, rnd256_10m and dna_100m, five comparisons each (default parameters;
  `--window-kb 1 --symbols-per-pass 2 --read-kb 64`; the cyclic file's tail and entry 0; ISA by DIST against the ISA reference; Psi from that ISA): **40 of 40 `yes`**, and every one of the 90 timed runs of the protocol compared its output with the reference as well (**90 of 90 `yes`**, `cmp_output`).
  The PH0 counts equal the BWT's counts at dna_100m (`counts_equal_ph0` yes, three runs: the `.pk.meta` of the E6 pipeline lists them).
* **Tests** (`tests/SuffixTree/SuffixTree.EmProto.Tests/Ph6`, 51 tests; the whole-suite numbers are in section 8): the BWT file equals the definition for all nine alphabets with read buffers of 4-4096 bytes and windows of 1-65536 bytes; Psi equals the in-RAM Psi for 108 random texts over the nine alphabets with rotating windows (1, 4, 6, 12, 64, 4096, 65536 bytes),
  groups (1, 3, 7, 64, 256 symbols per pass) and both conventions of entry 0, for every combination of window, group and convention on three alphabets (sigma = 5, 65, 256, n = 700), for the default parameters at sigma = 64, 65, 129, 255, 256 (single pass and 2-4 passes), on all binary texts to length 9, ternary to 6 and 4-ary to 4,
  on the adversarial families at 10 lengths (1-400), on a^n, Fibonacci and Thue-Morse of 30,000 (4096-byte windows, 64 symbols per pass; 8-byte windows, one symbol per pass), the sentinel row at row 1, an inner row and row n (`abcdefg`, `gfedcba`, `banana`, ...), n = 0 and n = 1; the passes, windows and byte counts equal their formulas;
  a BWT that does not match its counts, a symbol outside the alphabet, a file of the wrong length, a window cap that is exceeded and an `ArenaBudget` that is released are checked; ISA by DIST equals the in-RAM ISA for chunk sizes 1 .. 2^30 and buffers of 8 bytes to 64 KB, its I/O volume is 6 words per entry when every tuple goes through a file and 1 word when none does,
  and a non-permutation is refused (duplicate, missing, out of range); Psi from the ISA equals the in-RAM Psi; the smoke test is run on the definition for all nodes of the adversarial families (9 lengths: 1, 2, 3, 5, 8, 17, 64, 100, 257), of all texts of the exhaustive enumeration (binary to 9, ternary to 5, 4-ary to 4: 1,725 texts), with the Psi of the stream, the Psi of the ISA route and the cyclic Psi, at two fan-outs, and it **fails** when two entries of Psi or of the ISA are swapped or when the definition is wrong
  (negative controls). The commands (`e9-bwt`, `e9-psi`, `e9-isa`, `e9-psi-isa`, `e9-oracle`, `e9-smoke`) are exercised through `Cli.Run`: exit codes 0 / 1 / 2, the rows the scripts read, the counters equal to their formulas.
* **G3 at scale** (`e9_identity.tsv`; one process per input, no step limit but the P1 budget of 3,000,000 per link; the whole smoke includes building the in-RAM index and the in-RAM tables of P1, so it runs without the 96 MiB heap limit): the SA file equals the SA of the in-RAM index, the ISA of DIST equals `LinkTables.BuildIsa`, the Psi of the stream and of the ISA route equal `LinkTables.BuildPsi`
  (entry 0 mapped), then 8 cells x (sample) links with both tables; **0 mismatches in every cell of every input** (7 inputs: dna_100m has no smoke) (the success flag, the target and the step count are all compared; no link exceeded the budget except the 4,968 of 5,000 descent links on a^n, in both runs, as P1 predicted).

| input | internal nodes (= the header statistic) | sample | cells | mismatches |
|---|---|---|---|---|
| dna_10m | 6,229,745 | 200,000 | 8 x 2 Psi routes | 0 |
| rep_10m | 8,568,247 | 200,000 | 8 x 2 | 0 |
| prot_1m | 295,123 | **145,561** | 8 x 2 | 0 |
| rnd256_10m (256 symbols) | 2,093,383 | 200,000 | 8 x 2 | 0 |
| a^n, Fibonacci, Thue-Morse, 10^6 | 999,999 / 999,995 / 999,902 | 5,000 (step budget 20,000) | 8 x 2 | 0 |

## 5. G4: memory, time, and the finding about window size

* **RssAnon.** Independent sampler (`sample_e9.py`, `/proc/<pid>/status` at 20 ms, whole process including the runtime and the thread stacks): Psi pass at dna_100m **10.9 MB** [10.87-10.97] with 64 KB windows and **10.0 MB** with 1 MiB windows (both include the in-process samplers' buffers; without them 5.28 and 4.38 MB); at dna_10m 10.0 and 9.6 MB, at prot_1m 5.3 and 4.3 MB;
  the 256-symbol alphabet 10.1 MB (64 cursors, 4 passes) and 10.2 MB (256 cursors x 64 KB). The runtime alone is 3.05 MB (`e9.base`). BWT pass 10.2 MB, ISA by DIST 41.3 MB (375 tuple buffers of 64 KB = 24 MB + 1 MB arrays), Psi from the ISA 11.3 MB. All far below 256 MB; and independent of n: the RssAnon of the Psi pass is 10.0 MB at n = 10^7 and 10.9 MB at 10^8 (64 KB windows; 9.6 and 10.0 MB with 1 MiB).
  The program's own sampler (10 ms) and the independent one (20 ms) agree to 0.1 MB at 10^7 and 10^8; at 10^6 they differ by up to 1.4 MB (the 10 ms sampler sees short peaks the 20 ms one misses).
* **Arena.** `psi_arena_peak` = 1,310,720 bytes with 64 KB windows and 5,242,880 with 1 MiB windows at both n = 10^7 and n = 10^8 (the BWT read buffer, 1 MiB, plus the reservation sigma x window = 4 windows): constant in n; at n = 10^6 the read buffer is its n + 1 bytes and sigma = 20 (prot_1m: 2,310,721 = 1,000,001 + 20 x 64 KB); 256 symbols: 5.2 MB (64 cursors x 64 KB), 17.8 MB (256 cursors x 64 KB), 68.2 MB (64 cursors x 1 MiB). The bound sigma x window (64 symbols x 64 KB = 4 MiB, x 1 MiB = 64 MiB, 256 x 64 KB = 16 MiB) is
  charged to the budget before the first window is opened and the pass refuses (`ArenaBudgetExceededException`) a configuration over the cap of 160 MiB (256 x 1 MiB is refused by a test).
* **Time.** Psi pass at dna_100m, cold cache, median of 3: 2.66 s [2.61-2.68] (64 KB), 1.10 s [1.08-1.13] (1 MiB), 10.4 ns per entry at 1 MiB without the in-process samplers. The BWT pass 2.43 s [2.36-2.97], ISA by DIST 1.94 s, Psi from the ISA 4.48 s. **ESTIMATE for 3x10^9, linear in n (x30), not measured, and assuming the device keeps
  the 400 MB/s it showed here (the 417 MB of the Psi pass reached the device within 1.04 s): the Psi pass 80 s with 64 KB windows, 33 s with 1 MiB windows (12.5 GB written, 3 GB read); the BWT pass 73 s; both together 106 s with 1 MiB windows; the ISA route 192 s (58 s ISA + 134 s gather, if the gather's random reads over a 12 GB file stayed as cheap as here, which they will not).**
  The abandon rule (projected > 5 min at 3x10^9) does not fire.
* **Finding: the window size decides the block-layer volume and the time.** `strace -c` of a Psi pass at dna_10m shows 613 `msync`, 855 `mmap`, 682 `munmap` for 612 windows: the .NET `MemoryMappedViewAccessor.Dispose` flushes the view, one synchronous write per window. With small windows the flush cost (about 250 microseconds per window here) and
  the extra device writes (journal and the partly filled pages at the window boundaries, which are written twice) dominate: 16 KB windows write 12.9 B/char for 4.0 of data, 64 KB windows 6.12, 256 KB 4.65, 1 MiB 4.17, 4 MiB 4.04. The same flush-on-dispose applies to every mapped writer of the prototype (E6's SA leaf slices are 133 KB views: not examined).
  The extra bytes were measured, their split between journal and double-written boundary pages was not. The ADR's sentence "sigma windows of 64 KB" (section 3.3) therefore costs 2 B/char at the block layer and 2.4x the time against the "sigma windows x 1 MB" of its PH6 row; **recommend windows of at least 1 MiB** (64 x 1 MiB = 64 MiB of mapped windows for a 64-symbol group; 256 symbols in 4 groups of 64).

## 6. The alphabet: more than 64 symbols

`rnd256_10m` (10^7 random bytes, sigma = 256, the BWT symbols are the 256 byte values and the sentinel is out of band): at the default of 64 symbols per pass the Psi pass takes **4 passes**, counted 8.000 B/char (4 x 1.000 + 4.000), 0.443 s [0.427-0.482]; with 256 cursors of 64 KB in **one pass** 5.000 B/char, 0.393 s [0.381-0.496]
(arena 16 MiB + 1 MiB); with 64 cursors of 1 MiB and 4 passes 8.000 B/char, 0.302 s (the window effect dominates the extra passes; rchar + wchar 4.38 B/char). The 6 B/char gate therefore holds for sigma <= 64 with one pass and for byte alphabets only with all 256 cursors in one pass (256 x 64 KB = 16 MiB, or 256 x 256 KB = 64 MiB, inside the 160 MiB cap);
at 4 passes it is 8 B/char. Identity holds in every variant (`rnd256_10m.id`: 5 of 5, and the 3 + 3 timed runs equal the reference).


## 7. Defects found and findings

1. **Product code (found by a test written for it, red before the fix):** `LinkSmoke.Run` with the P1 sampling drew from an empty node set when the text has no internal node (`a`, `ab`, `abc`): `ArgumentOutOfRangeException` in `drawn[..half]`. Fixed (the sample is taken only when nodes exist); the measured inputs are not affected.
2. **Measuring script:** the first run of the smoke group died with "Out of memory" under `DOTNET_GCHeapHardLimit` = 96 MiB (the smoke test holds the in-RAM index and the in-RAM tables of P1: it is a test of the kernels, not a bounded-memory run). `run_e9_identity.sh` now unsets the limit for the smoke command only; all rows of `e9_identity.tsv` are from the corrected script.
3. **Finding (not a defect of correctness):** windows below 1 MiB make the pass slower and amplify block-layer writes (section 5). The 64 KB windows of the task statement miss the block-layer part of gate G2.
4. **Finding:** the DIST tuple files of the ISA route (7.86 B/char written and read back) did not reach the device at 10^8 (page cache); the block layer under-reports the route's I/O by 15.7 B/char at this size, which is why the program counters are the primary figure.
5. **Observation, not changed:** `MemoryMappedViewAccessor.Dispose` calls `msync` for every view; E6's bounded dirty pages (< 13 MB for 400 MB of SA) are a consequence of that, not of the kernel's write-back policy.

## 8. Tests, builds, mutation audit

* **Tests first:** 41 tests written before any implementation against a stub whose methods throw `NotImplementedException`: **41 of 41 red** (`c9c84557`); after the implementation 41 of 41 green (`67a5dec7`). The stub is weak evidence (everything fails at the first call); the evidence of test strength is the audit below. Added afterwards, not red-first: the 5 tests of the commands and one more
  refused-SA case (`a7dacc81`), 3 tests with the smoke fix (`a995fa14`; `ATextWithoutInternalNodes_HasAnEmptySample` red before the fix), 1 test of `--no-sampler` (`8d79fc47`): **51 Ph6 tests** (4 more tests after the audit of round 1 and a stronger assertion after round 2: `TheKindsOfMismatch_AreSeparated_AndTargetOnlyDifferencesExist`, the overflow / window / `N - 1` cases and the longer-SA case in existing tests).
* **Existing code:** the e6 pipeline, the existing commands and `src/**` are unchanged. Changes to existing files: `Cli.cs` (usage lines and six `case` lines), `SuffixTree.EmProto.csproj` (six `Compile Include` lines linking P1 files unchanged), `SuffixTree.EmProto.Tests.csproj` (two `Compile Include` lines). `apps/SuffixTree.Esa.Bench` is untouched.
  Route (ii) was chosen because it needs no flag in the SA sink: there is nothing to prove "flag off unchanged" about, and the whole suite is the proof that nothing else changed.
* **Default build, whole EmProto suite:** **392 of 392** (341 existing + 51 Ph6; 12 m 42 s) on the tree of commit `a9c5ede2` (the final code: later commits change only documents).
* **Checked build** (`-p:EsaBounds=true`, `ESA_EXPECT_BOUNDS=true`): the new code allocates through `BudgetedArena` / `Buf<T>` (which the checked build compiles differently), so the Ph6 tests were run in it: **51 of 51 green** (1 m 18 s, run on the final tree). The rest of the suite was not run in the checked build (no code of the other stages changed).
* **Mutation audit** (hand-made single-site mutants of the new production logic: the sigma-way distribution, the C[] counts, the sentinel row, the multi-pass loop, the BWT definition and its counts, DIST, the gather, the file conventions, the smoke driver; `analysis/mutation/make_e9_mutants.py` writes them, `run_e9_mutants.py` builds and tests each with a watchdog in the Release configuration: stage 1 = the Ph6 tests without the slow ones, a survivor is run against all Ph6 tests but the 30-second exhaustive smoke; the source is restored after every mutant):
  * **70 mutants written and run** (round 1: 70 written against the sources of `8d79fc47`: `PsiStream.cs` 36, `BwtPass.cs` 10, `IsaByDist.cs` (DIST and the gather) 15, `PsiFiles.cs` 3, `LinkSmoke.cs` 6; 4 of them did not compile (analyzer errors S1125, S1481, S1656) and 1 pattern was stale: 65 were built and run; round 2: those 5 as replacements or corrected patterns plus the survivors that a test could kill, 11 runs; round 3: P16, 1 run).
  * **Final status: 64 killed, 4 equivalent (justified), 2 survivors not killed, 0 not run.** Round 1: **55 killed, 10 survived** (two of the kills are "killed(host crash)": a wrong cursor wrote outside its window and the test host died; a crash counts as a detection). Of the 10 survivors, 4 are equivalent (P20, B08, D07, D08: reasons in the table), **5 were killed after the tests were strengthened** (P14: the overflow message is asserted; P16: windows of one and two entries and the message of the end-of-pass check; P24: counts that sum to more than n; B10: an SA file longer than n + 1 entries; S03: the smoke cells now count target-only differences and a test finds them) and **2 stay alive**: S04 (an exception in the kernel is not counted as a mismatch by itself: it already differs from a successful reference link, the term matters only when the reference link gave up and the other threw, which no test constructs) and S02 (step counts compared with <= instead of !=: the corrupted tables of the tests never give a link that differs by more steps alone). Both are in the smoke driver, not in the pass.
  * The 4 equivalent mutants were not run against a changed test; the claim is by argument (table). No candidate was left written and not run.
  * Mutation caveat: the round-1 kills were obtained with the tests of `8d79fc47`; the later rounds only added assertions, and the whole suite and the Ph6 tests were re-run at the end (see above).
  * The full table of 70 is in section 11; the raw rows are `analysis/mutation/results_e9_*.tsv`, `results_e9_round2_*.tsv`, `results_e9_round3_*.tsv`.


## 9. Not verified / limits

* Everything at 3x10^9 is an ESTIMATE (linear extrapolation x30 of single-machine, 100 Mbp, cold-cache numbers; the 12 GB Psi file and the 24 GB of DIST tuples at that size would not fit the page cache; the device kept 400 MB/s at 400 MB).
* **Route (i), the BWT written by PH3 "for free", is not built and not measured.** The cost of the BWT here is a separate pass (5.0 B/char counted, 2.43 s at dna_100m); what the fused route would cost (one random text read and one byte written per suffix inside the leaf sorter, in the hottest loop of the pipeline) is an argument, not a number.
* The Psi pass is a standalone program over files; it was **not** fused with PH5 or with PH4's B3 hook, and was not run inside the pipeline's process (RssAnon of the whole pipeline with PH6 is not measured). Windows < 1 MiB were measured only on DNA (sigma = 4) at 10^8 and for 256 symbols at 10^7.
* The BWT is one byte per row: alphabets above 256 symbols (the 32-bit names of the DCS recursion) are not supported by this pass; positions above 2^32 - 3 (40-bit positions) are refused; the PSI file is a standalone file, **not a section of the container** (no section offset, no CRC).
* Only the P1 link kernels consume the Psi (the smoke test); `RankOf` / `ExtendLeft` and the navigation pack's RMQ / PSV / NSV are not built here. The in-RAM Psi of P1 is the oracle of the files (G1); the definition (rank interval of the label minus its first symbol) is the oracle of the kernels only on small texts.
* The block-layer figures are system-wide sector counters of `vda` minus a baseline run; below 10^7 symbols fixed reads of 1 MB or more distort them (prot_1m); the split of the extra device writes of small windows into journal and double-written pages was not measured; one disk (virtio, probably host-cached), one machine, 3 runs per cell.
* Mutation audit: hand-made single-site mutants; the equivalent ones are justified one by one, a "killed" means that the stage-1 or stage-2 test run failed (or hung or crashed the test host) for that mutant in this run.

## 10. Open points for 1b.7 (the container writer and the pack)

1. **Where PSI sits in the container.** ADR-003 kind `PSI` (id 12 proposed, 4 B/char, `PSI[r] = ISA[SA[r] + 1]`, `PSI[0] = 0`): its size is known after PH0 (4 (n + 1) bytes), so its offset is assigned up front with TEXT, SA, CHILD; the pass writes the symbol cursors at `sectionOffset + 4 C[c]` through windows of **at least 1 MiB** (section 5; `MemoryMappedViewAccessor` per cursor, or `RandomAccess.Write` of 1 MiB buffers, which this prototype did not try), the section CRC-32C is one more read pass. P1 keeps -1 in RAM for entry 0: the reader maps 0 to -1 (`PsiFiles.LoadForP1`), no kernel reads it.
2. **PSI[0] convention.** The file holds 0 (ADR-003 section 4). The cyclic variant (entry 0 = the row of the suffix 0, one cycle through all rows) costs nothing and is a free consistency check of the sentinel row (`--cyclic`); the BWT's `primary` row (the row of the suffix 0) must travel from PH3 to PH6 (one number).
3. **Multi-pass for large alphabets.** 64 cursors per pass means ceil(sigma / 64) scans of the BWT (8 B/char at sigma = 256, measured); with 1 MiB windows 64 cursors are 64 MiB. The one-pass variant with 256 cursors needs 256 x (64-256 KB) = 16-64 MiB and meets the 6 B/char gate for every byte alphabet; pick the group size from the window size so that group x window stays well under the 160 MiB arena cap. The 32-bit alphabets of the DCS names do not occur in PH6.
4. **How the BWT is produced in the streaming build.** This prototype: a separate pass, SA read 4 B/char + BWT written 1 B/char + random reads of the packed text (0.25 B/char file, resident), 2.4 s at 10^8, 5.0 B/char in the program counters, 4.3 + 1.0 at the block layer. Whether PH3 should write the BWT instead is open: the leaf sorter holds the sorted positions before `sink.Write`, so `T[pos - 1]` is one random read of the mapped text per suffix and the byte goes to the BWT scratch file at the same offset as the SA slice (the BWT is scratch: it can be deleted after PH6), at the price of one more mapped writer and 1 B/char of scratch at the PH3 peak (4.3 B/char + 1). The row of the suffix 0 is known there (the leaf that writes SA value 0).
5. **Can the pack be built fused with PH5?** In principle yes for alphabets of at most one group: PH6 consumes the BWT in rank order exactly as PH5 consumes the LCP in rank order, so both can hang on the chunks of PH4's B3 hook (the BWT chunk read next to the LCP chunk: +1 B/char of reads, +4 B/char of writes through 4-64 more windows, about 10 ns per entry on B3's single thread against 23 ns for the PH5 scan, ESTIMATE from 1.04 s per 10^8 here and E5's number; the hook already runs on one thread, so the two scans would share it or need the second thread that E6's open point (a) asks for). Not built, not measured. A larger alphabet would need the re-reading passes or a second hook consumer per group.
6. **The window size is a first-order parameter of every mapped writer** (section 5): 1b.7's `MappedWindowWriter` should take 1 MiB or more by default and the cost of flush-on-dispose should be known before the CHILD, SA and LCP windows are tuned.
7. **ISA stays optional** (`RankOf` / `ExtendLeft`): DIST costs 23.7 B/char and 1.9 s per 10^8 here; the gather Psi-from-ISA is only a comparison (31.7 B/char against 5.0), and at out-of-core sizes the 12 GB ISA cannot be gathered from cheaply (ESTIMATE: another DIST, about 48 B/char in all).

**Cost of the work (measured only as wall-clock between commits; tokens are an estimate, not a counter reading):** about 2.3 hours of work in two sessions (07:25-09:05 and 12:02-12:40 UTC; the 5-hour limit fell between them, no measurement was running then), about 0.5-0.6 M tokens (estimate); the economy cap of the brief (forecast 1-2 h, 0.3-0.6 M tokens, stop at 3 h or 0.8 M) was not reached.

## 11. The mutants, one by one (`tables_e9_mutation.py`)

| file | id | mutant | final status |
|---|---|---|---|
| BwtPass | B01 | the symbol before the suffix is taken at p instead of p - 1 | killed |
| BwtPass | B02 | the sentinel row is the row of suffix 1 | killed |
| BwtPass | B03 | counts add two | killed |
| BwtPass | B04 | row 0 is not checked to be the empty suffix | killed |
| BwtPass | B05 | primary row off by one | killed |
| BwtPass | B06 | symbol equal to sigma accepted | killed |
| BwtPass | B07 | an SA entry equal to n + 1 accepted | killed |
| BwtPass | B08 | two rows with suffix 0 accepted | equivalent (round 1: survived): row 0 must hold n (checked), so the first row that holds 0 is row >= 1 (for n = 0 there is one row): "primary >= 1" and "primary >= 0" agree whenever a second row holds 0 |
| BwtPass | B09 | the dummy byte of the sentinel row is 1 | killed |
| BwtPass | B10 | SA file of the wrong length accepted (<) [a longer file is tested now] | killed (round 2, after round 1: survived) |
| IsaByDist | D01 | chunk of a position off by one | killed |
| IsaByDist | D02 | rank in the tuple off by one | killed |
| IsaByDist | D03 | buffer entries from 4-byte tuples | killed |
| IsaByDist | D04 | chunk count by floor division | killed |
| IsaByDist | D05 | the chunk array is not initialised to the missing marker | killed |
| IsaByDist | D06 | the last chunk array is as long as a full chunk | killed |
| IsaByDist | D07 | duplicate positions are not detected when they occur | equivalent (round 1: survived): rows entries in [0, n] without a missing position are a permutation; a duplicate always leaves a position missing, which the missing-position check refuses (same exception type; only the message differs) |
| IsaByDist | D08 | missing positions are not detected | equivalent (round 1: survived): the converse: with the duplicate check in place, no duplicate and rows entries in [0, n] leave nothing missing |
| IsaByDist | D09 | entries beyond n accepted | killed |
| IsaByDist | D10 | the tuples of a flushed buffer are written at the start of the file | killed |
| IsaByDist | D11 | the chunk files are read from the wrong offset | killed |
| IsaByDist | D12 | the in-memory remainder of a chunk is not scattered | killed |
| IsaByDist | F01 | Psi reads ISA[SA[i]] instead of ISA[SA[i] + 1] | killed |
| IsaByDist | F02 | Psi entry 0 is 1 | killed |
| IsaByDist | F03 | SA entry n at a row >= 1 accepted | killed |
| LinkSmoke | S01 | sample A of P1: half without the -4000 | killed |
| LinkSmoke | S02 | step counts compared one way only (<=) | SURVIVOR, not killed: the step counts are compared with <= instead of !=; the tables the tests corrupt never give a link that differs only by more steps in the table under test |
| LinkSmoke | S03 | targets are not compared | killed (round 2, after round 1: survived) |
| LinkSmoke | S04 | kernel exceptions are not counted as mismatches | SURVIVOR, not killed: an exception leaves okB = false, which already differs from a successful reference run; the term matters only when the reference run gave up and the other run threw, which the tests could not construct |
| LinkSmoke | S05 | the empty-node guard of the P1 sampling is dropped | killed (round 2, after round 1: PATTERN-NOT-UNIQUE) |
| LinkSmoke | S06 | the definition is compared for equality instead of difference | killed (round 2, after round 1: BUILD-ERROR) |
| PsiFiles | L01 | LoadForP1 keeps entry 0 as the file has it (0): the mutant maps to 0 instead of -1 | killed (round 2, after round 1: BUILD-ERROR) |
| PsiFiles | L02 | the oracle file stores P1's -1 in entry 0 | killed (round 2, after round 1: BUILD-ERROR) |
| PsiFiles | L03 | a ragged file is accepted | killed |
| PsiStream | P01 | C[] base without the empty suffix: start[0] = 0 | killed |
| PsiStream | P02 | C[] ignores the counts: start[c+1] = start[c] | killed |
| PsiStream | P03 | C[] adds one per symbol: start[c+1] = start[c] + counts[c] + 1 | killed |
| PsiStream | P04 | sentinel row: skip index off by one (primary - first + 1) | killed |
| PsiStream | P05 | sentinel row at the start of a block is not skipped (skip > 0) | killed |
| PsiStream | P06 | rows after the sentinel get the row number of the sentinel (firstRow of the second scan without +1) | killed |
| PsiStream | P07 | the first scan includes the sentinel row (skip + 1) | killed |
| PsiStream | P08 | the second scan starts at the sentinel row (no +1 in the pointer) | killed |
| PsiStream | P09 | symbol of the pass off by one (bytes[j] - lo + 1) | killed |
| PsiStream | P10 | row written off by one (firstRow + j + 1) | killed |
| PsiStream | P11 | cursor advances by two entries | killed |
| PsiStream | P12 | the cursor window is the larger of the window and the rest (Max) | killed |
| PsiStream | P13 | the next window of a cursor starts one entry later | killed |
| PsiStream | P14 | a cursor that has used its range is not refused (remaining < 0) [message of the overflow is asserted now] | killed (round 2, after round 1: survived) |
| PsiStream | P15 | under-full cursors are not detected inside the last window | killed |
| PsiStream | P16 | windows not fully mapped are not detected [windows of one entry are tested now] | killed (round 3, after round 2: survived) |
| PsiStream | P17 | mapped bound without the factor 4 (entries instead of bytes) | killed |
| PsiStream | P18 | passes by floor division: the last partial group is dropped | killed |
| PsiStream | P19 | the last group runs past sigma (hi without Math.Min) | killed |
| PsiStream | P20 | the group start uses the open cursors instead of the limit (pass * group) | equivalent (round 1: survived): the group start pass * SymbolsPerPass and pass * group differ only when sigma < SymbolsPerPass, and then there is one pass (pass = 0) |
| PsiStream | P21 | convention swapped: entry 0 is the primary row unless cyclic | killed |
| PsiStream | P22 | cyclic entry 0 off by one | killed |
| PsiStream | P23 | primary row == n + 1 accepted | killed |
| PsiStream | P24 | counts that sum to more than n accepted [N - 1 is tested now] | killed (round 2, after round 1: survived) |
| PsiStream | P25 | a BWT file longer than n + 1 accepted | killed |
| PsiStream | P26 | rows = n | killed |
| PsiStream | P27 | zero symbols per pass accepted | killed |
| PsiStream | P28 | window below four bytes gives zero entries | killed |
| PsiStream | P29 | read buffer not bounded by the number of rows | killed |
| PsiStream | P30 | the BWT is read from the wrong offset | killed |
| PsiStream | P31 | window count of the passes is not accumulated | killed |
| PsiStream | P32 | BWT bytes read reported for one pass | killed |
| PsiStream | P33 | the window cap is not charged (reservation of a different size) | killed |
| PsiStream | P34 | the shared budget is not released | killed |
| PsiStream | P35 | symbols of other passes are counted as this pass (s > span) | killed |
| PsiStream | P36 | the first entry of the cursor range is entry start (without the +1 of the empty suffix) via counts offset | killed |

distinct mutants with a final status: 70; killed 64; equivalent (justified) 4; survivors not killed 2; never run (no compile / stale pattern, not replaced) 0

