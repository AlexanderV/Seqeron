# E4: PH4 (the LCP array from a given suffix-array file with bounded RAM) - measured results

Work package of the owner decision of 2026-11-11 (ADR-004 section 6c authorises E4-E9; section 5, row E4 defines the prototype: "PH4 (chunked Phi/PLCP via DIST) on a given SA file", gate "LCP equal to `LcpKernel`; I/O <= 50 B/char; throughput >= 5 MB/s on 4 threads", abandon "LCP phase > 3x the SA phase"). The gates G1-G4 below were fixed by the orchestrator before the work.
Code: `apps/SuffixTree.EmProto/Ph4` (`ExternalLcp.cs` entry point, options and stats; `LcpBuild.cs` the three phases in both modes; `PlcpScan.cs`; `Spool.cs` = the generalised DIST primitive; `LcpOutput.cs`; `LcpGeometry.cs`) and `Experiments/E4.cs` (commands `e4`, `lcp-oracle`, `check-lcp`); tests `tests/SuffixTree/SuffixTree.EmProto.Tests/Ph4`.
Raw rows: `e4_*.tsv` in this directory (`R <TAB> E4 <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes`; `e4_setup.tsv` also holds the E3 rows of the SA phase), scripts `run_e4_setup.sh`, `run_e4_protocol.sh`, `run_e4_big.sh`, `run_e4_cgroup.sh`, `run_e4_cgroup_sweep.sh`, `check_allsame.py`, `summarize_e4.py`; `e4_rounds8.tsv` (8 rounds of the bounded mode, run with the same binary and options by hand after the protocol), `e4_cgroup.tsv` with `e4_cgroup_stats.txt`.
Everything below is **measured here** (4 vCPU Intel Xeon, 15 GB, kernel 6.18.44, .NET 10.0.14, Release, workstation non-concurrent GC, nothing else running, `DOTNET_GCHeapHardLimit` = 96 MiB, every timed run in a fresh process after `sync; echo 3 > /proc/sys/vm/drop_caches`, i.e. a **cold page cache**, 4 threads unless noted) unless it is marked ESTIMATE. The binary under test of every measurement is the build of commit `5e93c92b` (the head of each protocol is in the rows as `git rev-parse HEAD`); the later commits change tests, scripts, results and one reporting helper of `E4.cs` (`E4.Healthy`, the exit-code condition, moved into a function so that it can be tested) - nothing on the measured path.

## 1. What PH4 does here

Input: the packed text (read through a read-only mapping) and the suffix array file of the E3 builder (`n + 1` little-endian `uint`, `SA[0] = n`). Output: the compact LCP in the form of `CompactLcp`: a base file (`n + 1` entries of 1, 2 or 4 bytes, little endian, escape value 255 / 65,535 in the narrow kinds) and an overflow file (sorted pairs `(rank, value)` of two 4-byte fields, ADR-003 section 4 kind 5), plus the counts of escapes (`ExtLcpStats`). Chunk size `R` = 8,000,000 positions (ADR-004 4.2), `P = ceil(n / R)` chunks (375 at 3x10^9).

**Base mode ("Distributed", the ADR's data flow).**
* **B1** scans the SA in rank order (explicit `pread`, 1 MiB blocks) and writes the tuple `(SA[i], SA[i-1], i)` (12 B) to the buffer of chunk `SA[i] / R` (division by `FastMod32` of E3, `P` buffers of 32 KiB, one arena block); a full buffer is appended to the chunk file `ph4phi-<chunk>.bin`. One thread.
* **B2** takes the chunks in any order, several in parallel: reads the chunk file, scatters `Phi` and the rank into two arrays of `R` entries, **deletes the chunk file**, runs the PLCP loop of `LcpKernel.ComputePlcp` on the chunk (`PlcpScan.Run`: carried length `l`, `PackedText.Lce` = the 256-bit vector compare, prefetch of the text line 16 positions ahead, counts of values >= 255 and >= 65,535 in the loop), then writes `(i, PLCP[j])` (8 B) to the buffer of rank chunk `i / R`.
* **Width.** When B2 has finished, `LcpWidth.Select(n, e8, e16)` (the rule of ADR-002 section 4 on 64-bit counts, pinned to `CompactLcp.SelectKind` by a test) chooses U8 / U16 / U32; **the output files are created only now**, with their final lengths (`(n + 1) x width` and `8 x overflow entries`), so the width is known before any LCP entry is written.
* **B3** takes the rank chunks in ascending order, scatters `(i, value)` into one dense `uint[R]`, calls `ILcpChunkConsumer.Consume(firstRank, values)` (the stream point of PH5: dense resolved values in strictly ascending rank order), and writes the base array and the overflow pairs through `MappedWindowWriter` (1 MiB sliding windows, no explicit flush; the owner's rule "stream into memory-mapped files"). One thread.

**Bounded mode ("Bounded", ADR-004 4.2, bounded-scratch paragraph).** The text positions are cut into rounds of `RoundChunks x R` positions. Per round: **B1'** reads the whole SA again and keeps the 8-byte pairs `(SA[i], SA[i-1])` of the positions of the round (no rank is needed), **B2'** is B2 without the rank array (a task holds `R` Phi entries only, so four tasks fit the cap) and writes the PLCP chunk (4 B per position) at its final offset of one text-order file by `pwrite`. After the last round the width is chosen as above; **B3'** reads the SA once more in rank order and gathers `LCP[i] = PLCP[SA[i]]` from the PLCP file through a read-only mapping (software prefetch 16 entries ahead), then the same consumer hook and writer.

Output independent of thread count, `R`, buffer sizes, mode and rounds (owner answer Q5): the files are a function of the SA and the text only; tested (section 3) and measured (section 4).

## 2. Design choices the ADR does not fix (reason, alternative)

| choice | reason | alternative not taken |
|---|---|---|
| Tuples as in the ADR: `(SA[i], SA[i-1], i)` 12 B, `(i, PLCP)` 8 B, positions and ranks as `uint` | the ADR's I/O sum; 32-bit positions as E3/E3b (n < 2^32 - 1) | 11-byte packed tuples, 40-bit positions (no 3x10^9+ text is in the prototype) |
| One scratch file per chunk, **deleted as soon as its consumer has it in memory** (B2 deletes the Phi file after the read, B3 the rank file after the chunk) | owner answer Q1: scratch released stage by stage; the peak is then exactly the end of B1 (12 B/char) | one file per distribution with a region per chunk (peak 12 + 8 = 20 B/char: the ADR's "12n + 8n") |
| Always through files, no "small chunk stays in memory" shortcut of `SaInverter` | one code path for every size; the tests exercise it with 12-byte buffers | the shortcut (a branch that only small inputs take) |
| B1 and B3 single-threaded, B2 parallel by chunk | B1 1.0-1.2 s and B3 0.8 s per 100 Mbp: far from a gate; no per-thread buffer sets (`P x 32 KiB` each) | B1 parallel over rank ranges with one buffer set per thread (+12 MB per thread at P = 375) |
| **B2 tasks = what the 160 MiB arena cap allows at the reference geometry**: a base-mode task holds `Phi` + rank (2 x 4 B x R = 64 MB at R = 8 M) + the tuple buffers (376 x 32 KiB) = 73.8 MiB, so **2 tasks** (bounded mode: 31.5 MiB, min(threads, 5) = 4) | the cap is enforced from inside (`ArenaBudget`); the scaled runs use the task count of the reference geometry (`--b2-workers`), not what their small arrays would allow | R = 4 M (P = 750) to run 4 tasks (the task fixes P = 375); a larger cap; the 4 threads are used by B2 only partly (2 of 4) |
| The PLCP loop restarts at every chunk start with `l = 0` (no carry between chunks) | chunks are independent tasks in any order; ADR: restart cost <= max PLCP per chunk | a serial carry (no parallelism), a fix-up pass |
| Width from the **global** escape counts after B2; per-rank-chunk counts are not kept | B3 is sequential, a running count is enough; the rule needs only the totals | per-chunk counts (needed only for a parallel B3 with disjoint overflow offsets) |
| Overflow file = interleaved pairs `(rank, value)`, 4 bytes each (ADR-003 kind 5); the base array little endian | the container layout, so 1b.7 can adopt the files | two parallel arrays (the in-RAM `CompactLcp` form; the test decoder turns pairs into it) |
| Bounded mode keeps the PLCP file as **U32, 4 B/char**, not the ADR's compact 1.1 B/char | the width of the final LCP is unknown until the last round, so the gather source must hold any value; U32 needs no escape structure and gives one random 4-byte read per rank | compact PLCP (U8 + escape + sorted escape list): a gather of an escaped value is a binary search, `n log n` random reads on a^n, where nearly every value escapes; **the price of the choice is 3.6x the page cache the ADR assumed for the gather (section 7)** |
| Bounded mode: pairs of 8 B (no rank), `--rounds r` = `ceil(P / r)` chunks per round | the gather re-reads the SA, so the rank is redundant: scratch `8 rho + 4 n` instead of the ADR's `12 rho + 1.1 n` | 12-byte tuples as in the ADR text |
| `ArenaGuards` option (guard bands and poison of the arenas) on `ExtLcpOptions` | test hook: a buffer read before it is written shows as 0xCDCDCDCD; all oracle tests run with it | none (null = production arena) |
| Scale-down emulation: `R = ceil(8,000,000 n / 3x10^9)`, so `P = 375` for every n; fixed buffers (32 KiB tuple, 1 MiB read, 1 MiB window) not scaled | as E3b: the count of chunks and rounds is that of 3x10^9, the arrays that grow with R scale with n; a fixed buffer does not depend on n | scaling the buffers too (they would no longer be the real ones) |
| No resume, no state file, 32-bit positions | as E3b; out of the prototype | ADR-004 4.3 crash safety (1b.7) |

## 3. Verification

| item | result |
|---|---|
| EmProto tests before the work | 112, all green |
| tests first, against the stub (`eba46d43`) | **23 of 135 red** (21 `NotImplementedException`, 2 expected `ArgumentException`), the 112 older green; the explicit fuzz class is not in the count |
| after the implementation | 23 green; with the component tests written after the code (labelled in the file), the tests added before the measurements and the 8 added after audit round 1: **180 green in the default run (1 m 49 s) and 180 green in the checked variant** (`-p:EsaBounds=true`, `ESA_EXPECT_BOUNDS=true`, 2 m 04 s); 112 -> 180, the 112 older tests among them unchanged |
| exhaustive (default run, both modes, R = 3 / 4 / 5 / 7 / 2, rounds of 1-2 chunks, threads 1 / 2 / 4, compared with `LcpKernel` in compact form and with `LcpReference` as decoded values) | every binary string of length 1-12 (8,190, base mode), 1-11 (4,094, bounded mode), every ternary string of length 1-6 (1,092, both modes): **14,468 builds**, all equal |
| differential, adversarial families at n = 1,000 and 5,000 (the families of `StringFamilies.Adversarial`, 2 random parameter sets each, more than 100 cases: R = 3..50, b in {min, 2, 3, 4, 8, 16, 32}, vector and scalar compare, threads 1 / 2 / 4, rounds 1-4, read / tuple / window buffers of 12-4,000 bytes, touched buffers) | all equal; asserted: > 90 % have several chunks, > 30 per mode, > 10 with several rounds, kinds U8 and U16 occur |
| widths (n = 6,400 / 140,000 / 150,000, each in both modes with chunk sizes 300-20,000, six cases) | U8 with 100+ overflow entries, U16 with 4,466, U32 without; equal to `LcpKernel`; the planted repeats are asserted to give those kinds |
| output independence (4 texts, threads 1 / 2 / 4, R in {3, 17, 200, 2^20}, buffers 12 / 100 / 4,096, both modes, rounds of 1 / 2 / 4 chunks) | base and overflow files byte-identical to one run and to each other |
| `Fuzz` class (explicit, texts 10^3..10^5 plus X X blocks and a^n of 140-170 k symbols that force U16 and U32; R = 3..50, many chunks, both modes) | seed 1, 120 s: **144 cases (78 bounded), kinds U8 35 / U16 80 / U32 29, all equal**; seed 2, 60 s: 89 cases (50 bounded), U8 21 / U16 51 / U32 17, all equal |
| written after the code (labelled in the test file) | chunk arithmetic against the division, the width rule against `CompactLcp.SelectKind` (grid and boundaries), the PLCP loop against the definition for every chunking and the restart cost of a^n, the spool, the writers, the output writer, the options and the task budget, the scratch peaks (12 / 8 / 6 B/char exactly), non-permutations, the three commands |
| expectations of my own tests that were wrong in the first version and corrected after the first run (the code was right) | `Select(160, 11, 20)` has more escapes16 than escapes8 (invalid); the bound on a^n's compare steps; the failure of a second escape is raised by `Write`, not by `Complete`; the bounded worker fit at the reference geometry is 5 (the command takes min(threads, 5)); the rank chunk count at 3x10^9 is 376 (n + 1 ranks); the corruption tests of `check-lcp` need a^70,000 (a^30,000 is U16 without overflow). One code change came out of a test: the mode of `e4` is validated before the text's meta file is read |

## 4. Gate G1: the LCP equals `LcpKernel`

Every protocol run compares the output with the files of the in-RAM `LcpKernel` (`lcp-oracle`: `LcpKernel.Build` then `CompactLcp` written in the same two-file layout; `--verify`, byte for byte). **360 of 360 verified runs equal, none different** (cold cache, fresh process each):
* the protocol groups on 23 inputs and families: 100 Mbp random DNA and genome-like (36 runs: threads 1 / 2 / 4, both modes, rounds 1 / 2 / 4, `R / 2` and `2 R`, i.e. P = 750 / 375 / 188), 10 M inputs and a^n 50 M (126 runs, the same set), the committed small families a^n, Fibonacci, Thue-Morse at 10^5 and 10^6, a 171-bp tandem, an N run, `dna_1m`, `prot_1m` (80 runs), the slope series 10 / 50 / 100 / 200 M random DNA and 10 / 50 / 100 M genome-like (84 runs): **326 runs**;
* 500 Mbp random DNA at the unscaled parameters (R = 8,000,000, P = 63): 11 of 11 equal to the kernel (the oracle took 36 s there and needs about 5 GB of RAM, not measured);
* a^n 500 Mbp: the kernel would need about 6.5 GB (estimate), so every run is compared with the **closed form** (`check_allsame.py`: the SA is `n, n-1, ..., 0` and `LCP[r] = r - 1` for every one of the 500,000,001 ranks): 11 of 11 equal; a first pass of the same configurations (`e4_a500_firstpass.tsv`, same timings within the run-to-run spread) was checked less strictly: the sampled `check-lcp` (2,002 ranks, structure scan) on the base output, a `cmp` of the 2-round bounded output against it, and the closed form on the last output, 0 differences;
* the cgroup sweep (section 7): 6 of 6 equal; 8 rounds of the bounded mode (10 M and 100 M, 3 runs each, `e4_rounds8.tsv`): 6 of 6 equal.

Widths seen: **U8 + overflow** (random DNA 0 escapes, genome-like 1.39 % of the values >= 255 and 0.50 % >= 65,535, 1,390,028 overflow pairs at 100 M), **U16** (`rep_10m`: 6.94 % >= 255, none >= 65,535, so U16 without overflow; the 1 M N run: 9.9975 %, U16 + 934,465 overflow pairs at 10 M), **U32** (a^n, Fibonacci 10 M, Thue-Morse 10^6, tandem: 99.5-99.9995 % of the values >= 255; the Fibonacci and Thue-Morse words of 10^5 symbols are U16).

**G1: PASS.**

## 5. Gate G2: RssAnon <= 256 MB including the runtime and thread stacks, 100 Mbp, 4 threads

Scale-down emulation (derivation: the 3x10^9 reference has `R = 8,000,000` positions per chunk and `P = 375` chunks; `R` is multiplied by n / 3x10^9 and rounded up, so `P = 375` for every n; 100 Mbp: `R = 266,667`; bounded mode: `--rounds r` gives `ceil(375 / r)` chunks per round, so the round count is that of 3x10^9; the number of parallel B2 tasks is the one the 160 MiB arena cap allows at the reference geometry, 2 in the base mode and 4 in the bounded mode; the fixed buffers are the real ones):

| input (3 runs each) | base mode | bounded, 2 rounds | bounded, 1 round | cap |
|---|---|---|---|---|
| dna_100m (random DNA) | **50.7 MB** [50.6-50.8] | 34.6 [34.4-35.7] | 33.7 [33.6-33.9] | 256 |
| genome-like 100 M | **50.8 MB** [50.8-50.8] | 34.5 [34.5-35.3] | 33.9 [33.7-33.9] | 256 |

The worst RssAnon of any of the 326 scale-down runs is **73.2 MB** (100 Mbp with `R / 2`, i.e. P = 750 chunks: the tuple buffers double). Mapped file pages (RssFile: 54-68 MB in the base mode; **450-470 MB in the bounded mode, the PLCP file read by the gather**) are reported, not counted (ADR-004 6a).

**Unscaled real-size runs** (R = 8,000,000, 32 KiB tuple buffers, 1 MiB windows, arena cap 160 MiB, 4 threads, cold cache; section 7): 500 Mbp random DNA, base mode **154.0 MB** [153.9-154.0] (arena peak 134.2 MB), bounded 152.8-155.8 MB; a^n 500 Mbp 153.8 MB [153.6-153.8] (base) and 153.3-155.3 MB (bounded, 4 threads). With 4 B2 tasks instead of 2 the cap is enforced from inside: `ArenaBudgetExceededException` "arena budget 167772160 bytes exceeded by 'b2.phi' (32000000 bytes; 162097144 live)" (`e4_big_enforce_*.txt`, exit code 134 of the unhandled exception).

**Slope of RssAnon against n** (least-squares fit of the RssAnon peak over log2 n, median of each size, in brackets the fit through every run; scale-down emulation, 4 threads, cold cache, 3 runs per size; "touched" = `--touch-buffers`, every arena block resident at its allocation):

| series | RssAnon peak per n (median [min-max], runs) | slope MB per doubling | step slopes |
|---|---|---|---|
| random DNA, base, plain, 10-100 M | 10 M: 37.9 [37.8-40.3]; 50 M: 48.7 [48.7-48.7]; 100 M: 50.8 [50.7-50.9] | **4.01** (3.76) | 4.64 / 2.08 |
| random DNA, base, plain, 10-200 M | + 200 M: 54.9 [54.9-55.0] | **3.88** (3.70) | 4.64 / 2.08 / 4.10 |
| random DNA, base, touched, 10-100 M / 10-200 M | 39.4 / 48.7 / 50.8 / 54.9 | 3.53 (3.27) / 3.53 (3.36) | 4.00 / 2.09 / 4.15 |
| random DNA, bounded (2 rounds), plain, 10-100 M | 10 M: 20.0; 50 M: 28.4; 100 M: 35.4 [34.4-35.4] | **4.44** (4.36) | 3.59 / 7.02 |
| random DNA, bounded (2 rounds), plain, 10-200 M | + 200 M: 38.9 [38.4-39.4] | **4.47** (4.44) | 3.59 / 7.02 / 3.53 |
| random DNA, bounded (2 rounds), touched, 10-100 M / 10-200 M | 22.4 / 28.4 / 35.3 / 37.4 | 3.65 (3.50) / 3.62 (3.59) | 2.59 / 6.86 / 2.15 |
| genome-like, base, plain / touched, 10-100 M | 38.3 / 48.7 / 50.8; 40.0 / 48.7 / 50.8 | 3.90 (3.81) / 3.33 (3.40) | 4.49 / 2.10; 3.73 / 2.12 |
| genome-like, bounded (2 rounds), plain / touched, 10-100 M | 20.0 / 30.6 / 34.4; 22.3 / 28.6 / 35.3 | 4.37 (4.40) / 3.70 (3.57) | 4.57 / 3.76; 2.71 / 6.72 |

All fits are below 5 MB per doubling (3.3-4.5), **but the margin of the bounded mode is thin (4.44-4.47) and the single step 50 -> 100 M of that mode is 7.0 MB per doubling** (the run-to-run spread of the managed arrays at 100 M is about 1 MB; the fit, which is what the gate asks for, is below 5). What grows with n in the emulation is the scaled part by construction: the arrays of R entries (base: 2 tasks x 8 B x R = 4.3 MB per doubling at 100 M, measured 4.1 for 100 -> 200 M; bounded: 4 tasks x 4 B x R = 4.3 MB, measured 3.5); the buffers that depend on the chunk count are constant (P = 375). A slope below 5 in this emulation therefore holds up to a doubling of n beyond which the scaled arrays alone exceed 5 MB (the 200 -> 400 M step would be 8.5 MB); at the real size those arrays are constant at their reference value (the real-size runs above), which is why the unscaled runs are the reference for memory.

**ESTIMATE of RssAnon at 3x10^9** (does not rest on the slope but on the phase model checked against the unscaled 500 Mbp runs): the worst phase is B2, the arena of 2 tasks (base) or 4 tasks (bounded), plus the runtime and threads (about 20 MB). Base mode: a task holds `2 x 32 MB + 1.05 MB + 376 x 32 KiB` = 77.4 MB, two tasks 154.7 MB (measured at 500 Mbp: 134.2 MB arena, 154.0 MB RssAnon, with P = 63 instead of 376); so **about 175 MB** (+ 313 more chunks x 32 KiB x 2 tasks = +20 MB over the 500 Mbp figure). Bounded mode: 4 x 33 MB + the round buffers (188 x 32 KiB) = 138 MB, **about 158 MB** (152.8-155.8 MB measured at 500 Mbp). Both below the cap by 80 MB or more. It rests on: the arena cap (160 MiB, enforced from inside, shown to fire), the independence of every other buffer of n, and the unscaled runs; it does not cover 40-bit positions (n >= 2^32) or a text of more than 4x10^9 symbols.

**G2: PASS.**

## 6. Gate G3 (ADR row E4): I/O, throughput and the abandon rule

**I/O** (explicit `pread`/`pwrite` of the SA and the tuple / PLCP files plus the bytes written through the mapped windows; the bounded mode's gather, 4 B/char of random reads through a mapping, is added; the cold read of the text, 0.25-0.4 B/char, comes on top and is in the block-layer figure):

| mode | scratch peak B/char | I/O B/char, U8 | U16 | U32 |
|---|---|---|---|---|
| base | **12.0** | **45.0** (genome-like 45.1: the overflow pairs) | 46.0 | 48.0 |
| bounded, 1 round | **8.0** | 33.0 | 34.0 | 36.0 |
| bounded, 2 rounds | 6.0 | 37.0 | 38.0 | 40.0 |
| bounded, 4 rounds | 5.0 | 45.0 | 46.0 | 48.0 |
| bounded, 8 rounds (100 Mbp, 3 runs: 8.41 s [8.19-8.71]) | 4.5 | 61.0 | | |

(the ADR table: 47 B/char with a 1.1 B/char LCP; the base mode is 45 + the LCP width - 1 + 8 x overflow share.) The ADR's gate is <= 50 B/char: **base mode 45.0-48.0, bounded mode 33-48 with up to 4 rounds** (more rounds exceed it: 8 rounds 61, and 8.4 s instead of 5.0 s at 100 Mbp).

**Throughput** (gate: >= 5 M symbols/s on 4 threads; PH4 only, from the SA file):

| input | mode | PH4 s | ns/symbol | Msym/s |
|---|---|---|---|---|
| dna_100m | base | 5.21 [4.65-5.27] | 52.1 [46.5-52.7] | **19.2** [19.0-21.5] |
| dna_100m | bounded, 2 rounds / 1 round | 4.98 [4.90-5.71] / 4.13 [4.01-4.21] | 49.8 / 41.3 | 20.1 / 24.2 |
| genome-like 100 M | base | 4.53 [4.24-5.49] | 45.3 [42.4-54.9] | **22.1** [18.2-23.6] |
| genome-like 100 M | bounded, 2 rounds / 1 round | 4.70 [4.53-4.86] / 3.60 [3.54-3.73] | 47.0 / 36.0 | 21.3 / 27.7 |

The slowest configuration at 100 Mbp (bounded mode, one thread) is 11.9 Msym/s; one thread in the base mode 13.4-13.6 Msym/s. At the unscaled real size the per-symbol time is higher in the base mode (the arrays of 32 MB miss the 2 MB L2): 500 Mbp random DNA, base 38.5 s [34.4-39.5] = 77 ns/symbol = **13.0 Msym/s**, bounded 2 rounds 27.3 s (54.6 ns, 18.3 Msym/s), 1 round 21.0 s (42 ns); a^n 500 Mbp: 19.6 s (39 ns, 25.5 Msym/s) base, 12.7 s bounded 2 rounds. All above 5 Msym/s.

**Abandon rule (LCP phase > 3x the SA phase), the SA phase re-measured in the same session** (E3 builder `e3`, K = n / 3000, 4 threads, cold cache, ranks from `mkranks`):

| input | SA phase s | PH4 base s (ratio) | PH4 bounded 2 rounds (ratio) | bounded 1 round (ratio) |
|---|---|---|---|---|
| dna_100m | 14.30 [14.13-14.41] (3 runs) | 5.21 (0.36) | 4.98 (0.35) | 4.13 (0.29) |
| genome-like 100 M | 13.78 [13.36-13.98] (3 runs) | 4.53 (0.33) | 4.70 (0.34) | 3.60 (0.26) |
| a^n 50 M | 18.65 [16.52-18.66] (3 runs) | 1.75 (0.09) | 1.44 (0.08) | 1.22 (0.07) |
| dna_10m / rep_10m | 1.71 / 2.10 | 0.59 (0.34) / 0.57 (0.27) | 0.71 (0.41) / 0.62 (0.30) | 0.87 (0.51) / 0.48 (0.23) |
| Fibonacci 10 M / tandem 171 x 10^4 / N run 1 M in 10 M / genome-like 10 M | 4.44 / 1.70 / 2.40 / 1.81 | 0.65 (0.15) / 0.22 (0.13) / 0.60 (0.25) / 0.73 (0.40) | 0.59 / 0.16 / 0.64 / 0.61 (all <= 0.34) | 0.60 / 0.16 / 0.64 / 0.60 (all <= 0.33) |
| 500 Mbp random DNA (real size) | 42.8 (1 run) | 38.5 (0.90) | 27.3 (0.64) | 21.0 (0.49) |
| a^n 500 Mbp (real size) | 172.9 (1 run) | 19.6 (0.11) | 12.7 (0.07) | 14.2 (0.08) |

Largest ratio 0.90 (500 Mbp random DNA, base mode); the rule's limit is 3. (The 500 Mbp SA phase was run once, with K = n / 3000 = 166,667; the PH4 time of the real-size runs is the one to use for the extrapolation.)

**G3: PASS** (I/O 45.0-48.0 B/char in the base mode, throughput 13-28 Msym/s, ratio <= 0.90).

## 7. G4 report (no gate)

Per phase, 4 threads, cold cache (B2 of the base mode uses 2 tasks; the scale-down runs have R = 266,667 at 100 Mbp):

| input | mode | B1 s | B2 s | B3 s | PH4 s | width | share of values >= 255 | compare steps per position | restart symbols (sum over the chunk starts) |
|---|---|---|---|---|---|---|---|---|---|
| dna_100m | base | 1.24 | 3.15 | 0.82 | 5.21 | U8 | 0 | 1.000 | 4,693 |
| dna_100m | bounded 2 rounds | 1.98 | 1.27 | 1.74 | 4.98 | U8 | 0 | 1.000 | 4,693 |
| genome-like 100 M | base | 0.91 | 2.90 | 0.79 | 4.53 | U8 | 1.39 % (>= 65,535: 0.50 %) | 1.000 | 856,432 |
| genome-like 100 M | bounded 2 rounds | 1.89 | 1.24 | 1.55 | 4.70 | U8 | 1.39 % | 1.000 | 856,432 |
| a^n 50 M | base | 0.48 | 0.65 | 0.60 | 1.75 | U32 | 99.9995 % (>= 65,535: 99.87 %) | 0.734 | **9,399,952,875** |
| a^n 50 M | bounded 2 rounds | 0.55 | 0.24 | 0.71 | 1.44 | U32 | 99.9995 % | 0.734 | 9,399,952,875 |
| rep_10m | base | 0.13 | 0.33 | 0.11 | 0.57 | U16 | 6.94 % | | 36,222 |
| Fibonacci 10 M | base | 0.16 | 0.38 | 0.14 | 0.65 | U32 | 99.995 % | | 958,122,772 |
| 171-bp tandem x 10^4 | base | 0.05 | 0.13 | 0.04 | 0.22 | U32 | 99.975 % | | 321,415,875 |
| N run 1 M in 10 M | base | 0.12 | 0.36 | 0.12 | 0.60 | U16 | 9.9975 % | | 19,001,188 |
| 500 Mbp random DNA, real size | base | 7.64 | 20.26 | 10.61 | 38.52 | U8 | 0 | 1.000 | 850 |
| 500 Mbp random DNA, real size | bounded 2 rounds | 10.71 | 7.01 | 8.63 | 27.30 | U8 | 0 | 1.000 | 850 |
| a^n 500 Mbp, real size | base | 5.56 | 6.65 | 7.42 | 19.63 | U32 | 100 % | | 1.57x10^10 |

* **Phases.** On random and genome-like DNA the base mode spends 60 % in B2 (the PLCP loop with its two scatters and the second distribution), 24 % in B1, 16 % in B3; the bounded mode shifts the time to the SA re-reads (B1, B3).
* **SIMD compare work on a^n.** The restart without a carry costs, per chunk, the length of the suffix at the chunk start: a^n 50 M with P = 375: 9.4x10^9 symbols summed over the chunk starts (on average n / 2 = 25 M symbols per chunk start; 188 per text position), 0.734 vector steps per position in total (36.7 M steps), B2 0.65 s: **negligible**. The ADR bound "restart cost <= max PLCP per chunk" holds (max PLCP 49,999,999, the largest start 49,999,999). ESTIMATE for a^n at 3x10^9 (not measured): `sum_c (n - c R)` = P n / 2 = 5.6x10^11 symbols = 4.4x10^9 vector steps at 128 two-bit symbols per 256-bit step (about 280 GB of sequential text reads), i.e. of the order of 5-30 s spread over the 2 B2 tasks; the carried length keeps every other position at one check. Measured at the real size: a^n 500 Mbp (P = 63, the text packed at 1 bit per symbol) 1.59x10^10 restart symbols, B2 6.65 s.
* **Escapes and width.** Random DNA has no value >= 255 in 100 Mbp (U8, 1.0 B/char). Genome-like text has 1.39 % and U8 with 11.1 MB of overflow pairs (0.111 B/char); texts with long repeats are U16 or U32 (above).
* **Peak scratch** 12.0 B/char (base), 8.0 / 6.0 / 5.0 B/char (bounded, 1 / 2 / 4 rounds); the ADR/6a limit is 8. **The base mode exceeds it; the bounded mode with one round is at it (8.000).**

### Both scratch modes (for the owner; the product choice is not made here)

| | base mode | bounded, 1 round | bounded, 2 rounds | bounded, 4 rounds |
|---|---|---|---|---|
| peak scratch, B/char (3x10^9 ESTIMATE) | 12.0 (36 GB) | 8.0 (24 GB) | 6.0 (18 GB) | 5.0 (15 GB) |
| I/O, B/char (U8 / U32) | 45 / 48 | 33 / 36 | 37 / 40 | 45 / 48 |
| PH4 at 100 Mbp, dna_100m / genome-like | 5.21 / 4.53 s | 4.13 / 3.60 s | 4.98 / 4.70 s | 6.29 / 6.37 s |
| PH4 at 500 Mbp unscaled, random DNA / a^n | 38.5 / 19.6 s | 21.0 / 14.2 s | 27.3 / 12.7 s | 33.4 / 14.3 s |
| RssAnon at 100 Mbp scaled / 500 Mbp unscaled | 50.7 / 154 MB | 33.7 / 153 MB | 34.6 / 153 MB | 31.4 / 156 MB |
| random access to a file | none (sequential) | PLCP file 4 B/char (400 MB at 100 M) | same | same |
| ADR gates (scratch <= 8, I/O <= 50) | scratch 12 > 8 | both met | both met | both met |

The base mode is purely sequential (every read and write is a stream); the bounded mode reads the PLCP file at random: **the bounded mode needs the PLCP file (4 B/char = 12 GB at 3x10^9, 3.6x the 1.1 B/char of the ADR) in the page cache for the gather** (the owner's Q4 accepts a dependence on the page cache; the floor of 6b is 1x(text + rank), this one is of another size). Measured (extra, report only; `run_e4_cgroup_sweep.sh`, cgroup-v1 limit = 160 MiB + rho x 4n, bounded 1 round, 100 Mbp random DNA, cold cache, one run each):

| rho (cache / PLCP file) | 3 | 1.5 | 1 | 0.75 | 0.5 |
|---|---|---|---|---|---|
| PH4 s (B3 s) | 3.84 (1.79) | 5.83 (2.42) | 5.59 (2.48) | 8.71 (3.54) | **killed after 240 s** |
| against the unlimited 4.13 s | 0.9x | 1.4x | 1.4x | 2.1x | > 58x |

The base mode as control under the same kind of limit (sequential, but scratch files of 12 B/char pass through the limited cache): 14.8 s at 160 MiB + 0.25 x 4n (2.8x) and 24.1 s at 160 MiB (4.6x), no cliff. So the bounded gather falls off a cliff between 0.75x and 0.5x of the PLCP file; the base mode degrades with the disk bandwidth only.

### Real-size runs, unscaled parameters (R = 8,000,000, P = 63, arena cap 160 MiB, 4 threads, cold cache, 3 runs where shown)

| text | n | mode | PH4 s | ns/symbol | RssAnon MB | arena peak MB | I/O B/char | scratch B/char | equal |
|---|---|---|---|---|---|---|---|---|---|
| random DNA | 5x10^8 | base | **38.5** [34.4-39.5] | 77.0 | **154.0** [153.9-154.0] | 134.2 | 45.00 | 12.000 | 3 of 3 kernel |
| random DNA | 5x10^8 | bounded 2 rounds | 27.3 [26.1-27.5] | 54.6 | 152.8 [152.8-153.5] | | 37.00 | 5.952 | 3 of 3 |
| random DNA | 5x10^8 | bounded 1 round / 4 rounds | 21.0 / 33.4 | 42.0 / 66.7 | 152.9 / 155.8 | | 33.00 / 45.00 | 8.000 / 4.928 | 1 of 1 each |
| random DNA | 5x10^8 | base, 1 / 2 threads | 52.4 / 38.2 | | 86.7 / 153.9 | | 45.00 | 12.000 | 1 of 1 each |
| a^n | 5x10^8 | base | 19.6 [19.3-19.9] | 39.3 | 153.8 [153.6-153.8] | | 48.00 | 12.000 | 3 of 3 closed form |
| a^n | 5x10^8 | bounded 2 rounds / 1 round / 4 rounds | 12.7 [11.8-13.1] / 14.2 / 14.3 | 25.5 / 28.5 / 28.6 | 154.8 / 153.3 / 155.3 | | 40.00 / 36.00 / 48.00 | 5.952 / 8.000 / 4.928 | 5 of 5 |

**ESTIMATE for 3x10^9 (labelled so; it rests on the 500 Mbp unscaled runs, linear scaling at fixed R, and files in the page cache):** base mode 77 ns/symbol x 3x10^9 = **about 230 s** [206-237] (random DNA; a^n about 120 s), bounded 1 round about 125 s, 2 rounds 165 s, 4 rounds 200 s; I/O 135 GB (base) or 99-135 GB (bounded) which at 3x10^9 no longer fits the page cache of a 16 GB machine (SA 12 GB, tuples 36 GB): the time is then `max(CPU, I/O / disk bandwidth)`, e.g. 270 s at a sustained 500 MB/s; not measured. The base-mode PH4 time is 0.5-0.9x the SA phase of the same size (SA phase 14.3 s per 100 Mbp here, 430 s linear to 3x10^9; 42.8 s at 500 Mbp, 257 s linear to 3x10^9), so far from the abandon limit of 3x.

## 8. What was not verified, and the limits of the evidence

* **No page-cache limit on the base mode's streams below the sweep above** and none on the text read; every run starts cold and then has all of the 15 GB for the page cache (the 500 Mbp runs write 6 GB of scratch and keep it cached: the disk of this VM is a virtio disk, probably host-cached, real NVMe or network latency is not represented).
* **Sizes**: the gates are at 100 Mbp; the unscaled runs reach 500 Mbp (random DNA and a^n); nothing above (the disk allows no more: the base mode needs 12 B/char of scratch + 4 B/char SA + the output), no genome-like text above 100 Mbp (the 3x10^9 numbers are estimates), no both-strand 6.2x10^9 text, no 40-bit positions (the prototype's positions are `uint`).
* **The oracle is `LcpKernel`** (the in-repository kernel, pinned to `LcpReference` by the earlier work and, in the exhaustive and differential tests, by direct comparison with the reference), not an independent implementation; for a^n the closed form is independent. It needs about 1 GB at 100 Mbp and about 5 GB at 500 Mbp (not measured).
* **Threads**: 1, 2 and 4 only (4 vCPU). In the base mode B2 runs 2 tasks (the arena cap), B1 and B3 one thread: the 4 threads are not all used; 2 and 4 threads give the same time (4.7 s and 5.2 s at 100 Mbp, the difference is noise).
* **No resume, no crash safety**, no consumer for PH5 (the hook is tested with a recording consumer only), no `fsync` of the outputs in the timing (reported separately in the rows: 0.4 ms at 10 M).
* **The time figures of the scale-down runs** are those of a build whose chunks are 30 times smaller than the real ones (the arrays fit the L2 cache); the unscaled runs are the reference for time (the base mode is 1.5x slower per symbol there).
* The bounded mode was measured under a cache limit only at 100 Mbp, one run per point, and the hang point was found by a 240 s kill, not located more finely.
* **The cgroup sweep and the real-size runs use the same SA files; the SA phase at 500 Mbp was run once.**

## 9. Open points for E5 (PH5) and for 1b.7

* **Interface for PH5.** `ILcpChunkConsumer.Consume(long firstRank, ReadOnlySpan<uint> lcp)` is called by B3 in strictly ascending rank order with the resolved values (never escapes), once per chunk of R ranks, before the chunk is written; PH5 (child table, external stack) can read the stream directly and needs no second pass over the LCP file. PH5 also needs `SA[i]` (NEXTSYM `text[SA[i] + LCP[i]]`): B3 of the base mode does not read the SA; the bounded mode's B3' does (it has the SA block of the chunk in memory at the time of the gather: handing it over would be 5 lines), the base mode would have to read the SA file again (4 B/char).
* **The mode for the product** (the owner's decision): the base mode meets the I/O gate but not the proposed 8 B/char scratch limit (12.0); the bounded mode with one round is at the limit (8.0) with 33 B/char of I/O but needs a 4 B/char PLCP file in the page cache; with 2 rounds 6.0 B/char and 37 B/char. Possible middle: a compact PLCP file with an escape structure (U8 + a sorted escape list searched with a moving cursor per rank chunk) to bring the cache need back to the ADR's 1.1 B/char; not built.
* **B2 uses 2 of 4 threads in the base mode** because two Phi + rank chunks of R = 8 M fill the 160 MiB arena. Options for 1b.7: R = 4 M (P = 750, four tasks), a rank array of 3 bytes (not possible for n > 2^24 ranks per chunk), or B2 with one array and the rank implicit (the bounded mode).
* **B1 and B3 are single-threaded** (1.2 s and 0.8 s per 100 Mbp, 7.6 s and 10.6 s at 500 Mbp): fine for the gates; a parallel B1 over rank ranges needs one buffer set per thread.
* **Container writer (1b.7c).** The two output files are the LCP and LCP_OVERFLOW section bodies with `COUNT = n + 1` / `PARAM0` = overflow entries; the base file is written at its final length through mapped windows, so the writer can place it at a section offset directly (the section table is written last); the CRC-32C of the sections needs one more read pass or a hash in the consumer hook.
* **Overflow pairs with `POS_WIDTH` 8** (n >= 2^32): the tuples would be 16 and 12 bytes; not in the prototype.

## 10. Mutation audit of the new code

Hand-made single-site mutants of the new production code: `LcpGeometry` (7), `Spool` (11), `PlcpScan` (9), `LcpOutput` (15), `LcpBuild` (33), `Experiments/E4` (11) = 86 in round 1, plus 12 added in round 2 (replacements of the mutants that did not compile and mutants of `E4.Healthy`): lists `analysis/mutation/e4_*_mutants.json` (generators `make_e4_mutants.py`, `make_e4_round2_mutants.py`, which also check that every pattern is unique), runner `run_e4_mutants.py` with a **watchdog** (build limit 300 s, test limit 240 s; the whole process group, test host included, is killed on expiry and the mutant counts as killed(hang)), drivers `run_e4_audit.sh`, `run_e4_audit_round2.sh`, `run_e4_audit_round3.sh` (run in a separate git worktree, so the checkout is never mutated; removed afterwards), raw results `results_e4_*.tsv`, merged table with verdicts and notes `results_e4.tsv` (`merge_e4.py`). Stage 1 runs the Ph4 tests without the slow exhaustive class; a survivor is run again against the whole Ph4 namespace.

| | mutants |
|---|---|
| listed (round 1) | 86 |
| did not compile (analyzer errors S4487, S1066, CS0649, S1854, S1481 twice), replaced by B12b, P09b, B26b, E07b, E08b, E09b in round 2 | 6 |
| **compiled in round 1** | **80** |
| killed in round 1 | 75, of them 1 by the watchdog (S05: a reader whose position never advances) and 1 by a test-host crash (S07: a buffer flushed one tuple late writes into its neighbour) |
| survived round 1 | 5: B06, B29, E05, E06, P06 |
| round 2 (14 mutants, head `1b07e175`): the survivors B06 / E06, the replacements, the mutants of `E4.Healthy` | 11 killed (E06, B12b, B26b, E07b, E09b, P09b, H01-H05), 3 survived (B06, B06b - the same mutation in the bounded worker - and E08b) |
| round 3 (3 mutants, head `94daa207`), after two strengthened tests | B06, B06b, E08b killed |
| **final: listed 98, compiled 92** | **killed 89**, closed by a refactoring and a test 1 (E05), equivalent argued 2 (P06, B29), **survivors 0**, not compiled 6 (each replaced by a mutant that compiled) |

The five survivors of round 1: **B06** (B2 accepts a chunk with too few tuples): the first test of non-permutations (one tuple moved between two chunks) is also caught at the chunk with one too many; the one-task version added after round 1 still passed in round 2, because the long chunk raises the same message (a build that skipped the short chunk ends at the long one); with the numbers in the message asserted ("holds 499 Phi tuples, expected 500" when the short chunk comes first) it is killed in round 3, as is B06b, the same mutation in the bounded worker. **E05** (the exit code ignored a scratch entry left behind): not testable inside the long command method, so the condition became `E4.Healthy` and a test pins each term; its mutations H01-H05 are killed. **E06** (the closing comparison of the structure scan of `check-lcp`): every corruption in the middle of the overflow list is also found by the pair that follows it; a test with the last escape slot overwritten found it (round 2); a pair value below the escape (E09b) was killed in round 2 and a pair ahead of its slot (E08b) in round 3, by a test on the last pair. **P06** and **B29** are equivalent:
* **P06** (the PLCP loop does not set `l = 0` at a position whose predecessor is the empty suffix): Kasai's invariant `PLCP[j] >= PLCP[j-1] - 1` with `PLCP[j] = 0` gives `PLCP[j-1] <= 1`, so the carried length `max(PLCP[j-1] - 1, 0)` is already 0; at the first position of a chunk it starts at 0 anyway. The assignment is a defensive reset that no valid suffix array can make necessary.
* **B29** (a failed B2 task does not stop the others): the others finish their chunks, the first error is stored and rethrown after all tasks ended, the spools remove the scratch files: only the time to the exception changes.

Not listed as mutants because they cannot change a result: the prefetch distances and addresses (16 positions / ranks ahead), the order of `Consume` and `Write` in `Deliver`, the slice size of the output conversion, the choice of `Math.Clamp` bound 0 or 1 in the task planning (a later `Math.Max(1, ...)`).

Limits of the audit: mutants are hand-made and single-site; the reporting rows of `e4` are checked by a few assertions only (`parameters`, `n`, the verification row, the exit code); the older E0-E3b code was not audited; the first round's stage-1 tests are the Ph4 tests without the exhaustive class (a mutant that survived them was run again against the whole namespace: B06, B29, E05, E06, P06 survived both stages); the tests of the commands run on small families (the commands' timing code is not mutated).
