# E5: PH5 (the lcp-interval scan with an external stack: child table in layouts A and F, NODE_COUNT, LRS) - measured results

Work package of the owner decision of 2026-11-11 (ADR-004 section 6c authorises E4-E9; the yes for E5 was given on 2026-11-13). ADR-004 section 3.2 designs the scan, section 4.2 row PH5 and section 5 row E5 define the prototype ("PH5 with the external stack on a^n (50 M), N run 30 M, Fibonacci; both child layouts", gate "child table equal to the reference; spill volume equals 12-16 B x (max LCP + 1); mapped dirty pages <= 64 MB sampled from `/proc/meminfo` Dirty; write amplification < 2", abandon "layout A thrashes: use layout F"). The gates G1-G4 below were fixed by the orchestrator before the work.
Code: `apps/SuffixTree.EmProto/Ph5` (`SpillStack.cs` the external stack; `AlignedWindowWriter.cs` the bounded mapped windows; `CompactLcpFile.cs` the E4 output form read forwards and backwards; `ChildScan.cs` the forward scan, the sinks of layout A (both codecs) and of the PSV half of layout F; `NsvScan.cs` the backward half of layout F; `ExternalChild.cs` the entry point on a given LCP file; `NextSymScan.cs` the optional NEXTSYM pass), `Harness/DirtySampler.cs` (dirty pages every 20 ms) and `Experiments/E5.cs` (commands `e5`, `e5-oracle`); tests `tests/SuffixTree/SuffixTree.EmProto.Tests/Ph5`.
Raw rows: `e5_*.tsv` in this directory (`R <TAB> E5 <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes`; the rows of the independent sampler start with `ext_`, `device_write_bytes` is taken by the run script), scripts `run_e5_setup.sh`, `run_e5_protocol.sh`, `sample_proc.py` (the independent /proc sampler), `summarize_e5.py`.
Everything below is **measured here** (4 vCPU Intel Xeon, 15 GB, kernel 6.18.44, .NET 10.0.14, Release, workstation non-concurrent GC, nothing else running, `DOTNET_GCHeapHardLimit` = 96 MiB, every timed run in a fresh process after `sync; echo 3 > /proc/sys/vm/drop_caches`, i.e. a **cold page cache**) unless it is marked ESTIMATE. PH5 is single-threaded (the optional concurrent NSV pass uses a second thread). **The block sizes are the real ones** (256 K entries per stack block, 1 MiB windows): no scale-down is needed for PH5, whose buffers do not depend on n.
The binary under test of each group is the build of the commit in its `head` row: `10de29a9` for `adv`, `shm`, `nosp`, `fused`, `slope`, `par`, `slots` and `ph4ref`, `5a1c7660` for `nextsym` and `a200`, `6fe1b628` for `stats`; the later commits add the NEXTSYM pass (`NextSymScan`, `--nextsym`, `.nsym` of the oracle), the `BackwardConcurrent` field of `ChildBuildStats`, the usage text and tests, nothing on the path of the measured groups before them.

## 1. What PH5 does here

Input: the compact LCP of PH4 (`<prefix>.lcp` base array of 1, 2 or 4 bytes per entry, `<prefix>.ovf` pairs `(rank, value)`; the output form of E4), or the dense resolved values of the PH4-B3 hook (`ILcpChunkConsumer.Consume(firstRank, values)`), which is how the scan runs fused with PH4. The text and the SA are not needed (NEXTSYM, report only, needs both).

**The forward scan** (`ChildScan`, one pass in rank order): the bottom-up scan of the lcp-interval tree of the in-RAM `LcpIntervalScan` with the stack in a `SpillStack<IntervalEntry>`. An entry is `(l, lb, first, last)`, **16 bytes**: `first` is needed when the interval closes (`up[rb]` or `down[lb]` = first l-index), `last` when the next l-index arrives (`next[last] = k`), so the ADR's 12 B `(l, lb, lastChildIndex)` is not enough for layout A; the 12-16 B of the ADR are 16 here. The scan reports three kinds of events to the sinks: `LIndex(previous, k)`, `Close(l, lb, rb, first, lastChild)` (the events of the in-RAM `IIntervalSink`) and `Position(i, psv)` (the previous smaller value of every rank, which the scan knows at no cost: the open interval's left boundary). The statistics (internal nodes `I`, `NODE_COUNT = 1 + I + n + 1`, `LRS_DEPTH`, `LRS_RANK` = left boundary of the first-closing deepest interval) come from the `Close` events.
* **Layout A** (`TaggedChildSink`, `PackedChildSink`): the writes of the in-RAM `TaggedChild` / `PackedChild` (next link at the previous l-index, up at `rb`, down at `lb`), into a mapped file of final length `4 (n + 1)` through the `AlignedWindowWriter` (8 windows of 1 MiB, least recently used released first, **flushed (msync) when released**). Both codecs are implemented and verified; which one the container takes is not decided here (the tagged codec holds ranks in 30 bits: texts up to 2^30 - 1 symbols).
* **Layout F, forward half** (`PsvSink`): `PSV[i]` (previous smaller value, strict, the definition of `LcpPyramid.PreviousSmaller`; 0xFFFFFFFF = none) written at entry `i`, strictly sequential.
* **Layout F, backward half** (`NsvScan`): `NSV[i]` (next smaller value, strict; `n + 1` = none) written in **descending** rank order by a scan of the LCP file read backwards (the overflow pairs are read by a cursor that moves down the pair file); the stack holds the distinct values seen so far (strictly increasing from the bottom, 8 B entries `(value, rank)`), in a `SpillStack<ValueEntry>`. Every write is sequential, forwards (PSV) or backwards (NSV). Layout F is therefore 8 B/char (two arrays) against 4 B/char of layout A, and costs a second pass over the LCP (reading 1.0-4.0 B/char, depending on the width).

**The external stack** (`SpillStack<T>`, ADR-004 section 3.2): blocks of B = 256 K entries (4 MiB at 16 B); **two block buffers are in RAM** (the block of the top entry and either the full block below it or an empty spare). When the top block is full and an entry is pushed, the block below is written to the spill file at its own offset `j x B x sizeof(T)` (the file is accessed only at its end: LIFO) and its buffer becomes the new top block; when a pop empties the top block the block below becomes the top, **read back** from the file when it is not resident. After a refill the block above is only a spare and the block below is not needed until the top has climbed through a whole block again, so a push / pop alternation at a block boundary costs no I/O (tested); a block that was read back and not modified (its file copy is current) is evicted without a write (`CleanEvictions`). The spill file is `fsync`ed after every 16 MB written (`SpillFlushBytes`), so that its dirty pages are bounded too (section 5: without it a^n 50 M leaves 801 MB dirty).

## 2. Design choices the ADR does not fix (reason, alternative)

| choice | reason | alternative not taken |
|---|---|---|
| Entry of 16 B `(l, lb, first, last)` | `Close` writes `first` into `up` / `down`, a next link needs `last`; nothing else is known when an interval closes | 12 B `(l, lb, lastChildIndex)` as in the ADR text: `first` would have to be read back from the child table (a random read) |
| Two block buffers, spare above, block below resident, clean blocks not rewritten | ADR: "at least one full block hot"; the hysteresis makes the amortised I/O O(1/B) per operation (tested: at most one write and one read in 2,000 alternations at every boundary) | a third buffer and an I/O thread (read-ahead, write-behind): would hide the spill latency (section 6: the spill costs a^n 50 M 1.1 s of 2.4 s) |
| `AlignedWindowWriter`: aligned windows, LRU, flush on release | layout A writes behind the cursor at `previous`, `lb`, `rb`: a forward cursor window, a backward sweep and a few older windows; the `MappedWindowWriter` of E3 / E4 (one sliding window, no flush) cannot jump back and forth and leaves the dirty pages to the kernel (system-wide Dirty = the whole file, 400 MB at 100 Mbp) | no flush (dirty pages unbounded: 1.5 GB is the default `dirty_background_bytes` of this VM at 10 % of 15 GB); `RandomAccess.Write` of a window buffer (a pwrite per window: the owner's rule says mapped files) |
| The spill file is written with `pwrite` and `fsync`ed every 16 MB | bounds Dirty to the 64 MB of the gate | no fsync (faster by 0.4-0.6 s on a^n 50 M, 801 MB dirty); `O_DIRECT` (not in the BCL) |
| Both codecs of layout A as sink structs that make the writes of the in-RAM `TaggedChild` / `PackedChild` | "do not decide the container codec": the output is verified by decoded content, and entry by entry against the raw arrays of both | one codec |
| PSV from the same scan, NSV from a separate backward pass | the forward scan knows the PSV for free (`lb` of the open interval); the NSV needs the future, so the ADR's backward scan | both from the backward pass (needs the PSV stack as well); NSV from the forward scan by writes behind the cursor (that is layout A again) |
| `NSV` and `PSV` strict, stored as `uint`: none = 0xFFFFFFFF on the left, `n + 1` on the right | the definitions of `LcpPyramid.PreviousSmaller` / `NextSmaller` (the in-RAM oracle) and of REVIEW-S1-0 MIN-1 | `PSV + 1` with 0 = none |
| The scan is a consumer of the PH4 hook and is created at the first chunk | the stack blocks (8 MB) and the windows are not allocated while PH4 is in its B2 peak (154 MB at 500 Mbp unscaled); the NSV pass needs the whole LCP and runs after PH4 from the LCP file PH4 wrote | creating the scan before PH4 starts |
| Child enumeration through layout F exists only in the test decoders (`NsvView`) | the prototype writes the structure, the reader belongs to 1b.x; the decoder proves that the structure is sufficient: first l-index = first position of the minimum of `LCP(lb..rb]`, reached by the chain of next smaller values from `lb + 1`; next l-index = the same chain from `k + 1` when its minimum equals `LCP[k]` | a production reader |
| Positions are `uint` (n < 2^32 - 2), `NotSupportedException` above | as E3b / E4 | 40-bit positions |

## 3. Verification

| item | result |
|---|---|
| EmProto tests before the work | 180, all green |
| tests first, against the stub (`823e5d63`) | **39 of 219 red, all 39 by `NotImplementedException`** (every test of the new classes needs one of the stubbed members; the tests that expect an exception type fail with `NotImplementedException` instead), the 180 older green; the explicit fuzz class is not in the count |
| after the implementation (`bc567333`) | 36 of the 39 green on the first run; the other three had wrong expectations of their own, and the production code was not changed for them: the width-coverage test (a^70,000 is U16, not U32: a^150,000 is), the spill-volume arithmetic at n = B (nothing spills, the bound `B + 1 .. 2 B` entries in RAM holds only above one block) and a Fibonacci text whose stack (39 entries) is too shallow to spill at B = 16. One test of the red commit had a sample without an overflow list (Fibonacci n = 3,000 is U16): corrected before the commit, so that it failed by `NotImplementedException` like the others |
| tests written after the code | labelled in `Ph5/ChildComponentTests.cs` (and the tests added for the mutants): configurations, limits, the concurrent NSV pass, the backward spill, the NEXTSYM pass, the file reader, the smaps parser, the sampler, the commands |
| after the audit | **255 default tests (3 m 34 s) and 255 in the checked variant** (`-p:EsaBounds=true`, `ESA_EXPECT_BOUNDS=true`, 3 m 50 s): 180 -> 255; the 180 older tests are unchanged except for the `global::` qualification of three names in three files (the test project now references the Ukkonen tree, whose type `SuffixTree.SuffixTree` shadows the namespace name `SuffixTree` in a fully qualified name) |
| exhaustive (default run, tiny blocks 2-5, windows 16-64 bytes, 1-3 slots, compared with the in-RAM structures) | every binary string of length 1-12 (8,190) and every ternary string of length 1-6 (1,092), **each built twice** (tagged codec with PSV and NSV, packed codec): **18,564 PH5 builds**, all equal entry by entry (CHILD against the raw array of the in-RAM tagged / packed layout, against the three plain arrays `next / up / down` of the reference layout; PSV / NSV against the `LcpPyramid`; NODE_COUNT, internal nodes, LRS_DEPTH, LRS_RANK against `EsaStatistics`) |
| child enumeration and the Ukkonen oracle (tiny texts) | binary strings of length 1-8 (510) and ternary of length 1-5 (363) plus seven named strings (880 texts): the event stream of the walk through layout A (tagged, packed) and through layout F equals the stream of `EsaIndex.Traverse`, and the logical hash (`SuffixTreeSerializer.CalculateLogicalHash`) of the three walks equals the hash of the Ukkonen tree |
| differential, adversarial families (`StringFamilies.Adversarial`) | n = 200, 1,000, 5,000 with two seeds (more than 150 cases), block 2..50 entries, window 16 B..1 KiB, 1..4 slots, both codecs, F with the tagged codec; the families at 10^5 (all, both codecs) and a^n, Fibonacci, Thue-Morse at 10^6: all equal |
| fuzz (`ChildScanFuzzTests`, explicit) | seed 1, 120 s: **317 cases** (164 with a spill, 38 with the concurrent NSV pass, 180 with the packed codec, deepest stack 100,000 entries); seed 2, 60 s: **183 cases** (96 / 21 / 112), all equal to the in-RAM structures; texts 10^3-10^5 symbols, blocks 2-50, windows 16 B-4 KiB, 1-4 slots |
| spill arithmetic | a^n at four (n, B) pairs: height = n = max LCP + 1, blocks written = `max(0, floor((n - 1) / B) - 1)`, bytes = blocks x B x 16, every block read back once, the entries left in RAM between B + 1 and 2 B; the backward stack on `a^k b` (8 B entries) the same formula; `SpillStack` against `Stack<T>` for every B from 2 to 50, 4,000 random operations each; hysteresis: at most one block written and one read in 2,000 push / pop alternations at each of four boundaries; a clean block is not written again |
| hook | the scan fed from random chunks of 1-40 ranks equals the file mode; the real `ExternalLcp.Build(..., consumer: scan)` (two threads, R = 61) equals it too, then the NSV pass over the LCP file PH4 wrote |
| real inputs (section 4) | 36 verified protocol runs on 9 inputs, 6 verified fused runs, 5 NEXTSYM runs, a^n 200 M against the closed form |

## 4. Gate G1: identity with the in-RAM structure

Every verified run compares the files PH5 wrote byte for byte with the files of `e5-oracle`, which come from the in-RAM structures of the product: `ChildTable.Build` of the tagged and of the packed layout (the raw field arrays, read through reflection), the `LcpPyramid` (PSV and NSV of every rank, fanout 32) and `EsaStatistics`, all built from the same LCP file (the `LcpKernel` file of `lcp-oracle`; its equality with PH4's output is the gate G1 of E4). The program's own `--verify` does the comparison, cold cache, fresh process; each verified run is the first of the runs of its configuration. The reference files exist only during the runs of their input (4 x 4 B/char).

| input | n | LCP width | NODE_COUNT | LRS_DEPTH | stack peak | A tagged | A packed | F: PSV, NSV | A + PSV + NSV in one scan | fused with PH4: A, F | NEXTSYM |
|---|---|---|---|---|---|---|---|---|---|---|---|
| a^n 10 M | 10 M | U32 | 20,000,001 | 9,999,999 | 10,000,000 | equal | equal | equal | equal | - | - |
| **a^n 50 M** | 50 M | U32 | 100,000,001 | 49,999,999 | 50,000,000 | equal | equal | equal | equal | equal, equal | equal |
| **N run 30 M** (run of 20 M) | 30 M | U32 | 56,230,641 | 19,999,999 | 20,000,000 | equal | equal | equal | equal | - | equal |
| Fibonacci 10 M | 10 M | U32 | 19,999,997 | 5,702,885 | 38 | equal | equal | equal | equal | - | - |
| **Fibonacci 30 M** | 30 M | U32 | 59,999,994 | 15,069,648 | 39 | equal | equal | equal | equal | - | - |
| rep_10m | 10 M | U16 | 18,568,249 | 790 | 64 | equal | equal | equal | equal | - | - |
| dna_10m | 10 M | U8 | 16,229,747 | 25 | 15 | equal | equal | equal | equal | - | equal |
| **dna_100m** | 100 M | U8 | 162,149,501 | 25 | 18 | equal | equal | equal | equal | equal, equal | equal |
| **genome-like 100 M** | 100 M | U8 + 1,390,028 overflow pairs | 166,050,648 | 499,999 | 500,000 | equal | equal | equal | equal | equal, equal | equal |
| a^n 200 M (no in-RAM reference: the stack of the in-RAM scan alone grows four int arrays to 2^28 entries) | 200 M | U32 | - | - | 200,000,000 | closed form: equal | - | - | - | - | - |

**36 of 36 verified protocol runs equal** (9 inputs x A tagged, F, both, A packed; every one also compares NODE_COUNT, the internal-node count, LRS_DEPTH and LRS_RANK), **6 of 6 fused runs** (PH4 base mode with the scan as the consumer of its B3 hook, scale-down of PH4 by n / 3x10^9, then the NSV pass over the LCP file PH4 wrote; dna_100m, genome-like 100 M, a^n 50 M, layouts A and F), **5 of 5 NEXTSYM runs**; the a^n 200 M child table equals its closed form (`check_child_allsame.py`: `Down | (k + 1)` at every rank `k < n`, 0 at rank `n`). No run differed. The statistics are the in-RAM ones in all 47 runs. (The remaining runs of the groups `adv`, `shm`, `slots`, `par`, `slope`, `a200` were not compared: they exercise the same code with other block, window, slot and thread parameters, whose independence is what the differential tests and the fuzz check.)

Tiny texts and the Ukkonen oracle, the exhaustive builds, the families and the fuzz are in section 3. **G1: PASS.**

## 5. Gate G2: RssAnon <= 256 MB including the runtime, 100 Mbp and a^n 50 M, real block sizes

PH5's buffers do not depend on n: two stack blocks of 256 K entries (8.39 MB at 16 B, 4.19 MB for the 8 B stack of the NSV pass), the writer windows (8 x 1 MiB for CHILD, 2 x 1 MiB for PSV and for NSV, mapped file pages, not anonymous), the LCP read buffer (1 MiB, 4 MB as `uint[]`) and the runtime. **Nothing is scaled down.** Each number is the maximum of the program's own sampler (`RssAnon` of `/proc/self/status` every 10 ms) and, in the second column, of the independent sampler (`sample_proc.py`, `/proc/<pid>/status` every 20 ms, a separate process); medians of 3 runs, the arena column is the peak of the native buffers charged to the `ArenaBudget`.

| input (3 runs each) | layout A | layout F | A + PSV + NSV | independent sampler (A / F / both) | arena peak |
|---|---|---|---|---|---|
| **dna_100m** | **14.3 MB** [14.3-14.3] | **19.4 MB** [19.4-19.4] | 19.5 MB [19.4-19.8] | 14.4 / 19.4 / 19.5 | 8.4 MB (12.6 MB with the concurrent NSV pass) |
| genome-like 100 M | 23.4 MB [23.1-23.4] | 23.1 MB [23.0-23.1] | 23.1 MB [23.1-23.1] | 23.4 / 23.1 / 23.1 | 8.4 MB |
| **a^n 50 M** (stack 800 MB, 793 MB spilled) | **19.5 MB** [19.5-19.5] | 19.2 MB [19.2-19.2] | 19.3 MB [19.3-19.6] | 19.5 / 19.2 / 19.3 | 8.4 MB |
| N run 30 M (stack 320 MB) | 19.5 MB | 19.2 MB | 19.6 MB | 19.5 / 19.2 / 19.6 | 8.4 MB |
| Fibonacci 30 M | 11.1 MB | 13.0 MB | 13.3 MB | 11.1 / 13.4 / 13.3 | 8.4 MB |
| rep_10m / dna_10m / a^n 10 M | 11.8 / 13.9 / 19.4 MB | 15.0 / 19.2 / 19.1 MB | 15.4 / 19.3 / 19.1 MB | | 8.4 MB |
| **a^n 200 M** (stack 3.2 GB, 3.19 GB spilled; layout A, no reference) | **19.6 MB** [19.6-19.6] | | | 19.6 | 8.4 MB |
| fused with PH4 (base mode, scale-down by n / 3x10^9): dna_100m / genome-like 100 M / a^n 50 M | 51.0 / 51.0 / 38.2 MB | 51.0 / 51.0 / 38.3 MB | | 51.0 / 51.0 / 38.2 | PH4's cap of 160 MiB shared with nothing: PH5 has its own |

The worst of the 217 runs of PH5 alone (files, layouts A / F / both, tmpfs spill, slots, concurrent NSV, slope, NEXTSYM, a^n 200 M) is **32.7 MB** (genome-like 100 M, layout F with the concurrent NSV pass; both samplers read the same value); the fused runs are the RssAnon of PH4 alone (50.8 MB at the same scale, E4): the scan is created at the first chunk of the hook, after B2 has released its arenas, so its 8 MB and windows are not on top of B2's peak.

**Slope of RssAnon against n** (least-squares fit of the median RssAnon peak over log2 n; layout A + PSV + NSV, 3 runs per size; "touched" = `--touch-buffers`, every arena block resident at its allocation):

| series | RssAnon peak per n (median) | slope MB per doubling |
|---|---|---|
| a^n 10 / 25 / 50 M, plain | 19.14 / 19.20 / 19.26 | **0.05** (independent sampler 0.05) |
| a^n 10 / 25 / 50 M, touched | 19.11 / 19.55 / 19.24 | 0.07 (0.06) |
| a^n 10 / 50 / 200 M, layout A only (adv and a200) | 19.43 / 19.53 / 19.62 | 0.05 |
| random DNA 10 / 50 / 100 M, plain | 19.26 / 19.65 / 19.48 | **0.08** (0.08) |
| random DNA 10 / 50 / 100 M, touched | 23.73 / 23.79 / 23.90 | 0.05 (0.05) |

All fits are below 0.1 MB per doubling (gate: 5), also with all runs in the fit (0.05-0.07) and with a^n 200 M (the stack there is 3.2 GB). **ESTIMATE for 3x10^9** (it does not rest on the slope but on the phase model: nothing in PH5 depends on n): 11-24 MB (PH5 alone, up to 33 MB with the concurrent NSV pass), about 175 MB fused with PH4 in the base mode (PH4's estimate of E4, whose B2 peak dominates; PH5 allocates 8-13 MB in B3 after it); not measured above 200 M (a^n) and 100 M (random and genome-like) for PH5 and not at all for 40-bit positions.

**G2: PASS.**

## 6. Gate G3 (row E5): spill volume, dirty pages, write amplification, time against PH4

### 6.1 Spill volume against 12-16 B x (max LCP + 1)

The stack is bounded by `1 + max LCP` entries (Claim K1 of the ADR), and the bound is attained when the lcp-interval tree contains a chain as deep as the longest LCP (a^n, a long run); the stack height is the depth of the lcp-interval tree, so a text whose longest repeat is not a chain (Fibonacci: max LCP 15 M, depth 39) stays small. Spill bytes = `16 x B x max(0, floor((height - 1) / B) - 1)`: whole blocks of the entries that do not fit the two resident blocks. The model is pinned by tests and equals the measurement exactly.

| input | n | max LCP | peak height | spill written (blocks of 4 MiB) | 16 B x (max LCP + 1) | **spill / bound** | read back | peak spill file |
|---|---|---|---|---|---|---|---|---|
| a^n 10 M | 10 M | 9,999,999 | 10,000,000 | 155,189,248 (37) | 160,000,000 | **0.9699** | 155,189,248 | 155,189,248 |
| **a^n 50 M** | 50 M | 49,999,999 | 50,000,000 | 792,723,456 (189) | 800,000,000 | **0.9909** | 792,723,456 | 792,723,456 |
| a^n 200 M | 200 M | 199,999,999 | 200,000,000 | 3,191,865,344 (761) | 3,200,000,000 | **0.9975** | 3,191,865,344 | 3,191,865,344 |
| **N run 30 M** (run of 20 M) | 30 M | 19,999,999 | 20,000,000 | 314,572,800 (75) | 320,000,000 | **0.9830** | 314,572,800 | 314,572,800 |
| **Fibonacci 30 M** / 10 M | 30 M | 15,069,648 / 5,702,885 | 39 / 38 | 0 | 241,114,384 / 91,246,176 | 0 (the chain is not there) | | |
| genome-like 100 M | 100 M | 499,999 | 500,000 | 0 (500,000 < 2 x 262,144 entries in RAM) | 8,000,000 | 0 | | |
| rep_10m, dna_10m, dna_100m | | 790, 25, 25 | 64, 15, 18 | 0 | 12,656 / 416 / 416 | 0 | | |

"Up to block granularity": the difference to the bound is the entries that stay in RAM at the peak, between `B + 1` and `2 B` = 262,145 to 524,288 entries (4.2-8.4 MB): 1.15 blocks (a^n 10 M), 1.73 (a^n 50 M), 1.28 (N run). The ratio per input is in the table; the 3x10^9 consequence (ESTIMATE): a^n or a run of 3x10^9 would spill 16 x 3x10^9 = 48 GB (the ADR's 36-48 GB), a text with its longest chain at the 22 M of the hg19-era Genome file 352 MB, random DNA and every text whose lcp-interval tree has fewer than 524,288 nested intervals nothing. The NSV stack of layout F holds 8 B entries and spills only when the LCP *decreases* to the left for more than 524,288 ranks (the text `a^k b`, tested: 8 B x blocks); on a^n it has a peak height of 1.

### 6.2 Dirty pages (system-wide Dirty, baseline subtracted, every 20 ms) and the per-process view

Two independent samplers: the program's (`/proc/meminfo` Dirty every 20 ms; baseline = median of five samples before the start, subtracted) and `sample_proc.py` (a separate process, 20 ms, own baseline); the per-process view is the sum of `Private_Dirty + Shared_Dirty` of the mappings of `/proc/self/smaps` whose path is an output file of the run (`smaps_rollup` would add anonymous memory). Windows are flushed (msync) when released, the spill file is `fsync`ed every 16 MB. Medians of 3 runs, nothing else running (the baselines were 0-0.4 MB).

| input | layout | Dirty peak: own / independent | mapped dirty, per process | WA 1 (pages dirtied) | device bytes | output bytes | WA 2, spill on disk: (device - spill) / output | WA 2, spill on tmpfs: device / output |
|---|---|---|---|---|---|---|---|---|
| a50 | A | 16.8 / 14.8 MB | 8.4 MB | 1.000 | 995 MB | 200 MB | 1.009 | 1.008 |
| a50 | F | 15.9 / 18.0 MB | 2.1 MB | 1.000 | 1195 MB | 400 MB | 1.007 | 1.006 |
| nrun30 | A | 21.8 / 21.9 MB | 8.4 MB | 1.004 | 436 MB | 120 MB | 1.012 | 1.011 |
| nrun30 | F | 14.7 / 14.7 MB | 2.1 MB | 1.000 | 556 MB | 240 MB | 1.006 | 1.006 |
| fib30 | A | 8.3 / 8.3 MB | 8.4 MB | 1.012 | 122 MB | 120 MB | 1.017 | 1.017 |
| fib30 | F | 2.1 / 2.2 MB | 2.1 MB | 1.000 | 241 MB | 240 MB | 1.006 | 1.006 |
| dna_100m | A | 8.4 / 8.5 MB | 8.4 MB | 1.005 | 404 MB | 400 MB | 1.010 | 1.010 |
| dna_100m | F | 2.1 / 2.2 MB | 2.1 MB | 1.000 | 805 MB | 800 MB | 1.006 | 1.006 |
| genome_100m | A | 8.6 / 8.7 MB | 8.4 MB | 1.005 | 404 MB | 400 MB | 1.010 | - |
| genome_100m | F | 2.1 / 2.2 MB | 2.1 MB | 1.000 | 805 MB | 800 MB | 1.006 | - |

Peak of the system-wide Dirty over every run of the groups `adv`, `shm`, `par`, `slots` (default or 2-8 slots), `a200`, `slope` and `nextsym`: **27.6 MB** (N run 30 M, layout A + PSV + NSV with the concurrent NSV pass; the independent sampler 25.5 MB; baselines 0-0.23 MB), layout A alone at most 23.9 MB (N run), a^n 200 M 17.0 MB. **Gate <= 64 MB: met with a factor of 2.3**, by layout A (8.4 MB of mapped dirty pages = 8 windows of 1 MiB), layout F (2.1 MB = 2 windows) and the spill (16 MB between fsyncs). Dirty is bounded by the design (`slots x window` + the fsync interval), not by the input. Without the spill fsync the spill file is the only thing that is dirty: **801.6 MB for a^n 50 M and 322.8 MB for the N run** (`e5_nosp.tsv`, no flush of the windows either: the same), because the kernel's `dirty_background_ratio` of this VM is 10 % of 15 GB.

### 6.3 Write amplification (two definitions)

* **Definition 1, pages dirtied:** `(write_bytes of /proc/self/io - the spill bytes written) / final size of the output file(s)`. `write_bytes` counts every page that goes from clean to dirty, so a page that is flushed and written again counts twice. 1.000-1.012 (table).
* **Definition 2, bytes at the block device:** the sectors written (`/sys/block/vda/stat`, field 7) from the start of the run to the end of the `sync` that follows it, minus the spill bytes written when the spill is on the disk, divided by the final size; the same without subtracting when the spill is on tmpfs (`/dev/shm`, the isolation runs: nothing but the output file reaches the device). 1.006-1.017 (table); the 0.6-1.7 % over 1 are the filesystem's metadata and journal.

**Both definitions are below 2 on every input** (gate: < 2). **No layout A thrash**: the abandon rule of the ADR row E5 does not fire; layout F is not "used as the ADR prescribes" because of a failure, it is measured next to it (section 7.3).

**How close layout A is to thrashing: the number of windows.** `--slots k` with k = 1, 2, 4 instead of 8 (one run each, cold, `e5_slots.tsv`; windows of 1 MiB, flushed when released):

| input | 1 window | 2 windows | 4 windows | 8 windows (default, 3 runs) |
|---|---|---|---|---|
| a^n 50 M | WA 1.00 (2.49 s) | 1.00 (2.65 s) | 1.00 (2.30 s) | 1.00 (2.43 s) |
| N run 30 M | **1.51** (1.33 s) | 1.02 (1.49 s) | 1.00 (1.33 s) | 1.004 (1.25 s) |
| Fibonacci 30 M | **2.83** (1.01 s) | 1.19 (0.65 s) | 1.05 (0.67 s) | 1.012 (0.68 s) |
| dna_100m | **2.56** (3.58 s) | 1.06 (2.41 s) | 1.02 (2.37 s) | 1.005 (2.50 s) |
| genome-like 100 M | **2.77** (3.79 s) | 1.07 (2.64 s) | 1.02 (2.61 s) | 1.005 (2.60 s) |

One window is not enough on every input except a^n (the writes alternate between the cursor and `previous` / `lb` a few windows back: each alternation flushes a page and dirties it again, up to 2.8 amplification); two windows keep the amplification at 1.02-1.19, four at 1.00-1.05. The hot region is therefore small (cursor + one or two older windows + the backward sweep of a pop burst), as the ADR argued for this layout; the eight windows of the default are a safety margin that costs 8 MB of mapped (not anonymous) memory.

### 6.4 Time against PH4 on the same input

"PH5 time must not exceed the PH4 time on the same input". PH4 here = the base mode of E4 re-measured in the same session on the same inputs (`e5_ph4ref.tsv`, scale-down by n / 3x10^9, 4 threads, cold, 3 runs; a^n 50 M 1.44 s, dna_100m 5.11 s: the E4 figures 1.75 and 5.21 s). PH5 = the forward scan including the reading of the LCP file from the cold page cache, plus the NSV pass for layout F (a second pass after it, or concurrent with it: `--parallel-nsv`). In brackets the ratio to PH4; **bold = above 1**.

| input | PH4 s | A | A, spill on tmpfs | F, NSV after | F, NSV concurrent | A + PSV + NSV (both), NSV after | both, concurrent |
|---|---|---|---|---|---|---|---|
| a10 | 0.41 | **0.63 (1.53)** | - | **0.71 (1.70)** | - | **0.85 (2.05)** | - |
| a50 | 1.44 | **2.43 (1.68)** | 1.30 (0.90) | **2.94 (2.04)** | **2.68 (1.86)** | **3.43 (2.38)** | **3.17 (2.20)** |
| nrun30 | 1.48 | 1.25 (0.84) | 0.92 (0.62) | **1.70 (1.15)** | 1.47 (0.99) | **1.89 (1.27)** | **1.70 (1.15)** |
| fib10 | 0.67 | 0.36 (0.53) | - | 0.48 (0.72) | - | 0.55 (0.82) | - |
| fib30 | 1.54 | 0.68 (0.44) | 0.67 (0.43) | 1.15 (0.75) | 0.91 (0.59) | 1.39 (0.90) | 1.14 (0.74) |
| rep10 | 0.68 | 0.39 (0.57) | - | 0.59 (0.87) | - | **0.70 (1.02)** | - |
| dna10 | 0.60 | 0.39 (0.66) | - | 0.59 (0.98) | - | **0.67 (1.13)** | - |
| dna_100m | 5.11 | 2.50 (0.49) | 2.36 (0.46) | 4.69 (0.92) | 3.29 (0.64) | **5.41 (1.06)** | 4.07 (0.80) |
| genome_100m | 4.92 | 2.60 (0.53) | - | 4.71 (0.96) | 3.37 (0.68) | **5.54 (1.12)** | 4.32 (0.88) |

* **Layout A: met on every input except a^n**: 0.44-0.84 of PH4 on N run, Fibonacci, rep_10m, dna_10m and the 100 Mbp inputs, **1.53-1.68 on a^n** (a^n 10 M, 50 M). The cause is measured: the spill (a^n 50 M: the scan alone, no output file and the spill on tmpfs, 0.87 s; the CHILD writes behind the cursor add 0.43 s (1.30 s); the spill on the disk with its fsyncs adds 1.13 s (2.43 s)). With the spill file on tmpfs (isolation runs) a^n 50 M takes 1.30 s (0.90 of PH4), with the spill on the disk not `fsync`ed 1.81 s (no flush of the windows either 1.51 s); the fsync every 16 MB that bounds the dirty pages costs 0.6-1.0 s on 793 MB of spill. PH4 on a^n is cheap (29 ns per symbol at this scale: the PLCP loop does one vector compare per position), PH5 on a^n moves 16 B per entry out and back (1.6 GB) and writes 4 B/char at the end behind the cursor.
* **Layout F: met on the 100 Mbp inputs (0.92 and 0.96, NSV after) and Fibonacci, not on a^n (2.04 of PH4) and the N run (1.15)**; with the NSV pass concurrent 0.64 and 0.68 at 100 Mbp, 0.59 on Fibonacci, 0.99 on the N run, still above 1 on a^n (1.86). The second pass is 14-26 ns per entry (a^n 50 M 0.71 s, dna_100m 2.29 s).
* **All three writes in one scan (A + PSV) and the NSV pass after it** are above PH4 on a^n (2.05-2.38), the N run (1.27), rep_10m (1.02), dna_10m (1.13) and at 100 Mbp (1.06, 1.12); with the NSV pass concurrent 0.80-0.88 at 100 Mbp, 0.74 on Fibonacci, still above 1 on a^n (2.20) and the N run (1.15).
* **Fused with PH4** (the hook): the scan adds its own time to B3 (PH4 alone 5.11 s at dna_100m, 7.50 s fused with layout A: +2.4 s = 0.47 of PH4; genome-like 4.92 -> 7.52 s; **a^n 50 M 1.44 -> 3.86 s fused, +2.4 s = 1.7 of PH4**). Fusing saves the LCP read (the scan reads no file) but PH4's 12 B/char of tuple scratch is still being written back by the kernel (1.2 GB dirty at the start of B3): the msync of the windows waits behind it, and one of three genome-like fused runs of layout F took 14.8 s instead of 7.1 s. The dirty pages of the scan alone (Dirty after the start of the scan minus Dirty at that moment) are 0.2-9.5 MB, the per-process mapped dirty pages 3.1-10.1 MB.

**G3: spill volume and write amplification PASS, dirty pages PASS, time PASS for layout A except a^n, for layout F on the 100 Mbp inputs and Fibonacci; time FAILS on a^n (A 1.5-1.7x, F 1.7-2.0x), on the N run for layout F (1.15x) and for the all-in-one pass (1.02-1.13x on the 10 M and 100 Mbp inputs, 1.27x on the N run, 2.1-2.4x on a^n), with the cause measured (spill I/O; a second pass over the LCP). The ADR's optional run-length encoding of arithmetic progressions on the stack, which makes a^n and N runs O(1) entries, would remove the spill on a^n; no gate except this time clause on a^n requires it, it is not built.**

## 7. G4 report (no gate)

### 7.1 Per-phase times (4 vCPU, cold cache, PH5 is one thread; medians of 3 runs)

"scan alone" = `--layout stats`: no output file, the spill on tmpfs: the cost of the stack and of the interval scan itself (the ADR estimate was 3-10 ns per entry). The reading of the LCP file from the cold page cache is the difference `forward_s - scan_s` (100 MB of U8 in 0.21 s, 200 MB of U32 in 0.64 s). Layout A adds the writes behind the cursor (8 windows, flushed when released) and, on a^n / N runs, the spill on the disk; layout F's forward half adds the PSV writes, the second half is the NSV pass.

| input | n | PH4 s | scan alone, no output (spill on tmpfs) | A: forward_s (read + scan + writes) | A: scan_s (ns per entry) | A: LCP read | F: forward_s | F: NSV pass | F total, NSV after / concurrent |
|---|---|---|---|---|---|---|---|---|---|
| a50 | 50 M | 1.44 | 0.87 (14.4 ns) | 2.43 | 1.79 (35.8 ns) | 0.64 | 2.11 | 0.71 | 2.94 / 2.68 |
| nrun30 | 30 M | 1.48 | 0.57 (15.8 ns) | 1.25 | 1.01 (33.6 ns) | 0.24 | 1.14 | 0.55 | 1.70 / 1.47 |
| fib30 | 30 M | 1.54 | 0.41 (10.4 ns) | 0.68 | 0.60 (19.9 ns) | 0.08 | 0.62 | 0.53 | 1.15 / 0.91 |
| dna_100m | 100 M | 5.11 | 1.56 (13.5 ns) | 2.50 | 2.28 (22.8 ns) | 0.21 | 2.38 | 2.28 | 4.69 / 3.29 |
| genome_100m | 100 M | 4.92 | 1.66 (14.6 ns) | 2.60 | 2.40 (24.0 ns) | 0.21 | 2.45 | 2.21 | 4.71 / 3.37 |

* **ns per entry:** the scan alone 10.4 (Fibonacci) to 15.8 (N run) ns, 13.5 ns on random DNA (whose stack is 18 deep and whose pops are data dependent); layout A 22.8-24.0 ns on the 100 Mbp inputs (scan) and 25.0-26.0 ns with the LCP read; layout F 21.6-22.3 (forward) + 22.1-22.9 (NSV) ns; the deep inputs 33-36 ns (a^n, N run: the spill). **Linear in n:** a^n 10 M / 50 M / 200 M (layout A): 63 / 49 / 47 ns per entry including the read (0.63 s, 2.43 s, 9.4 s).
* **Stack height maxima** (forward stack, 16 B): a^n 10 / 50 / 200 M 10 / 50 / 200 M entries (160 MB / 800 MB / 3.2 GB), N run 30 M 20 M (320 MB), genome-like 100 M 500,000 (8 MB), Fibonacci 10 M / 30 M 38 / 39, rep_10m 64, random DNA 15-18; NSV stack (8 B) 1 on a^n, 15-18 on the N run, Fibonacci and random DNA, 70 on rep_10m, 97 on genome-like 100 M.
* **Spill file size and I/O:** section 6.1. The explicit I/O of a^n 50 M, layout A: 993 MB read (LCP 200 MB + spill 793 MB) and 793 MB written, plus the 200 MB CHILD through the mapped windows; layout F: 1,193 MB read (two LCP reads of 200 MB + the spill) and 793 MB written plus 400 MB through the windows. On 100 Mbp random DNA: A reads 100 MB (the LCP), F 200 MB (twice), nothing is written explicitly.
* **RSS and arena peaks:** section 5 (11-24 MB RssAnon, up to 33 MB with the concurrent NSV pass; arena 8.39 MB, 12.58 MB with the concurrent NSV pass (two stacks)); mapped file pages resident (RssFile) 30-38 MB (windows + the text, not counted).
* **The cost of layout F's second pass:** 0.53-0.71 s on the deep and Fibonacci inputs (14-19 ns per entry), 2.2-2.3 s at 100 Mbp (22-23 ns): it is a second sequential read of the LCP (U8 / U16 / U32 plus the overflow pairs through a cursor that moves down the pair file) and 4 B/char of sequential descending writes, 0.9-1.0 of a forward pass; run on its own thread (`--parallel-nsv`) it saves 1.4 s of 4.7 s at 100 Mbp (4.69 -> 3.29 s; the two threads also share the disk and the memory bus), 0.3 s on a^n 50 M.
* **Escape / overflow handling:** the compact LCP is read through the base array and a cursor over the overflow pairs (forward and backward, 2 pairs per 16-byte window in the tests, 1 MiB at run time); genome-like 100 M has U8 with 1,390,028 overflow pairs (11.1 MB): its read share is 0.206 s against 0.212 s for random DNA of the same size, i.e. the escapes cost nothing measurable; a^n, the N run and Fibonacci are U32 (no overflow), rep_10m U16 (none). A pair that does not match its escape slot, a missing or an extra pair is refused (`InvalidDataException`, tested forward and backward).

### 7.2 NEXTSYM (optional pass, report only; `e5_nextsym.tsv`)

`NEXTSYM[r] = text[SA[r] + LCP[r]]` (the dense symbol, 0xFF at the end of the text): one sequential read of the SA and of the LCP, **one random read of the text per rank**, one byte written per rank through two mapped windows. Measured after the scan-only pass (`--layout stats --nextsym`), cold cache; equal to the in-RAM computation in 5 of 5 runs (`e5-oracle --sa` writes the reference). Medians of 3 runs, the share is NEXTSYM / (NEXTSYM + the scan-only pass of the same run):

| input | n | NEXTSYM s | ns per rank | share | RssAnon | scan-only pass s |
|---|---|---|---|---|---|---|
| dna10 | 10 M | 0.11 | 10.6 | 30 % | 20.9 MB | 0.25 |
| a50 | 50 M | 0.31 | 6.2 | 13 % | 14.6 MB | 2.12 |
| nrun30 | 30 M | 0.21 | 7.1 | 19 % | 14.6 MB | 0.97 |
| dna_100m | 100 M | 1.78 | 17.8 | 54 % | 20.9 MB | 1.49 |
| genome_100m | 100 M | 1.94 | 19.4 | 54 % | 23.0 MB | 1.65 |

On random DNA and the genome-like text of 100 Mbp the NEXTSYM pass costs as much as the whole scan (54 %: 17.8 ns per rank, the text of 25-37 MB is read at random after a cold start); on the deep inputs it is small (13-19 %: the text is small and the scan is slow there). It needs the SA file (4 B/char of sequential reads) next to the LCP: PH4-B3 of the base mode does not hold the SA, so fused it needs a second reading of the SA file (the bounded mode of E4 has the SA block in memory at the gather).

### 7.3 Both layouts side by side (a^n 50 M, N run 30 M, Fibonacci 30 M, dna_100m, genome-like 100 M; 3x10^9 columns are ESTIMATES)

| | layout A (tagged or packed) | layout F (PSV + NSV) |
|---|---|---|
| output | `<prefix>.cld`, 4 B/char (12 GB at 3x10^9) | `.psv` + `.nsv`, 8 B/char (24 GB) |
| writes | behind the cursor: next link at `previous`, up at `rb`, down at `lb`, 1 per rank; 8 windows | sequential: PSV forward, NSV backward; 2 + 2 windows |
| passes over the LCP | 1 (reads 1.0-4.0 B/char) | 2 (forward, and backward); the NSV pass can run concurrently |
| mapped dirty pages (per process) | 8.4 MB (8 windows) | 2.1 MB (2 windows; 4.2 MB with the concurrent pass) |
| system Dirty, own / independent (a^n 50 M, N run, Fibonacci, dna_100m, genome) | 16.8/14.8, 21.8/21.9, 8.3/8.3, 8.4/8.5, 8.6/8.7 MB | 15.9/18.0, 14.7/14.7, 2.1/2.2, 2.1/2.2, 2.1/2.2 MB |
| write amplification, definition 1 / 2 (range over the 5 inputs) | 1.000-1.012 / 1.008-1.017 (tmpfs spill 1.008-1.017) | 1.000 / 1.006 |
| with 1 / 2 / 4 windows (WA 1, worst input) | 2.83 / 1.19 / 1.05 (Fibonacci) | n/a (every write is sequential: 1 window would do) |
| spill | 16 B entries: 16 B x (height - resident) | forward stack 16 B (PSV pass) + NSV stack 8 B (0 on a^n) |
| time, 100 Mbp (dna_100m / genome), cold | 2.50 / 2.60 s (0.49 / 0.53 of PH4) | 4.69 / 4.71 s after, 3.29 / 3.37 s concurrent (0.92 / 0.96, 0.64 / 0.68 of PH4) |
| time, a^n 50 M / N run 30 M | 2.43 / 1.25 s (1.68 / 0.84 of PH4) | 2.94 / 1.70 s after (2.04 / 1.15), 2.68 / 1.47 s concurrent |
| RssAnon (all inputs) | 11-24 MB | 13-24 MB (arena 12.6 MB concurrent) |
| fused with PH4 (the hook) | yes: the scan is the consumer of PH4-B3; +2.4 s on dna_100m (0.47 of PH4) | forward half only (PSV); the NSV pass needs the LCP file PH4 wrote |
| reader of the structure | one field per rank, kind in bits 30-31 (tagged: texts <= 2^30 - 1) or told apart by LCP comparisons (packed) | needs the LCP values and, per child step, a chain of next smaller values (first l-index = first minimum of `LCP(lb..rb]`); 8 B/char |

ESTIMATE for 3x10^9 (linear in n at the measured ns per entry; the page cache of a 16 GB machine does not hold the files of that size, so time becomes max(CPU, I/O / disk bandwidth), not measured): layout A forward scan on random or genome-like text 75-78 s (25-26 ns x 3x10^9), CHILD 12 GB; layout F 141 s with the NSV pass after the forward scan, about 99 s with the concurrent pass, 24 GB; a^n or a long run: 146 s at 48.5 ns per entry with 48 GB of spill written and read back (a^n 200 M measured 9.4 s for 3.19 GB of spill, linear); dirty pages, RssAnon and the number of windows are constants (8.4-27.6 MB, 11-24 MB, 8 or 2 per file); fused with PH4 in the base mode PH5 adds its forward scan to B3: about 70-75 s (23-24 ns per entry in the hook) on top of the 230 s of PH4 (ESTIMATE of E4: 77 ns per symbol at 500 Mbp unscaled).

## 8. What was not verified, and the limits of the evidence

* **Sizes.** PH5 alone was run on random and genome-like text up to 100 Mbp, on a^n up to 200 M (against the closed form), on an N run of 30 M and on Fibonacci up to 30 M; nothing above (the prototype disk holds one 100 Mbp run with its reference files; no 500 Mbp run of PH5 was made), no text above 2^32 - 2 symbols (positions are `uint`; 40-bit positions would make the entry 20-32 B), every 3x10^9 figure is an ESTIMATE (linear extrapolation of the measured ns per entry, constant memory and dirty pages). No genome-like text above 100 Mbp, no repeat structure that was not generated here (the genome-like generator, a^n, an N run, Fibonacci, `rep_10m`; the tiny-text and family tests cover the rest at small sizes).
* **The oracle** is the in-RAM structure of the product (`ChildTable.Build`, the `LcpPyramid`, `EsaStatistics`, built from the same LCP file), not an independent implementation; the Ukkonen tree is the independent oracle, but only for the tiny texts (binary to length 8, ternary to length 5, seven named strings); a^n 200 M is checked against the closed form. The definitions of PSV and NSV (strict, none = -1 / n + 1) are those of `LcpPyramid`.
* **Layout F is written but not read by production code:** the enumeration of children through PSV / NSV exists in the test decoders only (it proves the structure sufficient); its cost per step (a chain of next smaller values, bounded by the node depth of the first child, not O(1)) is not measured. The ADR-001 candidate (c) (RMQ + PSV/NSV) needs an RMQ as well (PH6).
* **The disk is a virtio disk, probably host-cached**: the cost of `msync` and `fsync` (the a^n spill: 0.6-1.1 s on 793 MB) is that of this disk; a real SSD or a network disk would move it, not the structure. The page-cache dependence is not tested: every run is cold at the start and then has 15 GB (the 100 Mbp runs keep their 1.2 GB of output cached).
* **Dirty pages** are the whole system's (`/proc/meminfo`), with nothing else running; the fused runs cannot be read against the 64 MB gate because PH4 (E4, which does not flush its tuple files or its LCP windows) has 1.2 GB (100 Mbp) dirty at the start of B3: for them the Dirty above the level at the start of the scan (0-9.5 MB) and the per-process mapped dirty pages (3.1-10.1 MB, the PH4 LCP windows included) are given. The flush on release uses `MemoryMappedViewAccessor.Flush` (msync with `MS_SYNC` on this kernel: Dirty falls from 801 MB to 17 MB when the spill is fsynced and the windows flushed, section 6.2); other platforms were not tried.
* **Fused mode:** one of three genome-like fused runs of layout F took 14.8 s instead of 7.1 s (the msync of the windows waits for the kernel's writeback of PH4's 1.2 GB of scratch); the fused runs use PH4's base mode with the scale-down emulation of E4, not real-size chunks; B3 stays single-threaded (the scan runs on its thread).
* **Threads:** PH5 alone is one thread; the concurrent NSV pass is two (a 4 vCPU machine; 3 runs per point); no run with fewer than 4 vCPUs.
* **Independent sampler:** `sample_proc.py` (Python, 20 ms, a separate process) against the program's samplers over 198 runs: RssAnon peak ext - own between -0.05 and +6.5 MB (median 0.0, 89 % within 0.5 MB; the larger ones are start-up samples of the runtime that the program's sampler, which starts inside the process, does not see), Dirty peak between -4.0 and +4.3 MB (different sample instants and baselines); it shares no code with the program but reads the same kernel counters.
* **No resume / crash safety**, no `fsync` of the outputs in the timing (reported separately, below 1 ms), no container writer, **no decision between the tagged and the packed codec** (both verified; the tagged codec holds ranks in 30 bits, so texts up to 2^30 - 1 symbols), no PH6.
* **NEXTSYM** needs symbol codes below 255 (one byte per rank, 0xFF = end of text) and was measured with the scan-only pass before it; it was not fused.
* **Re-measured PH4** (the reference of the time clause) is the base mode of E4 with its scale-down emulation, 3 runs in this session; it is not the PH4 time at the real size (E4: 77 ns per symbol at 500 Mbp unscaled, 52 ns per symbol at 100 Mbp scaled): the PH5 / PH4 ratios would be 0.7 times smaller against real-size PH4.
* The first run of the protocol was aborted after the `a10` rows of the group `adv` because two inputs were symbolic links (a program that asks for `FileInfo.Length` sees the length of the link); the rows of that attempt and the `ph4ref` rows were removed and the groups repeated; the setup script now copies those files.

## 9. Open points for E6 and 1b.7

* **How the PH4 hook and PH5 compose.** `ChildScan` is an `ILcpChunkConsumer`: it receives the dense resolved values of B3, rank chunk by chunk, in ascending order, and is created at the first chunk so that its blocks and windows do not overlap PH4's B2 peak (fused RssAnon = PH4's: 51.0 MB at 100 Mbp scaled). It adds 23-24 ns per entry to B3's single thread: +2.4 s on 5.1 s of PH4 at 100 Mbp (+47 %), +2.4 s on 1.4 s on a^n 50 M. To hide it, B3 could hand the chunks to the scan on a second thread through a two-slot queue (the scan is sequential by definition, the B3 writes are not on its path); not built. Layout F's NSV half **cannot** be fused with an ascending B3 (it needs the LCP in descending order): it runs from the LCP file PH4 wrote (a second read of 1.0-4.0 B/char); if E6 wants it fused, B3 would have to visit the rank chunks in descending order (the rank chunk files are independent) and keep them for the forward scan (+8 B/char of scratch), or the forward scan runs from the file as well.
* **NEXTSYM needs `SA[i]`** (and the text at random): the base mode's B3 does not hold the SA, so fused it would read the SA file again (4 B/char sequential) next to the text reads (17.8 ns per rank at 100 Mbp, as much as the whole scan); the bounded mode's B3' has the SA block at hand (5 lines to pass it). NEXTSYM is also the only PH5 output that does not stream the LCP at all.
* **Codec of layout A** (open, owner): Tagged (kind in bits 30-31, no LCP comparison on a read, texts to 2^30 - 1; 4 B/char) against Packed (raw value; reads compare LCP values; 4 B/char, any text below 2^32 - 2; fewer next-link writes: 0.85-1.00 writes per rank against 1.000 for tagged on the nine inputs, `e5_adv.tsv` `*.apacked`). Both are produced by the same scan and verified; the windows touched, the dirty pages and the amplification are the same for both (the packed run of each input maps the same windows).
* **Layout A or F for the container:** A is 4 B/char, one pass, and its writes behind the cursor do not thrash with 8 windows (1.00-1.012) or even 2-4 (1.00-1.19); F is 8 B/char, two passes, every write sequential (2 windows, 2.1 MB dirty) and needs an RMQ and a reader that follows next-smaller chains. Both pass the dirty and the amplification gates; time against PH4 is better for A (0.44-0.84 except a^n); the choice is the owner's, no gate decides it.
* **The spill and a^n:** 16 B per stack entry; the ADR's run-length encoding of arithmetic progressions `(l0, lb0, first0, count)` would make a^n and N runs O(1) entries and take 0.6-1.1 s (the spill) out of the 2.4 s of a^n 50 M, and at 3x10^9 48 GB of spill and its fsyncs; it is not built because no gate except the time clause on a^n needs it. The spill is accessed only at the end of its file; a read-ahead thread would hide the refill latency.
* **Windows and flushes are parameters** of `ChildBuildOptions`: 8 windows of 1 MiB for CHILD (>= 2 needed), 2 for PSV and NSV, flush on release, fsync of the spill every 16 MB. 1b.7 should place the CHILD file at a section offset (the writer takes a file of its own length now), write it at its final length first (the section table is written last), and take the CRC of the section in a pass or from the scan (the scan does not see the written bytes).
* **40-bit positions:** the entries become 20-32 B (the spill 20-32 B x the chain depth); `uint` positions here.
* **PSV / NSV in the pack:** the PH6 row of the ADR lists PSV and NSV for layout F "2 B/char of I/O": the two arrays are 8 B/char on the disk here (two `uint` arrays); a 40- or 32-bit packing was not tried.

## 10. Mutation audit of the new code

Hand-made single-site mutants of the new production code: `SpillStack` (20), `AlignedWindowWriter` (14), `CompactLcpFile` (13), `ChildScan` (26), `NsvScan` (6), `ExternalChild` (4), `NextSymScan` (5), `Experiments/E5` (12), `Harness/DirtySampler` (4) = **104 in round 1**, plus 16 replacements of mutants that did not compile and 6 re-runs in round 2, 2 re-runs in round 3: lists `analysis/mutation/e5_*_mutants.json` (generators `make_e5_mutants.py`, `make_e5_round2_mutants.py`, which also check that every pattern is unique in its file), runner `run_e5_mutants.py` with a **watchdog** (build limit 300 s, test limit 240 s; the whole process group, test host included, is killed on expiry, and the mutant counts as killed(hang)), driver `run_e5_audit.sh` (two audits at a time, each in its own git worktree, so the checkout is never mutated; the worktrees were removed afterwards), raw results `results_e5_*.tsv`, merged table with verdicts and notes `results_e5.tsv` (`merge_e5.py`). Stage 1 runs the Ph5 tests without the slow ones (the exhaustive enumeration, the families at 10^5 and 10^6); a mutant that survives it is run again against the whole Ph5 namespace.

| | mutants |
|---|---|
| listed (round 1) | 104 |
| did not compile (analyzer errors: S1144 unused member 4, CS0649 field never assigned 6, constant conditions S1125 / CS0162 / S4487 4, nullable CS8604 / CS8602 2), replaced by the mutants with suffix `b` in round 2 | 16 |
| **compiled in round 1** | **88**: 80 killed (1 by a test-host crash, S01; none by the watchdog), **8 survived** |
| round 2 (head `602195bd`: the 16 replacements and the 6 survivors that the strengthened tests should now kill) | 22: 21 killed (S03b, S04b, S16b, S17b, S19b, W09b, W13b, W14b, A07b, A17b, X04b, E07b, E12b, D02b, D04b and the re-runs S05, C11, X03, E08, E11, M03), **1 survived (W07b)** |
| round 3 (the test for W07b, the test for M04) | W07b killed (by a host crash: a write that leaves its mapping), M04 killed |
| **final: listed 120, compiled 104** | **killed 103**, equivalent argued 1 (S06), **survivors 0**, not compiled 16 (each replaced) |

The eight survivors of round 1 and what closed them: **S05** (a block that was read back counts as dirty, so it is written again): the first test of the clean-block rule only asserted that *some* eviction happened; the new test counts the writes exactly (a refilled, unmodified block is evicted with `CleanEvictions` 1 and `BlocksWritten` unchanged, a modified one with a write and the modification surviving the next refill). **S06** is equivalent (the new top block's dirty flag is set by the push that always follows `Advance`). **C11** (the base file may be longer than `(n + 1) x width`): the tests had only files that were too short; files of 7, 11, 13, 25 bytes are refused now. **X03** (the volume of the LCP files forgot the overflow file): the sample texts had none; the test uses a text with 700 escapes. **E08** (LRS_RANK ignored by the comparison with the oracle) and **E11** (a differing NEXTSYM did not fail the run): one test per statistic with a wrong value, one test with a corrupted `.nsym` reference. **M03** (a 9-bit text accepted): the wrong length of the SA file in the test refused it first; the test gives the right length. **M04** (an SA file one entry short or long accepted): files of 24 and 32 bytes for n = 6. **W07b** (a piece that runs past the end of its window is not cut there) survived round 2 because windows smaller than a page are mapped with the page around them (the bytes land in the right place in the file); the test with windows of a page (4,096 bytes) crosses them, and the mutant now leaves its mapping (host crash).

Limits of the audit: mutants are hand-made and single-site; the reporting rows of `e5` (the times, I/O, write-amplification and dirty-page numbers) are checked by a few strings only (`parameters`, `n`, the `*_equals_oracle` rows, the flush counts, the exit code); the older E0-E4 code was not audited; stage 1 excludes the exhaustive builds and the families at 10^5 / 10^6, so a mutant that only they could kill would have shown as a survivor of stage 1 and been run against the whole namespace (none needed it: every mutant of the scan was killed in stage 1).
