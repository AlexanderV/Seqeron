# REVIEW-S1a-opt: independent review of the optimized (unsafe) kernels, iterations 1o.1b, 1o.1c, 1o.2, 1o.3a, 1o.3b

Reviewer: independent agent (did not write the code). Date 2026-11-04. Reviewed: branch `claude/stoic-maxwell-0olr2z` at `984f86ec` (clean, in sync with origin; the unsafe kernels are unchanged by this review: `git diff -- src` is empty). Read: [PROCESS.md](../PROCESS.md), [PROGRESS-S1.md](../PROGRESS-S1.md) Log 1o.0-1o.4a, [REVIEW-S1a.md](REVIEW-S1a.md), [ADR-002](ADR-002-sa-construction.md) section 5. Code: `src/SuffixTree/Algorithms/SuffixTree.Esa` (`SaisKernel*`, `SaisInduction`, `TypedText`, `PtrText`, `SaisPrefetch`, `LcpKernel*`, `Memory/Buf.cs`, `Memory/NativeArena.cs`), tests `tests/SuffixTree/SuffixTree.Esa.Tests/{Kernels,Review,Safety}`, harness `apps/SuffixTree.Esa.Bench`. Raw data of every number below: [`analysis/review_opt/`](../analysis/review_opt/) (TSV/log files, named in the text). Machine: 4 vCPU Xeon (KVM), 15 GB, .NET SDK 10.0.403, Release, tiered PGO default.

## Verdict

**0 BLOCKER, 2 MAJOR, 8 MINOR. No wrong SA, wrong LCP, out-of-bounds access, read-before-write or data race was found in the optimized kernels.** Every suffix array and LCP of the optimized chain that I compared with libsais (28 texts at 10 Mbp and 26 at 1 Mbp, 378 builds through three chains) is byte-identical (SHA-256 of SA and LCP); three new fuzz seeds (101-103, about 69 million comparisons) found nothing; five effective runtime configurations (hardware intrinsics off, AVX2 off, SSE4.2/SSE4.1/POPCNT off, tiered PGO off, tiering off) plus three knobs that turned out to be no-ops pass default, Exhaustive, Vectors and Fuzz. The Log's numbers reproduce (section 1). What the review did find is (a) a performance cliff of the 1o.3b LCP kernel that makes its headline numbers valid only for repeat-free text (MAJ-1), and (b) holes in the safety net itself: four of the seven bounds checks of the checked `Buf<T>` survived every test, including the checked build (MAJ-2). The holes are closed by added tests; the cliff is not a defect but a decision input.

**Verdict: further optimization iterations and the S1b work MAY build on these kernels**, on four conditions: (C1) the 16 new tests of this review stay in the default suite (they are: 273 tests, raw and checked); (C2) S1b sizes memory and time from the repeat-rich rows of section 5 (peak 18.4 B/char for `U32` LCP, 20.2 B/char in the worst measured fallback), not from the 12.4 B/char / 1.28x random-DNA figures; (C3) the next optimization iterations are judged on a corpus that contains a repeat-bearing text (`rep_10m` and an embedded-repeat text) next to `dna_10m`, because on those texts 1o.3b contributes nothing and the LCP phase is 25-40 % of the build; (C4) the checked-pointer run (`-p:EsaBounds=true`) and the `Exhaustive`/`Vectors`/`Fuzz` categories are wired into an automatic run before the first release (S8), since the default CI job (`dotnet test Seqeron.sln -c Release`) runs none of them.

| id | one line |
|---|---|
| MAJ-1 | The direct LCP scan (1o.3b) is abandoned for any text that contains ONE repeat longer than about 8 sqrt(n) symbols (10,000 letters at 1 Mbp, 25,000 at 10 Mbp); the LCP phase then runs 1.6-2.1x slower than the three-pass kernel it falls back to (whole build 1.17-1.33x) and 1.1-1.6x slower than the safe reference. The Log's -23 % / 1.28x libsais / 12.4 B/char hold for repeat-free random DNA only |
| MAJ-2 | The checked accessor is the only detector of out-of-bounds READS, and four of its seven bounds checks (`ReadUnaligned` upper bound and lower bound, `LoadVector128` upper and lower bound) were wrong or missing without any test or the checked build noticing (mutants R82-R85) |
| MIN-1 | CI runs only `dotnet test -c Release`: no checked build, no Exhaustive/Vectors/Fuzz |
| MIN-2 | Performance-only mutations (prefetch address, distance, give-up fraction) survive every correctness test; 7 of 17 known ones are now pinned by new tests, 10 cannot be |
| MIN-3 | The managed output arrays (`sa`, `lcp`, `overflow*`) are neither guarded by canaries nor poisoned by the harness; the SA span is validated once and re-read later (contract: caller must not mutate it) |
| MIN-4 | The fallback of 1o.3b allocates `n/4 + n/8` ints and a discarded `n`-byte LCP array first: +2.2 B/char peak (18.6-20.2 B/char) on texts that fall back |
| MIN-5 | The Log's headline ratios are single-corpus: ESA / libsais is 1.23-1.35x on random DNA and DNA with a small repeat but 1.5-2.5x on the repeat-rich 10 Mbp texts I built |
| MIN-6 | The documented way to reduce the instruction set does not work on .NET 10: `DOTNET_EnableSSE41`, `EnableBMI1`, `EnableLZCNT`, `EnablePOPCNT`, `EnableSSE` are ignored; only `EnableHWIntrinsic`, `EnableAVX/AVX2` and `EnableSSE42` take effect |
| MIN-7 | ARM64, 32-bit and big-endian paths are untested here (prefetch is a no-op there; `CommonPrefix` has an endian branch); `LcpKernel.cs:119` claims `AggressiveOptimization` that no method carries |
| MIN-8 | An invalid but in-range suffix array is memory safe and linear time in both LCP kernels (verified, new test) but gives garbage: `Build(PreparedText, sa)` is not a validator |

## 1. Reproduction of the Log (item 6)

Quiet machine: `uptime` load 0.20 at the start of both sessions (the runner refuses above 0.30), nothing else running, build servers shut down; libsais built from the pinned commit `1aae6344` by a driver outside the repository (`analysis/review_opt/libsais_driver.c`, only calls the public API), same-session SA+PLCP+LCP time. Paired interleaved A/B runner of the harness (`ab`), 15 pairs, 9 A/A calibration pairs, fresh process per run, every child reproduced the libsais vectors (`vector_match = 1`). Raw: `ab_A_reference_vs_optimized_dna_10m.txt`, `ab_B_threepass_vs_direct_dna_10m.txt`.

| quantity (`dna_10m`) | Log | session A: A = safe reference chain, B = optimized chain | session B: A = three-pass LCP kernel, B = direct LCP kernel |
|---|---|---|---|
| SA+LCP total, B median (min-max) | 0.741 s | **0.746 s** (0.707-0.840) | **0.795 s** (0.731-0.844) |
| SA+LCP total, A median | - | 1.495 s (1.464-1.585) | 1.006 s (0.956-1.045) |
| paired B/A median [95 % interval] | 0.7707 [0.7596, 0.7768] (1o.3b) | **0.4986** [0.4917, 0.5009]; A/A band [0.973, 1.027] | **0.7890** [0.7532, 0.8038]; A/A band [0.978, 1.022] |
| LCP phase B/A | 0.2420 | 0.1994 | 0.2640 |
| SA phase B/A | 0.9894 (same code) | 0.5937 | 1.0131 (same code) |
| same-session libsais SA+PLCP+LCP | 0.581-0.601 s | 0.589 s (pre 0.589, post 0.592) | 0.583 s (pre 0.583, post 0.584) |
| ESA (B) / libsais | 1.28x | **1.27x** | **1.36x** |
| peak RSS (VmHWM) | 12.4 B/char | 123.4 MB = **12.3 B/char** (A: 153.6 MB = 15.4) | 123.4 MB = **12.3 B/char** |
| verdict of the runner | ACCEPT | ACCEPT | ACCEPT |

The numbers reproduce. Session B's 1.36x comes from a slower session for the managed code (SA phase 0.711 s against 0.66 s in the Log's 1o.3b session; libsais itself ran at 0.583 s), the known between-session spread of this VM (5.9 %, REVIEW-S1a MAJ-2); within a session the A/B ratios agree with the Log within 2 points (0.789 against 0.771). Cold single-process builds from my worst-case campaign give the same picture: `dna_10m` direct 0.788 s, three-pass 0.990 s, libsais 0.635 s (1.24x) in that session (section 5).

Machine note for later sessions: libsais ran at 0.58-0.64 s in my sessions against 0.543-0.601 s in the Log's; compare ratios, not absolute seconds.

## 2. Findings

### MAJ-1 The direct LCP scan is lost on any text with a long repeat; the 1o.3b headline is valid only for repeat-free text

* **Claim (PROGRESS-S1 1o.3b).** "Bounded, then linear in the worst case ... LCP phase 0.0675 s, SA+LCP -23 %, 12.4 B/char"; "`rep_10m` ... no regression beyond about 5 %".
* **Mechanism (`LcpKernel.Direct.cs:41,77,172`).** A rank whose LCP exceeds the vector width is resolved vector by vector; the total number of full-vector steps is limited to `2 n`. A repeat of length L that occurs twice contributes about L ranks with LCP between 16 and L, i.e. about `L^2 / 32` steps. The budget is exhausted as soon as `L^2 / 32 > 2 n`, i.e. `L > 8 sqrt(n)`: 8,000 symbols at 1 Mbp, 25,000 at 10 Mbp, 80,000 at 100 Mbp. Then `TryDirect` returns null after the scan and the budget are spent and the three-pass kernel runs from the start.
* **Evidence.** Random DNA with ONE block of 10,000 letters (1 %) occurring twice, n = 10^6: the attempt gives up (explicit failing test `Direct_Succeeds_OnRandomDnaWithOneEmbeddedRepeatOfOnePercent`). At 10^7 the threshold is exactly where predicted: a 15,000-letter repeat (0.15 %) succeeds (LCP phase 0.113 s, whole build 0.831 s), a 40,000-letter repeat (0.40 %) fails (LCP 0.470 s against 0.286 s for the three-pass kernel, whole build 1.177 s against 0.995 s). Cold fresh-process rows, 10 Mbp, median of 3 (`cold_embedded_10m.tsv`, `cold_10m.tsv`):

| text (10 Mbp) | libsais SA+PLCP+LCP | ESA direct (shipped) | ESA three-pass | ESA / libsais | direct / three-pass | LCP phase: direct / three-pass / safe ref | LCP kind | peak RSS B/char | SA+LCP = libsais |
|---|---|---|---|---|---|---|---|---|---|
| random-dna | 0.605 | 0.746 | 0.951 | 1.23 | 0.78 | 0.073 / 0.277 / 0.437 | U8 | 12.3 | yes |
| dna_10m | 0.635 | 0.788 | 0.990 | 1.24 | 0.80 | 0.079 / 0.281 / 0.419 | U8 | 12.2 | yes |
| embedded-015pct | 0.614 | 0.831 | 1.049 | 1.35 | 0.79 | 0.113 / 0.311 / 0.387 | U8 | 12.8 | yes |
| dup400-p10 | 0.565 | 0.842 | 0.940 | 1.49 | 0.90 | 0.170 / 0.261 / 0.380 | U8 | 12.6 | yes |
| dup400-p15 | 0.565 | 0.937 | 1.019 | 1.66 | 0.92 | 0.248 / 0.287 / 0.348 | U8 | 12.9 | yes |
| rep_10m | 0.516 | 0.957 | 0.935 | 1.85 | 1.02 | 0.237 / 0.225 / 0.335 | U16 | 16.5 | yes |
| embedded-040pct | 0.609 | 1.177 | 0.995 | 1.93 | 1.18 | 0.470 / 0.286 / 0.376 | U8 | 16.6 | yes |
| rep-100x-1pct | 0.481 | 0.936 | 0.897 | 1.95 | 1.04 | 0.230 / 0.226 / 0.346 | U16 | 16.4 | yes |
| embedded-1pct | 0.600 | 1.179 | 1.007 | 1.97 | 1.17 | 0.464 / 0.296 / 0.363 | U8 | 16.6 | yes |
| embedded-20pct | 0.604 | 1.221 | 1.091 | 2.02 | 1.12 | 0.494 / 0.311 / 0.410 | U32 | 20.2 | yes |
| embedded-5pct | 0.630 | 1.301 | 1.083 | 2.07 | 1.20 | 0.506 / 0.315 / 0.459 | U8 | 17.0 | yes |
| dup400-p25 | 0.583 | 1.238 | 1.004 | 2.12 | 1.23 | 0.568 / 0.271 / 0.411 | U16 | 18.8 | yes |
| dup400-p20 | 0.564 | 1.210 | 0.911 | 2.15 | 1.33 | 0.518 / 0.249 / 0.336 | U16 | 18.6 | yes |

  The worst measured loss is `dup400-p20/p25` (random DNA in which 20-25 % of the 400-letter blocks are copies of earlier blocks): LCP phase 0.518-0.568 s for the shipped kernel against 0.249-0.271 s three-pass (2.1x) and 0.336-0.411 s for the safe `LcpReference` (1.5x); whole build 1.21-1.24 s against 0.91-1.00 s (1.23-1.33x). The bare failed attempt costs 0.19-0.33 s at 10 Mbp (the scan, 0.06-0.09 s, plus the exhausted budget: 2n sequential vector steps and up to n/4 random first touches), then the three-pass kernel runs. It is bounded (about 2.1x the three-pass LCP phase; I found nothing worse) and not super-linear in the sense of section 5, but it is a regression against the previous kernel and against the reference.
* **Why it matters.** `dna_10m`, the acceptance corpus of the stopping rule, is the one text on which 1o.3b works. Real genomes contain repeats far above 8 sqrt(n) (segmental duplications, satellite arrays, transposon-rich regions); `rep_10m` (100 kbp block x 100, 1 % substitutions) already falls back, and there the LCP phase is 0.237 s (25 % of the build) and ESA / libsais is 1.85x. So the next iterations (1o.4+, S1b memory budget) should not assume 12.4 B/char, 0.07 s LCP or 1.28x libsais on real data.
* **Fix proposals (not applied).** (1) Decide early from a sample: charge the budget per long rank and give up after the first ~4,096 long ranks if they averaged more than ~8 steps; that bounds the waste to the scan (about 0.07 s at 10 Mbp, direct/three-pass <= 1.3 instead of 2.1). (2) Make the resolution cheaper for long matches (32-byte vectors where `Avx2` is available; Kasai-style reuse restricted to the long ranks sorted by text position). (3) Keep the shipped design, state the limitation in the Log and add `rep_10m` plus an embedded-repeat text to the acceptance corpus (C3). Option (1) is a few lines and removes most of the loss.
* **Tests.** `Direct_OnAnEmbeddedLongRepeat_IsLinear_AndEqualsTheReference` (default, n = 4 million, 20 % repeat; unbounded resolution would take about 2e10 steps) and `Direct_OnDuplicatedBlocks_JustBelowTheQuarterLimit_...` pin the contract "result equals the reference, time linear"; the explicit test above documents the cliff.
* **Resolution (iteration 1o.3c, [PROGRESS-S1.md](../PROGRESS-S1.md) Log).** Fixed by resolving the long ranks in text order with the Kasai carry (O(n / 16 + long ranks) vector steps for any repeat length) when a cumulative step allowance of the per-rank loop is exhausted; the explicit test is now the default test `Direct_Succeeds_OnRandomDnaWithOneEmbeddedRepeatOfOnePercent` (unchanged assertion). Measured paired against the 1o.3b chain: whole build 0.67-0.74 on the embedded-repeat rows and dup400-p20/p25, LCP phase 0.17-0.35, dna_10m / dna_1m / rep_10m unchanged within the A/A band. The first attempt of 1o.3c (text order for everything) was rejected because its JIT cost made dna_1m 5.8 % slower.

### MAJ-2 Four of the seven bounds checks of the checked `Buf<T>` were unverified and wrong in mutants

* **Claim (REVIEW-S1a MAJ-1, 1o.0).** The checked build catches every out-of-bounds read; "the checked build is a required acceptance step" because Release survivors such as K14/K15/D05/D14 die only there.
* **Evidence.** I mutated the checked section of `Buf.cs` (mutants R80-R86). Against the full pre-existing suite, in the checked configuration: R80 (indexer admits `index == length`), R81 (`Slice` admits one element too many) and R86 (`AsSpan` admits one more) were killed (4, 1, 1 tests); **R82 (`ReadUnaligned` upper bound +1 byte), R83 (`ReadUnaligned` accepts negative indices), R84 (`LoadVector128` upper bound +1 element) and R85 (`LoadVector128` accepts negative indices) survived** (`mutants_summary.tsv`, column "pre-existing suite only"). These two accessors are exactly the ones used by `CommonPrefix` (8-byte words) and by `ScanDirect` / `ResolveLong` (16-byte vectors); with R84 in place, the kernel mutant D05 (vector load with 15 symbols left) would no longer be detected by the checked build.
* **Why it matters.** The checked build is the only detector of out-of-bounds reads (canaries see writes only). A silent slack in a bound check turns "the checked build is green" into a weaker statement than the Log makes.
* **Fix (applied as tests).** `CheckedBuf_ReadUnaligned_...`, `CheckedBuf_LoadVector128_...`, `CheckedBuf_Slice_AndAsSpan_...` pin every bound at the exact element/byte for byte and 16-bit buffers, exact fit accepted, one past rejected, negative rejected; all four mutants are now killed (1 test each). Remaining work: C4.

### MIN-1 Nothing automatic runs the checked build or the heavy categories
`.github/workflows/dotnet.yml:53-56` is `dotnet build Seqeron.sln -c Release` and `dotnet test Seqeron.sln -c Release`: the 273 default Esa tests in the raw configuration only. Ten of my mutants (R05, R13, R44 and R80-R86) survive that run and die only in `-p:EsaBounds=true`; Exhaustive, Vectors and Fuzz are `[Explicit]`. **Fix:** a CI job (or the S8 job) that runs `analysis/run_esa_checked.sh`, and nightly `Category=Exhaustive|Vectors|Fuzz`.

### MIN-2 Performance-only mutations are invisible to correctness tests
Mutants that change only where or when a prefetch is issued, or when the scan gives up, return the same arrays. Of the 17 such mutants known from the Log and this review (R20, R24, R25, R40, R43, R64, R87; P04-P07 of 1o.2; Q06-Q08 of 1o.3a; D23, D29, D30 of 1o.3b), the old suite killed none. `ReviewS1aOptPrefetchTests` now pins the prefetch target of `InduceL`, `InduceS` and `NameLmsSubstrings` (a recording `ITypedText` shows that the position prefetched in iteration i is exactly the one read `Induce`=32 / `Name`=16 iterations later, including the clamps at both ends and the plateau range of the constants): R20, R24, R40, R43 and P04, P05, P07 die (7 of 17; P01-P03 die in Release now, no checked build needed). What cannot be pinned because the prefetch goes through the concrete `Buf<T>` or only the moment of giving up changes: the LCP prefetches (Q06-Q08, D23), the name-slot prefetch (R25, P06), the `nta` hint (R87), the give-up fractions (R64, D29, D30). A silent 10-15 % slowdown from such a change is caught only by the paired A/B runner. **Fix:** a nightly A/B of the current tree against the last accepted commit; or counters in `LcpTrace` for the give-up point.

### MIN-3 Managed output arrays are unguarded and uninitialised; the SA span is trusted after one validation
`GC.AllocateUninitializedArray` returns stale heap memory above about 2 KB for the SA (`SaisKernel.cs:57`) and the compact LCP (`LcpKernel.cs:184,190,196`, `Direct.cs:62`). The registry poisons only arena blocks and puts canaries only around arena blocks, so (a) a cell the kernel forgets to write is invisible on a clean heap (R60, `dst[0]` never written, is killed by one pre-existing test and by two of the new ones, `ProductionEntryPoints_DoNotDependOnTheInitialContentOfTheirOutputArrays` among them), and (b) a Release write one past these arrays corrupts the managed heap without a canary (only the checked build sees it). Separately, `ScanDirect:116` / `FillPhi:130` validate `sa[r]` and `ResolveLong:159` / `WriteU*` read it again: if the caller mutated the array concurrently the second read could be out of range (contract: not allowed; every current caller passes the freshly built array). **Fix:** document the contract on `LcpKernel.Build(PreparedText, ReadOnlySpan<int>)`; optionally copy the SA into the arena.

### MIN-4 The fallback of 1o.3b costs memory
`TryDirect` allocates `longRanks` (`n/4 + 1` ints), `escapeRanks` and `escapeValues` (`n/16 + 1` ints each) and a `n + 1` byte LCP array before the scan (`Direct.cs:57-62`); on give-up all of it is garbage and the three-pass kernel allocates its own. Measured peak RSS of the shipped chain: embedded repeats 16.6-17.0 B/char (`U8` result), `dup400-p20/p25` 18.6-18.8, `embedded-20pct` (`U32`) **20.2 B/char**, against 16.4 for the same `U16` kind when the attempt aborts after 1,024 ranks and leaves its pages untouched (`rep-100x-1pct`) and 18.4 for `U32` texts that give up immediately (`a^n`). All below the 24 B/char gate; the Log's 12.4 B/char is the random-DNA figure. **Fix:** grow `longRanks` lazily, or allocate the output array after the scan. **Resolution (1o.3c):** the heavy-text give-up allocates nothing that is touched (the lists are virtual, the byte array is touched for the first 1,024 ranks only); texts with a repeat no longer fall back, their peak RSS is 12.8-13.0 B/char (embedded repeats up to 5 %), 15.9-16.5 B/char (dup400-p20/p25), 18.0 B/char (the `U32` row: 40 MB of result), against 16.6-20.2 before.

### MIN-5 Headline ratios are single-corpus
Cold fresh-process ESA / libsais (SA+PLCP+LCP): 1.23x (random DNA), 1.24x (`dna_10m`), 1.35x (15 kb embedded repeat), 1.47x (protein), then 1.5-2.1x for repeat-rich texts and 2.47x for `a^n` (section 5 table; on `a^n` the SA phase alone, 0.337 s, is slower than libsais SA+PLCP+LCP, 0.175 s). This is not a defect: libsais is the stopping reference (`<= 1.1x`) and the Log states `dna_10m`; but any claim about "the" ratio should name the corpus.

### MIN-6 The recipe for reduced-instruction-set runs silently tests nothing on .NET 10
`ReportInstructionSets` (explicit, `Category=ReviewIntrinsics`) prints what a configuration really selects (`analysis/review_opt/`, section 4): `DOTNET_EnableSSE41=0`, `EnableBMI1/BMI2/LZCNT/POPCNT=0`, `EnableSSE=0` leave `Sse41`, `Bmi1`, `Lzcnt`, `Popcnt` all `True`; they are ignored. Effective: `EnableHWIntrinsic=0` (everything off, `Vector128.IsHardwareAccelerated` false, software vectors, no prefetch), `EnableAVX2=0` or `EnableAVX=0` (AVX2, BMI1, LZCNT off), `EnableSSE42=0` (SSE4.1, SSE4.2, POPCNT, AVX2, BMI1, LZCNT off). I ran all of them (section 4), so the claim "fallback paths pass" holds for the three effective configurations, not for the three no-op ones. Anyone repeating the run must use the effective knobs.

### MIN-7 Portability claims are untested; one stale comment
ARM64 `Prefetch` is a deliberate no-op (the SVE intrinsics are unverified), 32-bit processes and big-endian are not exercised here (`CommonPrefix` has a `BitConverter.IsLittleEndian` branch; the `Vector128` code is endian-neutral). `LcpKernel.cs:119` says the hot loops are "marked AggressiveOptimization"; `grep` finds no such attribute, and `jit_check.sh` shows `Tier1-OSR` listings for `ScanDirect`, `PlcpLoop`, `FillPhi`, `WriteU16`: with 0 range-check helpers and 0 accessor calls (reproduced here). The comment is stale, the code is what was measured.

### MIN-8 An invalid suffix array is processed safely and in linear time, but silently
`InvalidSuffixArrays_InRange_AreProcessedInLinearTime_ByBothLcpKernels` (new): a random permutation and an all-zero array over `a^n` and a period-3 text, n = 400,000, both LCP kernels: no exception, no out-of-range access (guard by `(uint)s >= (uint)n`), < 5 ms each. The Kasai loop is O(n) whatever the entries are (l <= n and drops by 1 per step), the direct scan gives up through its caps. So `Build(PreparedText, sa)` is not a validator; S1b must validate an SA read from a container itself (`SuffixArrayVerifier` is test-only).

## 3. Memory-safety audit (item 1)

Method: read every unsafe access of the five iterations, state the index-range invariant that keeps it inside its buffer, then try to break each invariant with mutants (section 8) and with inputs the corpora do not contain (section 9). All paths are in `src/SuffixTree/Algorithms/SuffixTree.Esa`. No kernel dereferences a raw pointer: `grep` for `.Pointer`, `Unsafe.` and `stackalloc` finds only the accessor itself (`Memory/Buf.cs:86,134`) and the stack maps of the reference-level `TextPreparation`. Every kernel access therefore goes through `Buf<T>`, i.e. it is bounds-checked in the `-p:EsaBounds=true` build.

### 3.1 Invariants per kernel

| kernel (file:lines) | accesses | invariant that keeps them inside the buffer | where it is established / evidence |
|---|---|---|---|
| `CountSymbols` (`SaisKernel.cs:198-208`) | `t[0..n-1]`, `counts[c]` | `0 <= c < sigma`; `counts` has `sigma` entries (`:109`) | explicit `(uint)c >= (uint)sigma` check at `:204` (`ArgumentException`); every later bucket index `c` is a symbol that passed it. Mutants S20/M62 (`>`), R03 killed. |
| `ClassifyTypes` (`:211-233`) | `t[n-1]`, `t[i]` for `i <= n-2`; sink writes `0..n` | `n >= 1` (`Sais` returns for `n == 0` at `:101`); type block / packed block has `n + 1` bytes (`:117`, `:124`); every index `0..n` is written, so the `zero: false` blocks are fully initialised | `n = 0, 1, 2` in the default corpus; poison sweep, section 9. R30/R31 (hostile sentinel byte `0xFF`) survive the whole suite: `types[n]` / `packed[n]` is never read. |
| bucket setup (`:235-253`) | `counts[c]`, `bkt[c]`, `c < sigma` | `bkt` has `sigma` entries (`:138`); all `sigma` entries are rewritten before use (`zero: false` is safe) | R04 (tail of the largest symbol left stale) is output-equivalent: the largest symbol is always L-type (a run of it ends at a smaller symbol or at the text end), so its tail is never used; hostile value `int.MinValue` passes 8.4 M + 2.4 M + 5.6 M exhaustive strings (E04h). The state oracle kills the mutant (9 tests): it compares the arrays. |
| `PlaceLmsSuffixes` (`:259-274`) | `sa[slot & -lms]`, `bkt[c]` | for an LMS position the slot lies in `[head(c), tail(c))` because the LMS count of a bucket is at most its S-type count, at most `counts[c]`; a non-LMS position stores to the dummy index 0, restored by `sa[0] = n` (`:273`) | exhaustive tests of the kernel in a guarded arena; checked build; R05 (`tt.Get(n-2)`, reads index -1 when `n = 1`) is invisible in Release and killed only by the checked build |
| `CompactSortedLms` (`:277-293`) | `sa[i]`, `i in [1,n]`; `tt.IsS(pc)`, `tt.IsS(before)`; write `sa[1+k]` | after round 1 every slot `1..n` holds a position `0..n-1` (all n suffixes are placed), so `pc in [0,n-1]`, `before in [0,n-1]`; the write index `1+k <= i` (behind the read cursor) | the clamp `p & -positive` is dead code: E06h (removed, and the pass THROWS if `p < 0` is ever read) passes the exhaustive tests; R06 survives. `k != m` is re-checked (`:291`). |
| `CompactNames` (`:299-308`) | read `sa[i]`, `i in (m,n]`; write `sa[j]`, `j >= i` | `j` starts at `n` and moves down at most once per `i`, so `j >= i > m`; `j >= n - m >= m + 1` | C13/C14 (cursor always moves / skips slot `1+m`) and S10/N08 of the earlier audits are killed; the stray entry at `n - m` lies outside the tail (S11 argument of the earlier audit) |
| `SolveReduced` (`:314-351`) | `sa1 = sa[0..m]`, `s1 = sa[n+1-m..n]`, `sa1[1+s1[i]]`, `sa[tailBase+sa[i]]` | `m <= (n-1)/2` (LMS positions lie in `[1,n-2]`, pairwise >= 2 apart), so `m + 1 <= n + 1 - m`: head and tail are disjoint; `s1[i] < names <= m`; reduced ranks `sa[i] in [0, m-1]` | R07 (tail one slot lower: still disjoint because of the spare slot) is output-equivalent: 16 M exhaustive strings (E07). R13 (`sa1` one entry short) killed by the checked build only. |
| `PlaceSortedLms` (`:354-366`) | `sa[r]`, `bkt[c]`, `sa[slot]` | `slot = bkt[c] - 1 >= r` (largest first), `slot <= n` | killed mutants N10/N16/N18 of the earlier audit; R10 killed |
| `InduceL` / `InduceS` (`SaisInduction.cs:26-66`) | `sa[i]`, look-ahead `sa[min(i+32,n)]` / `sa[max(i-32,1)]`, `tt.Get(jc)`, `bkt[c]`, `sa[slot & -take]` | `i` in `[0,n]` / `[1,n]`; look-ahead index clamped into the same range (`i + 32 <= n + 32 < 2^31` because `n <= Array.MaxLength - 1`); `jc = (p-1) & -positive in [0,n-1]`; cursor stays in its bucket (counts equal the number of induced entries); `sa[0]` is only a dummy target and is restored (`:44`, `:65`) | P01-P03, R41 (distance `int.MaxValue`: `i + d` overflows) killed; clamps now pinned by `ReviewS1aOptPrefetchTests` (section 8) |
| `NameLmsSubstrings` (`:71-92`) | look-ahead `sa[min(i+16,m)]`, `sa[1+m+(pos>>1)]`, `tt.Same(x,y)` with `x,y < n`, `tt.IsS(x-1)` with `x > a >= 1` | `pos <= n-2`, so the name slot `1+m+(pos>>1) <= 1 + (n-1)/2 + (n-2)/2 = n` | R22 (`x != b`), R21 killed; R44 (negative distance) killed by the checked build only |
| `FillPhi` (`LcpKernel.cs:124-135`) | `sa[1..n]`, `phi[s]` | `(uint)s >= (uint)n` is checked for every entry (`:130`), `phi` has `n` entries zero-initialised (`:102`) | `FillPhi`'s return value is turned into `ArgumentException`; arbitrary permutations and repeated entries are tested (`LcpKernelTests`) |
| `PlcpLoop` / `CommonPrefix` (`:138-166`, `:269-292`) | `phi[j]`, look-ahead `phi[min(j+16,n-1)]`, text reads at `j+l .. j+l+max`, `p+l .. p+l+max` | `limit = n - max(j,p)`, the call is made only for `l < limit` with `max = limit - l`, so `a+max <= j + limit <= n` and `b+max <= n`; word reads need `i + perWord <= max`; head loop is capped by `max` | R52 (head loop ignores `max`), R54 (word loop one word too far) killed (7 and 3 tests); R50 (`l <= limit`) is equivalent (`max = 0` returns 0) |
| `WriteU8/U16/U32` (`:204-256`) | `plcp[sa[r]]`, `dst[r]`, `ranks[k]`, `values[k]` | `sa[r] in [0,n)` (validated in pass 1 and not modified), `dst` has `n + 1` entries, `k < escapes` checked before each store (`:228`, `:250`), `k == escapes` at the end | the escape count comes from the same loop that produced PLCP; R53 (wrong announced count) killed |
| `ScanDirect` (`LcpKernel.Direct.cs:101-137`) | `sa[r]`, look-ahead `sa[min(r+32,n)]`, vector loads at `previous` and `s`, `longRanks[longs++]`, `dst[r]` | `(uint)s >= (uint)n` is checked (`:116`) before `s` is used; a vector is loaded only if `n - max(previous,s) >= width` (`:118-122`), so both 16-byte loads end at or before `n`; `longs < longCap` is checked before the store; `previous = n` only for `r = 1` (then `limit = 0`) | D05/D06/D14/R69/R70 killed (checked build or crash); the look-ahead value is never used as an index into memory (prefetch only) |
| `ResolveLong` (`:143-186`) | `longRanks[k+8]`, `sa[ahead-1]`, `sa[ahead]`, `sa[r-1]`, `sa[r]`, vector loads at `a+i`, `b+i`, scalar tail | long ranks have `r >= 2` (`r = 1` has `limit = 0`), both suffixes have `>= width` symbols left; the vector loop requires `i + width <= limit`, the tail `i < limit`; `escapes == escapeCap` is checked before each store (`escapeRanks` has `escapeCap + 1` entries) | R65/R66 (no spare slot) survive Release AND checked build: the spare slot is not needed; D13/D14/D21 killed |

### 3.2 n at the int limit, empty / 1 / 2 characters, alphabets to 65,536

* **Limit.** `EsaLimits.EnsureTextLength` rejects `n >= Array.MaxLength` (`EsaLimits.cs:16`), so `n + 1 <= Array.MaxLength = 2,147,483,591` and every int expression in the kernels is below `2^31`: the largest are `n + 1` (bucket sum, SA length), `i + 32` (`<= n + 32`), `r + 32`, `1 + m + (p >> 1) <= n`, `bkt[c] + take <= n + 1`, `(long)2 * n` (the only 64-bit expression, `Direct.cs:77`). `checked((nuint)length * sizeof(T))` guards the arena (`NativeArena.cs:65`). The case cannot be run (a string of 2^31 chars needs 4 GB plus SA 8.6 GB plus text 2 GB plus the packed copy: more than 15 GB), so this is argument plus the boundary tests of the guard itself (`IndexWidthTests`), not a run. The branch-free maximum `top = previous ^ ((previous ^ s) & ((previous - s) >> 31))` is correct for all `0 <= previous,s <= n < 2^31` (the difference cannot overflow).
* **n = 0, 1, 2.** `Sais` returns before the first allocation for `n = 0` (`:101`); for `n = 1` there is no LMS position (`m = 0`), `sa1[0] = m` and the loops are empty; `TryDirect` allocates `longCap + 1`, `escapeCap + 1` entries so `n < 4` works; all in the default corpus (every string over 2 symbols up to length 10, `n = 0..5` of every family) and in the sweep of section 9 in both build modes. `fixed` over an empty span yields a null pointer: it is never dereferenced (n = 0 has no loop iteration, and `Prefetch` is not called).
* **Alphabets.** `Width` is `Byte` up to 256 symbols, `UShort` up to 65,536 (`TextPreparation.cs:48`); packed symbol+type bytes are used for `sigma < 128` (`(symbol << 1) | isS <= 253`); `sigma = 128` would still fit (R08 is equivalent) but is not used; 129..256 use the split representation; `ushort` texts use `counts` of up to 65,536 ints; recursion levels use `IntPtrText` with symbols `< m < 2^30` (`d = c - next` cannot overflow). The sweep covers `sigma` in {1,2,3,4,5,16,20,126,127,128,129,200,255,256,257,1000,4096,65535,65536} at sizes 0..1025 (default test) and 4095..131073 (explicit), five text shapes each, through every optimized engine in a guarded arena with the checked and the raw build.

### 3.3 Reads of memory that was never written

* Workspace blocks with `zero: false`: `packed`/`types` (every index `0..n` written by `ClassifyTypes`), `bkt` (every entry written by `SetBucketHeads/Tails` before use), `longRanks`/`escapeRanks`/`escapeValues` (written before read: reads are bounded by the counters), `sa` (n + 1 entries: `sa[0] = n`, `Fill(Empty)` of `1..n` before the first use, `:138-142`). The registry harness fills these blocks with four poison bytes; the review test `EveryPoisonByte_GivesTheReferenceResult_ForAllOptimizedEngines` runs all 256 byte values over 47 recursion-heavy texts and all three optimized engines: green (a 4-pattern poison can hide a read-before-write whose value is inert for all four: `0x00` is "position 0", `0xCD`/`0xFF` are "empty").
* **Managed output arrays are NOT poisoned by the harness.** `GC.AllocateUninitializedArray` is used for the suffix array (`SaisKernel.cs:57`), the compact LCP arrays (`LcpKernel.cs:184,190,196`, `Direct.cs:62`). Above ~2 KB they return stale heap memory (large-object heap: never cleared), but a fresh process often hands out zero pages, so an unwritten cell (mutant R60: `dst[0]` never written) is invisible on a clean heap. `ProductionEntryPoints_DoNotDependOnTheInitialContentOfTheirOutputArrays` dirties the heap with 0xFF-filled arrays of the same sizes, forces compacting collections and rebuilds: it kills R60 (one pre-existing test kills R60 as well; with the added tests 3 tests do).
* **Stale-value dependence that exists by design** (the earlier audit's S06/S08/S11/S13/S17/S19/N06 and mine R04/R06/R07/R30/R31): all proven output-equivalent with hostile variants through the exhaustive tests (section 8).

### 3.4 Aliasing, lifetime, thread safety

* **Aliasing.** The only aliasing is inside the SA buffer (head `sa1 = [0..m]`, tail `s1 = [n+1-m..n]`, disjoint because `2m + 1 <= n`); text, `sa`, `phi`, `counts`, `bkt`, `packed` are distinct allocations; `LcpKernel.Build` pins the caller's `ReadOnlySpan<int>` and never writes it (the registry checks that the input is unchanged after the call).
* **Lifetime.** Every public entry point creates its `NativeArena` in a `using` (`SaisKernel.cs:46`, `LcpKernel.cs:40,51`); each recursion level `Mark()`s and `ReleaseTo()`s (`SaisKernel.cs:107,129`; N12/C20 killed by `TheWorkspaceOfEveryLevel_IsReleased_WhenBuildReturns`); exceptions leave blocks in the arena, freed by `Dispose` (and by the finalizer if a caller forgets). `fixed` blocks cover every use of a pointer into a managed array (`SaisKernel.cs:58`, `LcpKernel.cs:63,175-176,185,191,197`, `Direct.cs:66`); `BuildsUnderGcPressure_WithForcedCompactingCollections` (new) runs the kernels while another thread forces compacting collections: green. `NativeArena.Dispose` is idempotent; a `Buf<T>` outliving its arena is a use-after-free that the type system does not prevent (internal API: all uses are inside the call that created the arena).
* **Thread safety.** No mutable static state: `grep` finds no static field in the kernels (`SaisPrefetch`, `LcpPrefetch`, `EsaLimits` hold `const`s). `ConcurrentBuilds_OfDifferentTexts_GiveTheReferenceResults` (new): 12 builds of four different texts on four threads, all equal to the reference. A `NativeArena` is documented "not thread-safe" and is never shared.
* **TOCTOU.** `LcpKernel.Build(PreparedText, ReadOnlySpan<int> sa)` validates `sa[r]` once (`FillPhi:130`, `ScanDirect:116`) and reads it again later (`WriteU*`, `ResolveLong:159`). If another thread modified the caller's array between the two reads, the second read could be out of range. Contract: the caller must not mutate `sa` during the call (it is the kernel's own freshly built array in every current use). See MIN-3.

### 3.5 Prefetch address computations

`Buf<T>.Prefetch` is `Sse.Prefetch0(_p + index)` (`Buf.cs:149-153`): pointer arithmetic on an unvalidated `int`, no dereference, never faults on x86 (the hardware drops a prefetch of an unmapped address), no-op elsewhere. Overflow: `_p + index` is computed in `nint` (sign-extended), no wrap on 64-bit. The *look-ahead reads* that feed the prefetch are real loads and are what the clamps protect: `sa[Math.Min(i + 32, n)]` (`SaisInduction.cs:31`), `sa[Math.Max(i - 32, 1)]` (`:52`), `sa[Math.Min(i + 16, m)]` (`:79`), `phi[Math.Min(j + 16, n - 1)]` (`LcpKernel.cs:146`), `sa[Math.Min(r + 16, n)]` (`:209,220,242`), `sa[Math.Min(r + 32, n)]` (`Direct.cs:113`), `longRanks[k + 8]` guarded by `k + 8 < longs` (`:152`). Mutants that remove a clamp or make a distance overflow (P01-P03, Q01-Q05, D21, D22, R41, R42, R44, R57, R58, R69, R70) are killed (Release crash, or checked build, or the new pin test); the three that changed only the prefetch ADDRESS (P04, P05, P07, R20, R24, R40, R43) used to survive every test and are now killed by `ReviewS1aOptPrefetchTests`, which pins the target of every prefetch to the position read exactly `Induce` / `Name` iterations later.

### 3.6 Does the checked build / canary net cover each kernel?

* Checked build: yes for every kernel access, because no kernel bypasses `Buf<T>`. Not covered by design: `Prefetch` (never faults), `Fill/Clear` call `AsSpan(count)` which is checked in the checked build. In Release `Slice(start)` returns a Buf without a length: nothing to check.
* Canaries (256-byte bands, position-dependent pattern) cover writes into arena blocks: workspace, and, in the registry, also the SA and the text copy. They do **not** cover the managed output arrays (`sa` of the production entry, `lcp`, `overflow*`): a one-past write there corrupts the managed heap and is detected only by the checked build. All OOB mutants of the campaign that survived Release (R05, R13, R44; earlier K14/K15/P01-P03/Q01-Q05/D05/D14/D22) were caught by the checked build, so the combination works, but only if the checked run is part of every acceptance (`.github/workflows/dotnet.yml:53-56` runs only `dotnet test Seqeron.sln -c Release`: MIN-1).
* The safety net itself was mutated (R80-R86, R90-R94): see the mutant table and MAJ-2: the upper and the lower bound of `ReadUnaligned` and of `LoadVector128` could be weakened by one element (or removed) without any pre-existing test, even in the checked build, noticing, until `ReviewS1aOptTests` pinned the exact bounds.


## 4. Fallback paths (item 2)

Hardware-intrinsic uses: `Sse.Prefetch0` behind `Sse.IsSupported` (`Buf.cs:149-153`), `Vector128.Load/Equals/ExtractMostSignificantBits` (`Buf.cs:91,138`, `Direct.cs:193`), `BitOperations.TrailingZeroCount/LeadingZeroCount` (`LcpKernel.cs:285`, `Direct.cs:194`). The full Esa test project (default, `Category=Exhaustive` with 24 tests including 8.4 M binary, 2.4 M ternary and 5.6 M 4-ary strings through each engine, `Category=Vectors` with 474 libsais-hash tests at up to 10 Mbp in all five pipelines, `Category=Fuzz` with 18 tests, seed 7, scale 0.05) was run in each configuration on a binary of commit `b20f6609` (267 default tests); the final 273-test binary was run again in the default category (raw, and checked for the base configuration) at the end. Evidence of what each configuration really selects (`ReportInstructionSets`, output in `analysis/review_opt/intrinsics_report.txt`):

| configuration | what the runtime reports | default | Exhaustive | Vectors | Fuzz |
|---|---|---|---|---|---|
| `DOTNET_EnableHWIntrinsic=0` | **everything off**: `Sse.IsSupported` false (no prefetch), `Vector128` software | 267 pass (35 s) | 24 pass (7 m 59 s) | 474 pass | 18 pass |
| `DOTNET_EnableAVX2=0` | Avx2, Bmi1, Lzcnt off; Sse41, Popcnt on | 267 pass | 24 pass (6 m 17 s) | 474 pass | 18 pass |
| `DOTNET_EnableSSE42=0` | **Sse41, Popcnt, Avx2, Bmi1, Lzcnt off**, Sse and Vector128 on | 267 pass | 24 pass (5 m 53 s) | 474 pass | 18 pass |
| `DOTNET_TieredPGO=0` | all ISAs on, no dynamic PGO | 267 pass | 24 pass (6 m 39 s) | 474 pass | 18 pass |
| `DOTNET_TieredCompilation=0` | all ISAs on, everything fully optimised at first call | 267 pass | 24 pass (6 m 31 s) | 474 pass | 18 pass |
| `DOTNET_EnableSSE41=0` (**no-op**) | identical to base | 267 pass | 24 pass | 474 pass | 18 pass |
| `DOTNET_EnableBMI1=0 EnableBMI2=0 EnableLZCNT=0 EnablePOPCNT=0` (**no-op**) | identical to base | 267 pass | 24 pass | 474 pass | 18 pass |
| `DOTNET_EnableSSE=0` (**no-op**) | identical to base | 267 pass | 24 pass | 474 pass | 18 pass |

The final binary (commit `b20f6609` plus this review's 273 tests) was run again after the last test change: default (273 tests) and `Category=ReviewSweep` pass in the base, `EnableHWIntrinsic=0`, `EnableAVX2=0`, `EnableSSE42=0`, `TieredPGO=0` and `TieredCompilation=0` configurations in the raw build, and in the checked build (`-p:EsaBounds=true`, `ESA_EXPECT_BOUNDS=true`) for base, `EnableHWIntrinsic=0` and `EnableSSE42=0`. One failure occurred on the way and is a defect of one of my own tests, not of the kernels: `BuildsUnderGcPressure_...` asserted that the second thread had completed more than five forced collections, and a run of the avx2-off configuration that overlapped a mutant run did not get them; the test now builds until 20 collections happened (bounded by 15 s) and passed 3 runs next to six busy loops.

**Result: no difference in any category, in any configuration** (the three no-op rows are additional base-configuration repetitions). The software `Vector128` path (`EnableHWIntrinsic=0`) is slower (Exhaustive 7 m 59 s against 5 m 13 s - 6 m 46 s for the others, Vectors 75 s against 50-71 s) but produces identical results, including the 474 libsais hashes (corpora up to 10 Mbp). Raw: `fallback_runs.txt`.

## 5. Worst-case complexity of the LCP in suffix-array order (item 3)

Inputs: the families of the task at 1 Mbp and 10 Mbp plus de Bruijn, protein, 16-bit alphabets, and repeat shapes designed against the give-up rules (`ReviewS1aOptWorstCaseTests.Family`): `a^n`, `abab`, Fibonacci, Thue-Morse, period 2 and period 1000 with 0.1 % noise, `rep_10m` (the repository's file) and a 100x-1 % generator, de Bruijn B(2,24) and B(4,12), `dup400-pXX` (random DNA, XX % of the 400-letter blocks copy an earlier block: XX % of the ranks are long with LCP up to 400), `embedded-X` (one block of X % of the length occurs twice), two and three copies of a random block, tandem period 1000 (exact), `akb-k100`, runs of length 256, random protein, 16-bit alphabets (`u16-random-1000`, `u16-bmp` with 49,152 symbols, `u16-all-same`, `u16-period1000-noise`). Two instruments: **cold** = fresh process of the harness (`build`, JIT included, median of 3, the libsais driver on the same file, SA and LCP hashes compared with libsais for every row), and **warm** = the explicit test `ReviewS1aOptWorstCaseTests.WorstCaseTimings` (SA kernel, safe reference LCP, three-pass kernel, shipped direct kernel and the bare direct attempt in one process, best of 3).

**Cold, 10 Mbp, sorted by ESA / libsais** (every row: SA and LCP SHA-256 equal to libsais; `analysis/review_opt/cold_10m.tsv`, `cold_embedded_10m.tsv`):

| text (10 Mbp) | libsais SA+PLCP+LCP | ESA direct (shipped) | ESA three-pass | ESA / libsais | direct / three-pass | LCP phase: direct / three-pass / safe ref | LCP kind | peak RSS B/char | SA+LCP = libsais |
|---|---|---|---|---|---|---|---|---|---|
| random-dna | 0.605 | 0.746 | 0.951 | 1.23 | 0.78 | 0.073 / 0.277 / 0.437 | U8 | 12.3 | yes |
| dna_10m | 0.635 | 0.788 | 0.990 | 1.24 | 0.80 | 0.079 / 0.281 / 0.419 | U8 | 12.2 | yes |
| embedded-015pct | 0.614 | 0.831 | 1.049 | 1.35 | 0.79 | 0.113 / 0.311 / 0.387 | U8 | 12.8 | yes |
| random-protein | 0.629 | 0.927 | 1.078 | 1.47 | 0.86 | 0.068 / 0.244 / 0.337 | U8 | 13.1 | yes |
| dup400-p10 | 0.565 | 0.842 | 0.940 | 1.49 | 0.90 | 0.170 / 0.261 / 0.380 | U8 | 12.6 | yes |
| debruijn-4 | 0.318 | 0.495 | 0.613 | 1.56 | 0.81 | 0.044 / 0.140 / 0.180 | U8 | 12.3 | yes |
| random-binary | 0.566 | 0.887 | 0.869 | 1.57 | 1.02 | 0.300 / 0.286 / 0.369 | U8 | 15.4 | yes |
| fibonacci | 0.480 | 0.791 | 0.820 | 1.65 | 0.96 | 0.192 / 0.192 / 0.230 | U32 | 18.4 | yes |
| dup400-p15 | 0.565 | 0.937 | 1.019 | 1.66 | 0.92 | 0.248 / 0.287 / 0.348 | U8 | 12.9 | yes |
| period1000-noise | 0.463 | 0.795 | 0.837 | 1.72 | 0.95 | 0.187 / 0.184 / 0.222 | U16 | 16.4 | yes |
| thue-morse | 0.580 | 1.007 | 0.914 | 1.74 | 1.10 | 0.304 / 0.266 / 0.417 | U32 | 18.4 | yes |
| two-copies | 0.554 | 0.975 | 0.970 | 1.76 | 1.01 | 0.247 / 0.237 / 0.296 | U32 | 18.4 | yes |
| three-copies | 0.523 | 0.940 | 0.914 | 1.80 | 1.03 | 0.216 / 0.210 / 0.292 | U32 | 18.4 | yes |
| runs-L256 | 0.351 | 0.631 | 0.637 | 1.80 | 0.99 | 0.180 / 0.184 / 0.223 | U16 | 16.4 | yes |
| debruijn-2 | 0.330 | 0.602 | 0.608 | 1.82 | 0.99 | 0.177 / 0.166 / 0.214 | U8 | 15.4 | yes |
| rep_10m | 0.516 | 0.957 | 0.935 | 1.85 | 1.02 | 0.237 / 0.225 / 0.335 | U16 | 16.5 | yes |
| tandem-1000 | 0.422 | 0.815 | 0.831 | 1.93 | 0.98 | 0.189 / 0.179 / 0.241 | U32 | 18.4 | yes |
| embedded-040pct | 0.609 | 1.177 | 0.995 | 1.93 | 1.18 | 0.470 / 0.286 / 0.376 | U8 | 16.6 | yes |
| rep-100x-1pct | 0.481 | 0.936 | 0.897 | 1.95 | 1.04 | 0.230 / 0.226 / 0.346 | U16 | 16.4 | yes |
| embedded-1pct | 0.600 | 1.179 | 1.007 | 1.97 | 1.17 | 0.464 / 0.296 / 0.363 | U8 | 16.6 | yes |
| akb-k100 | 0.293 | 0.584 | 0.537 | 1.99 | 1.09 | 0.191 / 0.164 / 0.286 | U32 | 18.4 | yes |
| embedded-20pct | 0.604 | 1.221 | 1.091 | 2.02 | 1.12 | 0.494 / 0.311 / 0.410 | U32 | 20.2 | yes |
| embedded-5pct | 0.630 | 1.301 | 1.083 | 2.07 | 1.20 | 0.506 / 0.315 / 0.459 | U8 | 17.0 | yes |
| period2-noise | 0.356 | 0.742 | 0.730 | 2.08 | 1.02 | 0.150 / 0.145 / 0.176 | U16 | 16.4 | yes |
| dup400-p25 | 0.583 | 1.238 | 1.004 | 2.12 | 1.23 | 0.568 / 0.271 / 0.411 | U16 | 18.8 | yes |
| ab-repeat | 0.253 | 0.541 | 0.513 | 2.14 | 1.05 | 0.086 / 0.078 / 0.103 | U32 | 18.4 | yes |
| dup400-p20 | 0.564 | 1.210 | 0.911 | 2.15 | 1.33 | 0.518 / 0.249 / 0.336 | U16 | 18.6 | yes |
| all-same | 0.175 | 0.432 | 0.426 | 2.47 | 1.01 | 0.095 / 0.091 / 0.109 | U32 | 18.4 | yes |


**Growth with n (1 Mbp -> 10 Mbp, same instrument, cold).** The ratio of times for a 10x larger text is 6.9-13.5x for the shipped chain (26 rows) and 9.4-18.5x for libsais: no text grows faster than the memory hierarchy explains, and the shipped chain grows no faster than libsais on any of the 26 rows (`analysis/review_opt/cold_1m.tsv`). Selected rows:

| text | libsais 1 Mbp -> 10 Mbp | ESA direct 1 Mbp -> 10 Mbp | growth libsais | growth ESA direct | growth ESA three-pass |
|---|---|---|---|---|---|
| random-dna | 0.042 -> 0.605 | 0.073 -> 0.746 | 14.4x | 10.2x | 10.8x |
| embedded-015pct | 0.043 -> 0.614 | 0.079 -> 0.831 | 14.3x | 10.5x | 11.2x |
| random-protein | 0.041 -> 0.629 | 0.073 -> 0.927 | 15.3x | 12.7x | 12.1x |
| dup400-p10 | 0.041 -> 0.565 | 0.085 -> 0.842 | 13.8x | 9.9x | 10.6x |
| debruijn-4 | 0.024 -> 0.318 | 0.072 -> 0.495 | 13.2x | 6.9x | 7.9x |
| random-binary | 0.042 -> 0.566 | 0.091 -> 0.887 | 13.5x | 9.7x | 9.8x |
| fibonacci | 0.028 -> 0.480 | 0.093 -> 0.791 | 17.1x | 8.5x | 9.2x |
| dup400-p15 | 0.040 -> 0.565 | 0.092 -> 0.937 | 14.1x | 10.2x | 10.8x |
| period1000-noise | 0.026 -> 0.463 | 0.084 -> 0.795 | 17.8x | 9.5x | 10.1x |
| thue-morse | 0.036 -> 0.580 | 0.097 -> 1.007 | 16.1x | 10.4x | 9.6x |
| two-copies | 0.035 -> 0.554 | 0.093 -> 0.975 | 15.8x | 10.5x | 11.0x |
| three-copies | 0.032 -> 0.523 | 0.089 -> 0.940 | 16.3x | 10.6x | 10.5x |
| runs-L256 | 0.019 -> 0.351 | 0.074 -> 0.631 | 18.5x | 8.5x | 9.0x |
| debruijn-2 | 0.035 -> 0.330 | 0.084 -> 0.602 | 9.4x | 7.2x | 7.5x |
| tandem-1000 | 0.024 -> 0.422 | 0.086 -> 0.815 | 17.6x | 9.5x | 10.5x |
| embedded-040pct | 0.042 -> 0.609 | 0.087 -> 1.177 | 14.5x | 13.5x | 10.3x |
| rep-100x-1pct | 0.028 -> 0.481 | 0.084 -> 0.936 | 17.2x | 11.1x | 10.8x |
| embedded-1pct | 0.043 -> 0.600 | 0.111 -> 1.179 | 14.0x | 10.6x | 10.7x |
| akb-k100 | 0.017 -> 0.293 | 0.077 -> 0.584 | 17.2x | 7.6x | 7.1x |
| embedded-20pct | 0.043 -> 0.604 | 0.115 -> 1.221 | 14.0x | 10.6x | 11.2x |
| embedded-5pct | 0.045 -> 0.630 | 0.115 -> 1.301 | 14.0x | 11.3x | 11.2x |
| period2-noise | 0.025 -> 0.356 | 0.085 -> 0.742 | 14.2x | 8.7x | 8.9x |
| dup400-p25 | 0.039 -> 0.583 | 0.113 -> 1.238 | 14.9x | 11.0x | 11.4x |
| ab-repeat | 0.021 -> 0.253 | 0.076 -> 0.541 | 12.0x | 7.1x | 6.8x |
| dup400-p20 | 0.040 -> 0.564 | 0.111 -> 1.210 | 14.1x | 10.9x | 10.5x |
| all-same | 0.017 -> 0.175 | 0.059 -> 0.432 | 10.3x | 7.3x | 7.5x |


**Warm (one process), 10 Mbp.** Phase times and whether the bare attempt succeeded (`worst_warm_10m.tsv`, `worst_warm_embedded_10m.tsv`; the last column is warm LCP growth 1 -> 10 Mbp):

| text (10 Mbp) | sigma | SA-IS pointer | LCP safe ref | LCP three-pass | LCP direct | bare direct attempt | attempt succeeded | kind | direct / three-pass | direct / ref | LCP direct growth 1->10 Mbp |
|---|---|---|---|---|---|---|---|---|---|---|---|
| random-dna | 4 | 0.613 | 0.337 | 0.243 | 0.061 | 0.087 | yes | U8 | 0.25 | 0.18 | 30.6x |
| random-binary | 2 | 0.489 | 0.354 | 0.266 | 0.276 | 0.025 | no | U8 | 1.04 | 0.78 | 16.8x |
| random-protein | 20 | 0.822 | 0.338 | 0.216 | 0.060 | 0.082 | yes | U8 | 0.28 | 0.18 | 29.8x |
| all-same | 1 | 0.279 | 0.065 | 0.077 | 0.079 | 0.024 | no | U32 | 1.03 | 1.22 | 24.0x |
| ab-repeat | 2 | 0.351 | 0.063 | 0.060 | 0.066 | 0.021 | no | U32 | 1.11 | 1.05 | 20.8x |
| fibonacci | 2 | 0.563 | 0.195 | 0.167 | 0.181 | 0.028 | no | U32 | 1.08 | 0.93 | 17.0x |
| thue-morse | 2 | 0.622 | 0.383 | 0.256 | 0.275 | 0.026 | no | U32 | 1.08 | 0.72 | 16.0x |
| period2-noise | 4 | 0.406 | 0.192 | 0.131 | 0.129 | 0.023 | no | U16 | 0.98 | 0.67 | 26.3x |
| period1000-noise | 4 | 0.516 | 0.201 | 0.165 | 0.169 | 0.023 | no | U16 | 1.02 | 0.84 | 19.6x |
| rep-100x-1pct | 4 | 0.560 | 0.254 | 0.194 | 0.202 | 0.025 | no | U16 | 1.04 | 0.79 | 21.2x |
| debruijn-2 | 2 | 0.334 | 0.195 | 0.163 | 0.162 | 0.020 | no | U8 | 0.99 | 0.83 | 15.3x |
| debruijn-4 | 4 | 0.385 | 0.163 | 0.136 | 0.037 | 0.062 | yes | U8 | 0.27 | 0.22 | 19.3x |
| dup400-p10 | 4 | 0.616 | 0.331 | 0.253 | 0.165 | 0.192 | yes | U8 | 0.65 | 0.50 | 15.5x |
| dup400-p15 | 4 | 0.612 | 0.332 | 0.243 | 0.241 | 0.256 | yes | U8 | 0.99 | 0.73 | 15.4x |
| dup400-p20 | 4 | 0.603 | 0.346 | 0.260 | 0.545 | 0.310 | no | U16 | 2.10 | 1.58 | 17.8x |
| dup400-p25 | 4 | 0.655 | 0.345 | 0.248 | 0.532 | 0.326 | no | U16 | 2.15 | 1.54 | 17.3x |
| two-copies | 4 | 0.652 | 0.280 | 0.201 | 0.207 | 0.025 | no | U32 | 1.03 | 0.74 | 17.2x |
| three-copies | 4 | 0.601 | 0.236 | 0.186 | 0.195 | 0.026 | no | U32 | 1.04 | 0.82 | 18.7x |
| tandem-1000 | 4 | 0.543 | 0.177 | 0.151 | 0.164 | 0.024 | no | U32 | 1.08 | 0.92 | 21.6x |
| u16-random-1000 | 1000 | 1.039 | 0.286 | 0.267 | 0.069 | 0.091 | yes | U8 | 0.26 | 0.24 | 32.8x |
| u16-bmp | 49152 | 0.800 | 0.271 | 0.223 | 0.069 | 0.093 | yes | U8 | 0.31 | 0.25 | 32.8x |
| u16-all-same | 301 | 0.268 | 0.072 | 0.079 | 0.080 | 0.023 | no | U32 | 1.01 | 1.10 | 25.6x |
| u16-period1000-noise | 1000 | 0.676 | 0.183 | 0.158 | 0.175 | 0.025 | no | U16 | 1.11 | 0.96 | 19.9x |
| akb-k100 | 2 | 0.303 | 0.218 | 0.139 | 0.151 | 0.020 | no | U32 | 1.08 | 0.69 | 26.5x |
| runs-L256 | 2 | 0.277 | 0.182 | 0.150 | 0.151 | 0.019 | no | U16 | 1.01 | 0.83 | 21.6x |
| embedded-015pct | 4 | 0.609 | 0.333 | 0.287 | 0.108 | 0.137 | yes | U8 | 0.38 | 0.32 | 45.0x |
| embedded-040pct | 4 | 0.660 | 0.372 | 0.289 | 0.466 | 0.193 | no | U8 | 1.61 | 1.25 | 99.2x |
| embedded-1pct | 4 | 0.666 | 0.362 | 0.294 | 0.468 | 0.190 | no | U8 | 1.59 | 1.29 | 16.6x |
| embedded-5pct | 4 | 0.655 | 0.330 | 0.261 | 0.424 | 0.191 | no | U8 | 1.63 | 1.29 | 14.8x |
| embedded-20pct | 4 | 0.623 | 0.307 | 0.246 | 0.418 | 0.192 | no | U32 | 1.70 | 1.36 | 14.7x |


Reading: (1) **no super-linear growth** of the work: the SA-IS kernel needs 0.27-1.04 s at 10 Mbp, 11-24x the 1 Mbp time (a cache effect: the 1 Mbp working set fits in L2/L3; the 16-bit and protein alphabets are at the top); the warm LCP phase grows 15-33x (three-pass kernel 14-28x) on every row except the embedded ones, where the PATH changes with n: `embedded-040pct` (a 4,000-letter repeat at 1 Mbp is below 8 sqrt(n) = 8,000 and the attempt succeeds in 0.005 s; a 40,000-letter repeat at 10 Mbp is above 25,000 and the attempt fails: 0.47 s, ratio 99x) and `embedded-015pct` (45x: the resolution work grows with the square of the repeat length until it reaches the 2n cap). Both are the cliff of MAJ-1, not a growing constant; (2) **the direct kernel is within 0.25-1.11x of the three-pass LCP phase on every text except those that make the attempt fail late**, where it is 1.6-2.15x (embedded 0.4-20 %, `dup400-p20/p25`; whole build 1.12-1.33x); (3) **against the safe reference** the worst case is 1.5-1.6x slower in the LCP phase (`dup400-p20/p25`; the embedded rows 1.1-1.3x); (4) the SA phase on repeat-rich texts is where ESA loses to libsais (2.5x on `a^n`), not the LCP phase; (5) an invalid in-range SA is processed in linear time (MIN-8).

**Mutating the bound constants of the fallback** (`bound_constants_timing_1m.tsv`; direct-kernel LCP seconds at n = 10^6 on a loaded machine, `*` = the attempt gave up; the tests column is the mutant campaign of section 6):

| variant | tests | random-dna | random-binary | embedded-1pct | embedded-5pct | embedded-20pct | dup400-p10 | dup400-p20 | rep-100x-1pct |
|---|---|---|---|---|---|---|---|---|---|
| shipped constants (baseline) | - | 0.002 | 0.015* | 0.028* | 0.029* | 0.029* | 0.015 | 0.031* | 0.009* |
| R61: `ResolveStepsPerChar` 2 -> 1000 | killed (2 tests: 1 pre-existing, 1 new) | 0.002 | 0.015* | 0.018 | 0.401 | 1.871* | 0.011 | 0.038* | 0.010* |
| R62: `LongFractionDivisor` 4 -> 2^30 | killed (3) | 0.017* | 0.016* | 0.016* | 0.016* | 0.015* | 0.025* | 0.014* | 0.015* |
| early exit never fires (`probe = int.MaxValue`, D29 class) | survives (performance only) | 0.002 | 0.016* | 0.029* | 0.029* | 0.042* | 0.015 | 0.031* | 0.011* |
| R63: `probe = Math.Min(...)` | killed (1) | 0.002 | 0.015* | 0.028* | 0.029* | 0.029* | 0.011 | 0.032* | 0.010* |
| D01: long cap n/4 -> n/2 | killed (old suite) | 0.002 | 0.030* | 0.031* | 0.054* | 0.052* | 0.011 | 0.031* | 0.009* |
| R67-like: budget 2n -> n/16 | killed (3) | 0.002 | 0.015* | 0.018* | 0.019* | 0.019* | 0.032* | 0.019* | 0.010* |

Notices: (a) `ResolveStepsPerChar` 2 -> 1000 (R61) turns the embedded-repeat rows quadratic (0.029 -> 0.40 s at 5 %, **1.87 s at 20 %**, 65x), yet the output is identical (the fallback and the direct path give the same array): only a time-bounded test sees it; the new `Direct_OnAnEmbeddedLongRepeat_IsLinear_AndEqualsTheReference` (n = 4 million, limit 3 s, normally 0.1 s) fails on it together with one pre-existing test (killed(2)). (b) `LongFractionDivisor` -> `1 << 30` (R62: no attempt ever succeeds) and the quarter limit n/2 (D01) are killed by the pre-existing tests (3 and 1 tests) and show as an 8x loss on random DNA (0.002 -> 0.017 s). (c) The early exit (`probe`): never firing, `Min` instead of `Max` (R63) and the three-quarters variant (R64) change no row noticeably at 1 Mbp; R63 is killed by one test, R64 and "never fires" survive (they change only when the attempt gives up: MIN-2). (d) Budget `2n -> n/16` (R67) makes `dup400-p10` fall back (0.015 -> 0.032 s); killed by 3 tests.

## 6. Fresh mutant campaign (item 4)

64 new single-line mutants (ids R01-R94, `analysis/mutation/review_opt_*.json`, runner `run_kernel_mutants.py` with `ESA_ROOT` pointing at a git worktree so the production tree stayed untouched; `git diff -- src` is empty), none of them in the earlier lists. By iteration: **13 in the safety infrastructure** the kernels depend on (1o.0: the checked section of `Buf.cs`, `NativeArena.cs`), **22 in 1o.1b/1o.1c** (14 in `SaisKernel.cs`: classification, counting, buckets, placement, compaction, reduced problem, representation choice, recursion; 4 in `TypedText.cs`; 4 in the naming/induction logic of `SaisInduction.cs`), **8 in 1o.2** (prefetch distance constants, clamp of the look-ahead reads, prefetch address, wrong-side read: R20, R24, R25, R40-R44), **9 in 1o.3a** (`LcpKernel.cs`: PLCP, `CommonPrefix`, escape counting, `FillPhi`, prefetch constants R57, R58) and **12 in 1o.3b** (10 in `LcpKernel.Direct.cs`: scan, long-match resolution, caps, budget, early exit, spare slots, placeholder; R69, R70 for the direct-scan and resolution prefetch distances). Each ran against the Release (raw pointer) build with the default Kernels/Safety/Review/SaisPasses/CompactLcp tests and, when it survived, against the checked build (`-p:EsaBounds=true`, `ESA_EXPECT_BOUNDS=true`); a crashed test host counts as a kill (a corrupted pointer loop leaves its buffers). Every killed mutant, and the four whose kill depended on the checked build, was then re-run against the PRE-EXISTING suite only (Kernels, Safety, `ReviewS1aTests`, `SaisPasses`, `CompactLcp`; column "pre-existing suite only") to tell which kills are owed to the tests added by this review; the survivors of every run are marked "-".

**Result: 64 mutants. 52 killed (44 by the pre-existing suite, 8 survived it and are killed by the tests added in this review), 12 survive everything: 9 equivalent (each proven concretely below), 3 performance-only (not killable by any output test).**

| id | site | change | Release | checked | pre-existing suite only | final | note |
|---|---|---|---|---|---|---|---|
| R80 | Memory/Buf.cs:61 | checked indexer admits index == length | survived | killed(5) | survived / killed(4) | killed |  |
| R81 | Memory/Buf.cs:69 | checked Slice(start, length) admits one element too many | survived | killed(2) | survived / killed(1) | killed |  |
| R82 | Memory/Buf.cs:85 | checked ReadUnaligned admits a word that ends one byte past the buffer | survived | killed(1) | survived / survived | **survived; killed by added test** | Buf edge test |
| R83 | Memory/Buf.cs:85 | checked ReadUnaligned admits negative indices | survived | killed(1) | survived / survived | **survived; killed by added test** | Buf edge test |
| R84 | Memory/Buf.cs:93 | checked LoadVector128 admits a vector that ends one element past the buffer | survived | killed(1) | survived / survived | **survived; killed by added test** | Buf edge test |
| R85 | Memory/Buf.cs:93 | checked LoadVector128 admits negative indices | survived | killed(1) | survived / survived | **survived; killed by added test** | Buf edge test |
| R86 | Memory/Buf.cs:100 | checked AsSpan admits one element too many | survived | killed(2) | survived / killed(1) | killed |  |
| R87 | Memory/Buf.cs:152 | prefetch uses the non-temporal hint (correct results, evicts early) | survived | survived | - (survivor of every run) | **survived; performance only** | hint only (nta) |
| R92 | Memory/NativeArena.cs:102 | ReleaseTo leaves the oldest block of a level allocated and unverified | killed(2) | - | killed(2) | killed |  |
| R93 | Memory/NativeArena.cs:163 | the guard band after the payload is never verified | killed(7) | - | killed(6) | killed |  |
| R94 | Memory/NativeArena.cs:72 | zeroed blocks of a guarded arena are not cleared | killed(host crash) | - | killed(host crash) | killed |  |
| R90 | Memory/NativeArena.cs:181 | canary pattern is constant zero (a kernel that writes zeros into a guard band is not detected) | killed(2) | - | killed(1) | killed |  |
| R91 | Memory/NativeArena.cs:73 | poison does not cover the last byte of an uninitialised block (payload > 0) | killed(2) | - | killed(1) | killed |  |
| R60 | LcpKernel.Direct.cs:109 | dst[0] is not written (uninitialised managed array: invisible unless the heap is dirty) | killed(4) | - | killed(1) | killed |  |
| R61 | LcpKernel.Direct.cs:41 | step budget 1000 n instead of 2 n (quadratic worst case on a text with an embedded long repeat) | killed(2) | - | killed(1) | killed |  |
| R62 | LcpKernel.Direct.cs:32 | no cap on the long ranks (divisor 1 << 30: longCap = 0 for small n, huge otherwise) | killed(3) | - | killed(3) | killed |  |
| R63 | LcpKernel.Direct.cs:108 | early exit threshold: min instead of max (aborts after 1024 ranks even for large n: performance only) | killed(1) | - | killed(1) | killed |  |
| R64 | LcpKernel.Direct.cs:125 | early exit needs three quarters long instead of half (performance only) | survived | survived | - (survivor of every run) | **survived; performance only** | only the moment of giving up changes |
| R65 | LcpKernel.Direct.cs:59 | long ranks array has no spare slot (claimed equivalent: the cap check precedes the store) | survived | survived | - (survivor of every run) | **survived; equivalent** | checked build, exact size |
| R66 | LcpKernel.Direct.cs:60 | escape arrays have no spare slot (claimed equivalent: the cap check precedes the store) | survived | survived | - (survivor of every run) | **survived; equivalent** | checked build, exact size |
| R67 | LcpKernel.Direct.cs:172 | resolution budget counted per long rank instead of per vector step (performance only on repeats) | killed(3) | - | killed(3) | killed |  |
| R68 | LcpKernel.Direct.cs:166 | resolution: the scalar tail stops one symbol early (limit - 1) | killed(5) | - | killed(4) | killed |  |
| R72 | LcpKernel.Direct.cs:133 | the placeholder for a long rank is width + 1 (overwritten by the resolution: claimed equivalent unless the scan result is used raw) | survived | survived | - (survivor of every run) | **survived; equivalent** | placeholder always overwritten |
| R50 | LcpKernel.cs:156 | PLCP: l <= limit also calls the comparison (max = 0 returns 0: claimed equivalent) | survived | survived | - (survivor of every run) | **survived; equivalent** | max = 0 returns 0 |
| R51 | LcpKernel.cs:163 | PLCP: l is decremented without the l > 0 guard | killed(18) | - | killed(12) | killed |  |
| R52 | LcpKernel.cs:273 | CommonPrefix: the head loop ignores max (reads up to 7 symbols past the allowed range) | killed(7) | - | killed(4) | killed |  |
| R53 | LcpKernel.cs:81 | kind U16: the escape count announced from the 8-bit threshold | killed(10) | - | killed(6) | killed |  |
| R54 | LcpKernel.cs:280 | CommonPrefix: the word loop allows one word past max | killed(3) | - | killed(1) | killed |  |
| R55 | LcpKernel.cs:286 | CommonPrefix: the mismatch position of a word is counted in bits instead of symbols for bytes | killed(16) | - | killed(10) | killed |  |
| R56 | LcpKernel.cs:126 | FillPhi: the first predecessor is 0 instead of n (the sentinel) | killed(14) | - | killed(10) | killed |  |
| R57 | LcpKernel.cs:310 | pass 3 U32: the look-ahead distance overflows (int.MaxValue) | killed(host crash) | - | killed(host crash) | killed |  |
| R58 | LcpKernel.cs:307 | PLCP prefetch distance int.MaxValue (the clamp Math.Min(j + d, n - 1) sees a negative index) | killed(host crash) | - | killed(host crash) | killed |  |
| R69 | LcpKernel.cs:313 | scan: prefetch distance int.MaxValue (r + distance overflows) | killed(host crash) | - | killed(host crash) | killed |  |
| R70 | LcpKernel.cs:316 | resolve: prefetch distance negative (k + distance < longs reads long ranks before index 0) | killed(host crash) | - | killed(host crash) | killed |  |
| R01 | SaisKernel.cs:223 | classification: the distance of the symbols is reversed (next - c) | killed(host crash) | - | killed(host crash) | killed |  |
| R02 | SaisKernel.cs:226 | classification: equal symbols inherit the opposite type | killed(28) | - | killed(22) | killed |  |
| R03 | SaisKernel.cs:206 | counting stores presence instead of a count | killed(host crash) | - | killed(host crash) | killed |  |
| R04 | SaisKernel.cs:248 | SetBucketTails leaves the bucket of the largest symbol stale (claim: never used, the largest symbol is always L-type) | killed(9) | - | killed(9) | killed |  |
| R05 | SaisKernel.cs:262 | round 1 placement starts from the symbol of position n - 2 (reads index -1 for n = 1) | survived | killed(9) | survived / killed(7) | killed |  |
| R06 | SaisKernel.cs:285 | CompactSortedLms: the clamp of an empty slot is dropped (reads packed[-1] for p = -1) | survived | survived | - (survivor of every run) | **survived; equivalent** | E06h |
| R07 | SaisKernel.cs:336 | reduced problem: the LMS positions in text order are written one slot lower (spare slot: claimed equivalent because m <= (n - 1) / 2) | survived | survived | - (survivor of every run) | **survived; equivalent** | E07 |
| R08 | SaisKernel.cs:115 | representation choice: packed bytes also for 128 symbols (127 << 1 / 1 = 255 still fits: claimed equivalent) | survived | survived | - (survivor of every run) | **survived; equivalent** | E08 |
| R09 | SaisKernel.cs:348 | translation of the reduced ranks reads the tail one entry too high | killed(host crash) | - | killed(host crash) | killed |  |
| R10 | SaisKernel.cs:363 | round 2 placement: the bucket cursor is not stored back | killed(host crash) | - | killed(host crash) | killed |  |
| R11 | SaisKernel.cs:152 | stage 1 sorts with the heads of the buckets in the S pass (tails not re-set) | killed(host crash) | - | killed(host crash) | killed |  |
| R12 | SaisKernel.cs:324 | recursion receives the wrong alphabet size (m instead of names) | killed(8) | - | killed(7) | killed |  |
| R13 | SaisKernel.cs:319 | the reduced SA block is one entry short (m instead of m + 1) | survived | killed(31) | survived / killed(25) | killed |  |
| R14 | SaisKernel.cs:219 | classification: the last position is typed S-type for the equal-run rule only (nextIsS starts as 1) | killed(28) | - | killed(22) | killed |  |
| R20 | SaisInduction.cs:32 | L pass: the clamped prefetch index is dropped (prefetch of index -2 for empty entries: no architectural effect, claimed equivalent) | survived | survived | - (survivor of every run) | **survived; killed by added test** | prefetch pin |
| R21 | SaisInduction.cs:76 | naming: previous starts at 0 instead of -1 | killed(6) | - | killed(5) | killed |  |
| R22 | SaisInduction.cs:101 | naming: the end test compares with b instead of a | killed(27) | - | killed(22) | killed |  |
| R23 | SaisInduction.cs:87 | naming returns the last name instead of the count | killed(31) | - | killed(25) | killed |  |
| R24 | SaisInduction.cs:52 | S pass: look-ahead entry from the wrong end (Math.Min instead of Math.Max, reads past n) | survived | survived | - (survivor of every run) | **survived; killed by added test** | prefetch pin |
| R25 | SaisInduction.cs:81 | naming: the name slot of the look-ahead is computed from the position instead of the entry (no architectural effect, claimed equivalent) | survived | survived | - (survivor of every run) | **survived; performance only** | Buf-level prefetch, not interceptable |
| R26 | SaisInduction.cs:29 | L pass: loop starts at 1 (the sentinel suffix n is never scanned) | killed(30) | - | killed(24) | killed |  |
| R30 | TypedText.cs:95 | packed sentinel byte is hostile (0xFF): claim the sentinel is never read | survived | survived | - (survivor of every run) | **survived; equivalent** | E30 |
| R31 | TypedText.cs:62 | split sink writes the sentinel as L-type (claim: never read) | survived | survived | - (survivor of every run) | **survived; equivalent** | E30 |
| R32 | TypedText.cs:97 | packed: the type bit is OR-ed into the symbol position (isS << 1) | killed(26) | - | killed(20) | killed |  |
| R33 | TypedText.cs:43 | split Get reads the type of the neighbour | killed(host crash) | - | killed(host crash) | killed |  |
| R40 | SaisPrefetch.cs:13 | induction prefetch distance 0 (look-ahead is the current entry: no architectural effect, performance only) | survived | survived | - (survivor of every run) | **survived; killed by added test** | prefetch pin |
| R41 | SaisPrefetch.cs:13 | induction prefetch distance int.MaxValue (i + distance overflows to a negative index) | killed(host crash) | - | killed(host crash) | killed |  |
| R42 | SaisPrefetch.cs:16 | naming prefetch distance int.MaxValue (overflow) | killed(host crash) | - | killed(host crash) | killed |  |
| R43 | SaisPrefetch.cs:13 | induction prefetch distance 2^30 (i + 2^30 does not overflow below n = 2^30, clamped to n: claimed equivalent) | survived | survived | - (survivor of every run) | **survived; killed by added test** | prefetch pin |
| R44 | SaisPrefetch.cs:16 | naming prefetch distance negative (look-ahead reads before the first sorted LMS entry) | survived | killed(30) | survived / killed(24) | killed |  |


R80-R86 live in the `#if ESA_BOUNDS` section of `Buf.cs`: they cannot change the Release build ("survived" there by construction); their verdict is the checked column.


**Survivors of the pre-existing suite that are now killed (8).** R82-R85 (checked `Buf` bounds, MAJ-2): `ReviewS1aOptTests.CheckedBuf_*`. R20 (clamp of the prefetch index dropped), R24 (S-pass look-ahead from the wrong end), R40 (induction distance 0), R43 (distance 2^30): `ReviewS1aOptPrefetchTests`; each verified both ways (the new tests pass on the real code and fail on the mutant; of the earlier lists P01-P05 and P07 die as well).

**Survivors proven equivalent (9).** R04, R06, R07, R08, R30 and R31 by a hostile-variant RUN of the exhaustive engine-level tests; R50, R65, R66, R72 by inspection plus survival of the whole suite in both builds (R65/R66 in the checked build, where a store at index `cap` would throw).
* R04 / E04h (tail of the largest symbol left stale): the largest symbol is L-type everywhere (a run of it ends at a smaller symbol or at the text end), so `bkt[sigma-1]` after `SetBucketTails` is never used. The mutant is *killed* (9 tests) only by the intermediate-state oracle, which compares the bucket array; the output-equivalence is proven by writing `int.MinValue` there and running the engine-level exhaustive tests (8.4 M binary + 2.4 M ternary + 5.6 M 4-ary strings, each SA against the naive SA, guarded arena): 3/3 pass (`equivalence_probes.log`).
* R06 / E06h (clamp `p & -positive` in `CompactSortedLms` removed): E06h additionally THROWS if an empty slot (`p < 0`) is ever read there: 3/3 exhaustive pass, so after round 1 no empty slot exists and the clamp is dead code. The same holds in the checked build (R06 survives it too: no index -1 is read).
* R07 / E07 (`tailBase = n - m` instead of `n + 1 - m`): the head `[0..m]` and the tail stay disjoint because `m <= (n-1)/2`; exhaustive 3/3 pass.
* R08 / E08 (`sigma <= 128` selects the packed representation): `(127 << 1) | 1 = 255` still fits the byte; exhaustive 3/3 pass and the boundary sweep (sigma 127, 128, 129) passes in both representations.
* R30, R31 / E30 (sentinel byte of the packed / split text set to a hostile `0xFF` / `0`): `types[n]` and `packed[n]` are never read (all reads have index <= n-1); exhaustive 3/3 pass with both hostile.
* R50 (`l <= limit` calls `CommonPrefix` with `max = 0`): returns 0 by inspection (`head = 0`, no word, no tail).
* R65, R66 (no spare slot in `longRanks`, `escapeRanks`): the cap check precedes every store, the blocks of exact size survive the whole suite in the checked build (a store at index `cap` would throw there).
* R72 (placeholder `width + 1` instead of `width` for a long rank): every long rank is rewritten by `ResolveLong` (or the whole array is dropped on give-up), so the placeholder is never observable.

**Survivors that are performance-only (3).** R25 (name-slot prefetch address, a `Buf` prefetch), R64 (give-up at three quarters long instead of half), R87 (`prefetchnta` instead of `prefetcht0`; the Log measured the SA phase 75 % slower for `nta`). No output, state or access pattern changes. Not killable by a correctness test (MIN-2).

**Other observations.** Crash-kills (R01, R03, R09-R11, R33, R41, R42, R57, R58, R69, R70, R94) are real kills but a poor diagnosis; the checked build reports the exact access. R60 (`dst[0]` not written) is killed by 1 pre-existing test and by 2 of the new ones. R05, R13 and R44 survive Release and die only in the checked build (K14/K15 class). R91's first formulation (`payload - 1` for an empty block) wrapped to `ulong.Max` and crashed the host; the second formulation (`payload > 0 ? payload - 1 : 0`) is killed by a pre-existing test.

## 7. Extra correctness stress (item 5)

* **Differential fuzz, three new seeds** (`ESA_FUZZ_SEED` 101, 102, 103; five pipelines: SA pointer + LCP reference, SA reference + LCP three-pass, SA reference + LCP direct, pointer + three-pass, pointer + direct; each against the reference pipeline entry by entry, guarded arenas; scale 5, time budget 80 s per pipeline and tier, about 17 minutes per seed): per seed about 20.1-20.4 M cases with n <= 64 (acceptance: 1 M x 5), 2.0-2.1 M with n <= 2000 (acceptance 0.1 M x 5) and 25,000 with n <= 10^5 (acceptance 1,000 x 5), (0.62-0.65 G code units in the n <= 64 and n <= 10^5 tiers, 1.1-1.2 G in the n <= 2000 tier); **0 mismatches in 45 pipeline-tier runs** (`fuzz_summary.tsv`).
* **Vectors against libsais at 10 Mbp:** `Category=Vectors` (474 tests: 79 corpora x 6 pipelines, including `dna_10m`, `rep_10m`, every adversarial family at 10^5 and 10^6) green in all eight configurations of section 4. Independent of the manifest, my 23 adversarial 10 Mbp texts and 23 1 Mbp texts were compared with libsais SHA-256 of SA and LCP through the shipped chain, the three-pass chain and the reference LCP: **all 378 builds equal** (28 texts at 10 Mbp and 26 at 1 Mbp, 7 ESA builds each: 3 shipped, 3 three-pass, 1 reference LCP; `hash_ok = 1` in every row of the cold tables).
* **Sizes around block, unroll and vector boundaries**, every alphabet class: `BoundarySweep_SmallSizes_...` (default) and `BoundarySweep_LargeSizes_...` (`Category=ReviewSweep`): n in {0,1,2,3,4,5,7,8,9,15,16,17,31,32,33,47,48,49,63,64,65,127,128,129,255,256,257,511,512,513,1023,1024,1025} and {4095,4096,4097,65535,65536,65537,131071,131072,131073}, sigma in {1,2,3,4,5,16,20,126,127,128,129,200,255,256,257,1000,4096,65535,65536} (exactly sigma symbols realised, so the packed / split / 16-bit boundaries and the full 65,536-symbol alphabet are hit), five text shapes each (shuffled alphabet then random, then a long run with a different last letter, then a period-3 tail, the alphabet repeated, descending), every optimized SA and LCP engine against the reference, the naive oracle (n <= 2000) or the verifier and Kasai, guarded arenas, canaries verified, one or four poison patterns: green in the raw and the checked build, and the raw build under `EnableHWIntrinsic=0`, `EnableAVX2=0`, `EnableSSE42=0`, `TieredPGO=0`, `TieredCompilation=0`, and the checked build under `EnableHWIntrinsic=0` and `EnableSSE42=0` (section 4).
* **Other new checks:** all 256 poison bytes (`EveryPoisonByte_...`, 47 texts); heap dirtied with 0xFF arrays of the sizes the kernels allocate uninitialised (`ProductionEntryPoints_...`, n up to 200,000); four threads (`ConcurrentBuilds_...`); compacting GCs forced by a second thread while the kernels run (`BuildsUnderGcPressure_...`).

## 8. What I verified as correct (with evidence)

1. All Log numbers of 1o.3b reproduce (section 1): paired ratio 0.789 against 0.771, LCP phase 0.075 against 0.068 s, 12.3 B/char against 12.4, ESA / libsais 1.27-1.36x against 1.28x.
2. Every unsafe access has a stated invariant (section 3) and the invariants held under 64 mutants, 256 poison bytes, 45 fuzz runs of about 69 M cases in total, the exhaustive enumeration of 16 M strings, and the boundary sweep.
3. The kernels keep no static state and are safe for concurrent builds of different texts; all pointers into managed memory are pinned for their whole use.
4. The prefetch look-ahead reads are clamped everywhere (P01-P03, Q01-Q05, D21, D22, R41, R42, R44, R57, R58, R69, R70 all die), and now the targets and constants are pinned.
5. Output equals libsais (SA and LCP hashes) on 23 + 23 adversarial inputs and 79 corpora; equals the reference on every fuzz case.
6. The JIT claim of the Log reproduces: `jit_check.sh` on `ScanDirect`, `PlcpLoop`, `FillPhi`, `WriteU16` (Tier1-OSR listings): 0 range-check helpers, 0 accessor calls.
7. The earlier audit's six "equivalent" 1o.1b survivors are consistent with my independent probes (S06 = R30/R31 here, S17/S19 stale `SA[n]` = E04h-style reasoning on the largest symbol); I did not find a counter-example.
8. Linear behaviour on invalid input (MIN-8) and on every adversarial family (section 5).

## 9. Added files and tests

| file | content |
|---|---|
| `tests/.../Review/ReviewS1aOptTests.cs` | `CheckedBuf_*` (3), `GuardBands_DetectWrittenZeroAndFfByte_...`, `PoisonCoversTheWholeBlock_...`, `EveryPoisonByte_...`, `ProductionEntryPoints_DoNotDependOnTheInitialContentOfTheirOutputArrays`, `ConcurrentBuilds_...`, `BuildsUnderGcPressure_...`, `BoundarySweep_SmallSizes_...` (default); `BoundarySweep_LargeSizes_...` (`Category=ReviewSweep`, explicit); `ReportInstructionSets` (explicit) |
| `tests/.../Review/ReviewS1aOptPrefetchTests.cs` | the three prefetch-target pins (L pass, S pass, naming) |
| `tests/.../Review/ReviewS1aOptWorstCaseTests.cs` | `Direct_OnAnEmbeddedLongRepeat_IsLinear_...`, `Direct_OnDuplicatedBlocks_...`, `InvalidSuffixArrays_InRange_...` (default); `WorstCaseTimings` (explicit, `Category=ReviewTiming`); `Direct_Succeeds_OnRandomDnaWithOneEmbeddedRepeatOfOnePercent` (explicit, **fails today**: documents MAJ-1) |
| `analysis/mutation/review_opt_*.json` | the 64 mutants |
| `analysis/review_opt/` | raw data: cold/warm timings, A/B sessions, mutant table, equivalence probes, fallback runs, fuzz summary, libsais driver, campaign script |

Default Esa tests: 273 (was 257), raw and checked, 0 warnings; the older suites were not touched.
