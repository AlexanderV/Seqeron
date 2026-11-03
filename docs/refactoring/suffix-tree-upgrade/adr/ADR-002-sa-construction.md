# ADR-002: Suffix-array construction (SA-IS), LCP, memory plan and optimization roadmap

Status: **proposed** (work package S1-0, 2026-11-03). Rules: [PROCESS.md](../PROCESS.md) (managed only, no copied libsais code, optimization stopping rule at `PROCESS.md:25-29`). Contract and ESA model: [ADR-001](ADR-001-index-contract.md). Container: [ADR-003](ADR-003-esa-container-format.md). Numbers: [BENCHMARK.md](../BENCHMARK.md). **UNVERIFIED** marks anything not checked against code or a measurement; no benchmark was run for this ADR (the machine is busy), so every performance figure below is an estimate or a quoted measurement of another tool.

## 1. Gate and reference numbers

S1a gate (10 Mbp random DNA `dna_10m`, 1 thread, this machine class): SA + LCP build <= 3.3 s, index <= 16 B/char including text, peak RSS during build <= 24 B/char (about 240 MB). Ceiling goal 1.0 s, aspiration 0.5 s ([PROGRESS-S1.md](../PROGRESS-S1.md)). Reference points measured on this machine (`BENCHMARK.md:63-65`): libsais 0.328 s, libdivsufsort 0.662 s, Go `index/suffixarray` (SA-IS, bounds-checked) 0.73 s, Ukkonen persistent 14.9 s. A managed SA-IS in the same class as the Go implementation would land near 0.7-1.2 s (UNVERIFIED), so the 3.3 s gate leaves margin for a safe-code first version and the 1.0 s goal needs part of the roadmap in section 5.

## 2. Algorithm choice: SA-IS, int32 first, int64-ready

Chosen: **SA-IS** (induced sorting of LMS substrings with recursion on the reduced string; Nong, Zhang, Chan, "Two efficient algorithms for linear time suffix array construction", IEEE Trans. Computers 60(10), 2011), implemented from the paper. libsais, sais-lite and libdivsufsort are used only as binary oracles and as sources of ideas that are described publicly (prefetching, radix first pass, narrow types); **no source of them is read-and-copied into the repository** (PROCESS.md rule 1).

| Alternative | Why not |
|---|---|
| libdivsufsort-style (multikey quicksort + induction) | not linear, much more code, 2x libsais time (`BENCHMARK.md:65`) |
| DC3 / skew | 2-3x slower than SA-IS in practice (UNVERIFIED), more memory |
| prefix doubling (the internal one in `SequenceComplexity.cs:1909`, O(n log^2 n), cited in the K0 gaps, C7) | wrong complexity for 10-100 Mbp |
| SA from the Ukkonen tree | the thing S1 replaces (14.9 s) |

Concrete variant (every point is a design decision; the verification oracle for each is in [PROGRESS-S1.md](../PROGRESS-S1.md)):
1. **Output layout.** `SA` has `n + 1` entries and `SA[0] = n` (virtual smallest sentinel, never stored in the text). This is the natural output of SA-IS with an explicit sentinel suffix, equals the libsais result shifted by one (`SA[1..n]` = libsais `SA[0..n-1]`), and gives the terminated-tree model of ADR-001 section 4 without any copy.
2. **Alphabet preparation.** One pass counts distinct UTF-16 code units (65,536-bit map). Level-0 symbols are read through a struct reader `ISymbolReader` with three implementations: `byte` (alphabet <= 256 after dense remapping, order-preserving), `ushort` (larger alphabets), `int` (recursion levels). The JIT specialises each instantiation; no virtual dispatch in the kernels. The remap is order-preserving, so the SA is the SA of the original text; the stored text keeps original code units ([ADR-003](ADR-003-esa-container-format.md), TEXT codecs).
3. **Types.** 1 bit per position (`ulong[]`, n/8 bytes), classified right to left in one pass with word-wise operations; no per-position byte. (Alternative to evaluate in the roadmap, item R6: LMS flag in the sign bit of SA entries, saves the bit array and a cache stream.)
4. **Buckets.** `count`, `bucketStart`, `bucketEnd` arrays of size sigma. sigma <= 4 (DNA after remap): kept in locals/registers (roadmap R4). sigma <= 256: `stackalloc` (3 KB). Larger: pooled arrays. Recursion levels: arrays of size `names` (<= m).
5. **Recursion** exactly as in the paper: LMS substrings are sorted by one induce round, named (`names` stored at `SA[m + (pos >> 1)]`, unique because two LMS positions differ by >= 2), compacted into the reduced string `s1` placed in the **tail** of `SA` (`SA[n+1-m .. n+1)`), `SA1` is built in the head `SA[0..m)`; if `names == m` the reduced SA is placed directly (`SA1[s1[i]] = i`). Depth <= log2(n) frames (m <= n/2), so no stack risk; an explicit-stack driver is unnecessary (UNVERIFIED: confirm max depth on `a^n` and Fibonacci in tests).
6. **Index width.** *(Amended 2026-11-03, work package 1a.5, review MAJ-5; the earlier text made the optimized kernels generic over an index type with prototype P4.)* The optimized kernels are **`int`-only**: maximum speed, no generic-math risk on the hot path, and no P4. Only the safe reference is generic: `SaisReference<TIdx>` and `LcpReference<TIdx>` with `TIdx : unmanaged, IBinaryInteger<TIdx>, ISignedNumber<TIdx>` (all index conversions checked, `SaisReference`/`LcpReference` are the `int` facades). The `long` instantiation is an oracle for the index arithmetic: it must produce the identical SA and LCP (equality-tested on all small corpora in `IndexWidthTests`). The `long` kernels for S2 are generated from the `int` source by a script (checked-in, with a CI diff check; no source generators, no runtime codegen) and equality-tested against `SaisReference<long>` there. The reference is limited by `int`-indexed spans (`n + 1 <= Array.MaxLength`) whatever `TIdx` is.
7. **Limits.** `int` path: `n + 1 <= 2,147,483,591` (`Array.MaxLength`); the managed `string` limit is lower (about 1.07 Gchar, K0 gaps, "Hard caps"). Larger inputs throw `NotSupportedException("use the 64-bit path")`; the 64-bit path (segmented/native buffers, a BCL `NativeMemory`-style allocator is a decision for S2) is S2.

## 3. Workspace and memory plan (meeting peak RSS <= 24 B/char)

Per input character, random DNA, `n = 10^7`, all numbers B/char (estimates from the array sizes below; to be confirmed by the harness):

| Phase | Live data | B/char |
|---|---|---|
| input | caller's `string` (UTF-16) 2.0 (counted by RSS if the harness holds it; the persistent/stream build in S2 removes it) | 2.0 |
| A: SA-IS | byte text copy 1.0; `SA` `int[n+1]` 4.0; type bits 0.125; level-1 buckets <= `4*names` (random DNA: m about n/3, names about 0.9 m, so about 1.2; hard bound 2.0 when m = n/2) | 8.3 (bound 9.1) |
| B: LCP | + work `int[n+1]` for Phi, overwritten in place by PLCP (4.0); compact `LCP` bytes 1.0 (+ overflow list) | 12.0 |
| C: child table (S1b) | built **into the work buffer** after PLCP is dead: no new allocation | 12.0 |
| D: navigation pack (S1b, optional) | + ISA 4.0 + RMQ about 0.4; the work buffer may hold ISA | about 16.5 with the string held, 14.5 without |

Peak is phase B/C/D, about 12-17 B/char plus the .NET runtime baseline (30-60 MB, UNVERIFIED, i.e. 3-6 B/char at 10 Mbp), so the 24 B/char gate has about 1.5-4.5 B/char of headroom for the Full profile (24 - 16.5 - 3..6). Allocation rules: one `SA` and one work buffer allocated with `GC.AllocateUninitializedArray<int>(n + 1, pinned: true)` (pinned object heap, no zeroing, no LOH compaction, stable pointers for the unsafe kernels); no allocation inside loops or per recursion level (recursion reuses SA tail/head and one pooled bucket array per level); results are published as `ReadOnlyMemory`/pointer wrappers, never copied. When writing a container ([ADR-003](ADR-003-esa-container-format.md)) the SA is written section by section from these buffers (no second copy).

## 4. LCP algorithm

**Phi / PLCP** (Karkkainen, Manzini, Puglisi, "Permuted longest-common-prefix array", CPM 2009), 3 passes, linear:
1. `Phi[SA[i]] = SA[i-1]` for `i = 1..n` (for `i = 1` the predecessor is the empty suffix `n`): one sequential read of SA, one random write.
2. In text order `j = 0..n-1`: `p = Phi[j]`; if `p == n` then `l = 0`, `PLCP[j] = 0`; else extend `l` while `T[j+l] == T[p+l]` (bounded by `n`), `PLCP[j] = l`, `l = max(l - 1, 0)`. `PLCP` overwrites `Phi` in place (each `Phi[j]` is read exactly once before it is overwritten). Total compared characters <= 2n; the compare loop uses a vectorised common-prefix primitive (`MemoryExtensions.CommonPrefixLength` on byte spans, vectorised in recent .NET: UNVERIFIED for 10.0.403) for runs >= 16.
3. `LCP[r] = PLCP[SA[r]]` for `r = 1..n` (`LCP[0] = LCP[1] = 0`): sequential SA read, random PLCP read, written **directly into the compact LCP** (below), so no dense 4 B/char LCP ever exists.

Compact LCP (kinds recorded in the container, [ADR-003](ADR-003-esa-container-format.md)): `U8` with escape value 255 and a sorted overflow list `(rank, value)`; `U16` with escape 65535; `U32`. The builder picks the smallest representation: `U8+overflow` if `escapes * 8 <= n/2`, else `U16(+overflow)`, else `U32`. *(Clarified 2026-11-03, work package 1a.4.)* The escapes are counted on PLCP **before** pass 3 (the multiset of LCP values equals that of PLCP, so the counts are the same), because the representation must be known before the first entry is written. With `e8` = number of values `>= 255`, `e16` = number of values `>= 65535` and `n` the text length the rule is: U8 iff `16 * e8 <= n` (the literal rule: overflow list at most half of the base array); U16 iff `8 * e16 <= n` (the same condition one tier up, overflow list at most half of the 2n-byte base array: the ADR gave this tier only as "else U16(+overflow)"); else U32. Consequences measured at 10^6: random DNA is U8 with 0 escapes; `a^n`, `abab...`, Fibonacci and Thue-Morse words are U32 (most LCP values exceed 65,535); `a^70000` is U16 with 4,465 escapes. Random DNA: about 1.0 B/char (escapes are astronomically rare); `rep_10m` and genomes: a few percent escapes; `a^n`: U32. The sequential consumers (child table, RMQ build) read the compact LCP with a cursor into the overflow list (amortised O(1)); random access is binary search in the overflow list (used only on the rare >= 255 values).

Alternatives (not in S1a baseline): Gog and Ohlebusch's Phi variants (ALENEX 2011) and Fischer's inducing of LCP during SA-IS (WADS 2011) (details UNVERIFIED); consider only if LCP exceeds 15 % of the build time after R2/R3 (roadmap R7).

## 5. Ranked optimization roadmap

Order = expected gain / risk, to be run as iterations under the stopping rule of `PROCESS.md:25-29` (accept only if >= 3 % faster on the 10 Mbp row, median of 5, all tests and fuzz green; stop after 3 consecutive iterations < 3 % or at <= 0.36 s). Gains are **estimates** (UNVERIFIED), expressed on top of a correct safe-code baseline; the total headroom is bounded by the spread libsais 0.328 s vs Go 0.73 s (`BENCHMARK.md:63-64`), about 2.2x, plus managed overhead.

| Id | Idea | Expected gain | Risk | Verification |
|---|---|---|---|---|
| R0 | **Baseline hygiene** (part of 1a, not an iteration): pinned uninitialised arrays, one workspace, no per-level allocation, generic struct readers, `Span`-based safe kernels kept as `SaisReference` forever | n/a | low | SA equal to brute force and libsais hashes |
| R1 | **Unsafe pointer kernels / bounds-check removal** in classification, induce L/S, naming, compaction (`int*` over pinned arrays, `Unsafe.Add`, no `Span` indexers in the inner loops) | 1.2-1.5x | low | bit-identical SA vs `SaisReference` on all corpora and 10^6 fuzz cases; `DOTNET_JitDisasm` check that no `RngChkFail`/`CORINFO_HELP_RNGCHKFAIL` remains in the listed methods |
| R2 | **Software prefetch** in induce loops (prefetch the text symbol and bucket cell of `SA[i + d]`; `Sse.Prefetch0` guarded by `Sse.IsSupported`, portable fallback = no-op; ARM64 prefetch intrinsic availability is UNVERIFIED, then also no-op); distance `d` chosen by a small sweep | 1.2-1.8x on large inputs; smaller on this VM because it reports a 480 MiB L3 (`BENCHMARK.md:16`), so 10 Mbp structures largely sit in cache | medium (tuning is machine specific) | A/B on 1, 10, 50 Mbp; identical SA; accepted on the 10 Mbp rule, 50 Mbp recorded as information |
| R3 | **Prefetch in LCP passes** (Phi scatter, PLCP gather) and merge pass 3 with the compaction | 0.05-0.15 s of the total | low | A/B; identical LCP hash |
| R4 | **sigma <= 4 fast path** (DNA after remap): bucket pointers as four locals, 2-bit symbol codes, branch-free bucket select; instantiated by a `struct Sigma4 : IBuckets` so the generic code is reused | 1.1-1.3x on DNA | medium (duplication risk, mitigated by the generic struct) | DNA corpora + equality with the general path on every fuzz case |
| R5 | **Radix first pass for LMS prefixes**: bucket the LMS suffixes by their first 2-3 symbols (DNA: 16-64 buckets) before induction, so the first induce round starts from sorted prefixes | 1.1-1.2x | medium | equality; A/B |
| R6 | **LMS/type flag in the SA sign bit** (removes the type bit array and one random stream; valid because `n + 1 < 2^31`) | 1.05-1.15x | medium (bit tricks, easy to break near `n = 2^31`) | equality incl. a synthetic test at the limit with a fake reader (no large allocation) |
| R7 | **LCP alternatives** (Fischer induced LCP, Gog-Ohlebusch Phi variants) if LCP > 15 % of total | up to 0.1 s | medium | LCP hash |
| R8 | **JIT warm-up control** (section 6): `AggressiveOptimization` on hot kernels vs tiered PGO | 1 Mbp builds: tens of % ; 10 Mbp: < 3 % | low | A/B both ways, cold-process and steady-state rows |
| R9 | **Narrow types in recursion**: when `names < 65,536` use a `ushort` reader for the reduced string (half the memory traffic of the second level) | 1.02-1.05x | low | equality |
| R10 | Vectorised type classification (run detection with `Vector128/256` compare + prefix-xor on 64-bit blocks) | <= 3-5 % (classification is a small share) | low | equality; probably below the 3 % acceptance bar, kept last |

Deferred (named so they are not forgotten): **S2**: 2-bit packed DNA text and packed-word kernels; integer alphabets; the 64-bit SA path with segmented or native buffers (decision on the BCL `NativeMemory` API versus segmented managed arrays: it is not P/Invoke in our code, but is an owner call because of the rule "no native"); build from FASTA/stream without a managed string; external-memory build. **S5**: multi-threading (bucket-partitioned induction, parallel Phi/PLCP by text ranges, parallel child-table build; byte-identical output to the serial build). **S7**: BWT/FM structures derived from this SA.

## 6. JIT versus NativeAOT / ReadyToRun (design guidance only, no benchmarks run)

Position (owner decision on AOT as a build artifact is pending): the shipped library stays **managed IL**; no native binary is a dependency. NativeAOT is allowed only (a) as a measurement variant of the benchmark harness and (b) later for CLI/MCP packaging in S8. ReadyToRun is treated the same way (packaging only, IL stays in the assembly).
1. Steady-state throughput is expected to be about equal under the JIT (tiered PGO, hardware-detected `Vector128/256`) and under NativeAOT, so no large gain is expected from AOT for the kernels (UNVERIFIED expectation).
2. The real JIT cost is warm-up of the SA-IS hot loops (tiering and on-stack replacement of loop-heavy methods): this matters for 1 Mbp builds, which take tens of milliseconds (libsais 0.026 s, Go 0.055 s at 1 Mbp, `BENCHMARK.md:45-46`). Evaluate `[MethodImpl(MethodImplOptions.AggressiveOptimization)]` on the hot kernels: it skips tiering but also opts the method out of dynamic PGO, so it must be measured **both ways** on 1, 10 and 50 Mbp, in a fresh process (cold) and after warm-up (steady state).
3. The new `SuffixTree.Esa` project is marked `IsAotCompatible=true` (and `IsTrimmable=true`, as `SuffixTree.Core.csproj` and `SuffixTree.Persistent.csproj` already are) so the trimming/AOT analyzers run in every build; the ESA code uses no reflection and no runtime code generation (no `Reflection.Emit`, `Expression.Compile`, `Activator` by name, `dynamic`). Struct-constrained generics and static abstract interface members are AOT-compatible.
4. Measurement protocol ([PROGRESS-S1.md](../PROGRESS-S1.md)) gets an **optional variant row**: the same harness published with NativeAOT (`PublishAot`, instruction set `native`), recorded under a separate tool tag, never mixed into the JIT rows (needs the ILC package: availability offline is UNVERIFIED, so the row is optional). **An optimization is accepted on the JIT numbers only.**

## 7. Verification summary

Algorithmic correctness is established, in this order, by: brute-force SA/LCP for all strings up to a length bound and random strings up to n = 2000; equality with `SaisReference`; SA equality with libsais outputs on vectors generated outside the repository (hash manifest); the Ukkonen tree as an independent order oracle (`EnumerateSuffixes` is sorted); adversarial families; and the differential fuzz harness. Details and acceptance criteria: [PROGRESS-S1.md](../PROGRESS-S1.md), work packages 1a.1-1a.9 and 1o.*.

## 8. Risks

* R1 (high): the 1.0 s goal needs R1-R5 and may be limited by managed codegen; mitigated by the 3.3 s gate and the stopping rule; if the first safe version exceeds 3.3 s, stop and report (PLAN.md:81).
* R2 (closed 2026-11-03, 1a.5): generic-over-index kernels might have cost JIT quality (P4); resolved by making the optimized kernels `int`-only (section 2 item 6), so P4 is not run.
* R3 (medium): this VM hides cache effects (480 MiB L3); prefetch gains are measured here but could be larger elsewhere; the S2 50-100 Mbp runs are where they pay.
* R4 (low): `AggressiveOptimization` and tiered PGO interact; handled by the A/B rule in section 6.
