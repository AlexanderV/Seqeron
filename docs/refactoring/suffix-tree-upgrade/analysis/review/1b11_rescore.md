# 1b.11: re-score of the Seqeron column (K0) of COMPETITOR_MATRIX.md (2026-11-23)

Scorer: an independent agent that did not read PLAN.md from section 4 on, `analysis/plan_model.py` or `analysis/plan_model_output.txt`. The head at the start was `50936932`
(the production code of S1 is that of `12495c3e`; `git diff 12495c3e..HEAD -- src apps` is empty). Rubric: `analysis/SPEC.md` and `analysis/FUNCTIONS.md`, K0 only. Competitor cells, weights, function list
and groups were not touched. Evidence files changed: `analysis/audit/K0_audit.tsv` (C8, C14), `analysis/audit/bench_overrides.tsv` (C10-C13, C15), `analysis/audit/K0_func.tsv` (21 functions); each new row is
appended below the old one with the date 2026-11-23 and wins in `build_docs.py` and `merge_func.py` (later row of the same cell), so the old rows stay in the files as the history. The matrix was regenerated with
`python3 -I build_docs.py`.

## 1. Result

| | before (K0 at commit 24231842, scored 2026-11-02/03) | after (this re-score) |
|---|---|---|
| weighted property mean (24 properties, 0-4) | 1.80 (sum 97.2 of 54 weights) | **2.28** (sum 123.2) |
| function coverage, whole library (81 functions, 162 points) | 55.6 % (90 points) | **72.8 %** (118 points) |
| function coverage, tree layer only (rows whose evidence starts "Genomics:" set to 0) | 39.5 % (64 points) | **56.8 %** (92 points) |
| rank by property mean among the 13 tools | 9 | **4** (ropebwt3 2.48, GenomeTools 2.30, BWA/Bt2 2.30, Seqeron 2.28, r-index 2.13) |
| functions below the best competitor | 44 of 81 | 29 of 81 |
| ahead of every competitor | F25, F47, F68, F71 | the same four |

Cells: **28 up, 0 down, 77 unchanged** (105 cells: 24 properties + 81 functions). Properties: 7 up, 17 unchanged. Functions: 21 up, 60 unchanged.
Nothing that the old code had is gone: the Ukkonen tree and the v6 persistent tree are still there and unchanged in behaviour (S0 only added guards, the strict validator and the MCP size limit), so no cell moves down.
Where the new engine is narrower than the old one it is named below (the container holds byte texts only; the first `Open` is slower than the v6 cold load) and the old engine carries the cell.

## 2. Changed cells

Every figure is a row of `analysis/bench/results_1b9_*`, a line of `ADR-004` section 5a "1b.9" (and its "Orchestrator verification" paragraph) or an earlier gate row, or a line of code / a test, as named.
"best" and "libsais" are the existing rows of `BENCHMARK.md` / `bench.tsv` (same machine class); the same-session libsais of 1b.9 is `results_1b9_libsais.tsv`
(`build time` 0.3385 / 0.3454 / 0.3126 s at the start, middle and end, pooled 0.3442 s, against 0.328 s of 2026-11-03: no drift).

### 2.1 Properties (7 up)

| id | old -> new | evidence sentence (as written in the evidence file) | source | missing for the next level |
|---|---|---|---|---|
| C8 | 2 -> 3 | Public and tested against brute force: EsaSuffixTree implements ISuffixTreeNavigation (explicit node, parent, child by symbol, children, suffix link, edge symbol, string depth) and ISuffixArrayAccess (SA, ISA, LCP, Lcp(i,j), Lca) besides Traverse; no Weiner link, no level ancestor. | `src/SuffixTree/Algorithms/SuffixTree.Core/Indexing/ISuffixTreeNavigation.cs:6`, `ISuffixArrayAccess.cs:7`, `SuffixTree.Esa/EsaSuffixTree.Navigation.cs:97-167`; tests `tests/SuffixTree/SuffixTree.Esa.Tests/Contract/EsaContractChecker.cs`, `IndexConformanceTests.cs:488` (every member against the definitions, in-memory and file-backed) | 4 is the sdsl reference CST (K1 4): Weiner links, level ancestor, a node handle type with DFS / bottom-up iterators; the navigation pack (ISA, LCP pyramid, 4.13 B/char) is built on first use |
| C10 | 1 -> 3 | Measured 10 Mbp: SA+LCP+child table 1.138 s [1.118-1.174] = 3.5x libsais 0.328 s (36 cold builds 0.927 s = 2.8x; public CreateEsa 0.841 s); 50 Mbp 5.23 s = 2.5x libsais; single thread, was 44.6x | `bench/results_1b9_rows.tsv:9` (`build_total_s`), `:44` (`entry_build_s` 0.841 s), `:388` (`lean_build` 0.895 s), `:173` (`dna_50m` 5.227 s), ADR-004 1b.9 table A1, `bench.tsv` libsais 50 Mbp 2.087 s | 4 needs <= 1.5x of libsais (here 2.4-3.5x); a build below 0.5 s at 10 Mbp (aspiration of PROGRESS-S1 1.5x); the 4-thread container build is not public API |
| C11 | 1 -> 2 | Measured 10 Mbp: 10.00 B/char in memory (text 1 + SA 4 + LCP 1 + child 4), container Core 9.25 / Full 13.25 B/char, build peak VmHWM 15.4 B/char; was 169 B/char; SA 4, sdsl cst 1.65, csa 0.27 | `bench/results_1b9_rows.tsv:22` (`index_bytes_core_per_char`), `:12` (`peak_rss_bytes_per_char`), `results_1b9_sections.txt:11-13` | 3 = the 4-5 B/char of a plain suffix array (K5, K6, K8 score 3 at 4-5 B/char; K2 at about 10 B/char scores 2); 4 = within 1.5x of 0.27 B/char |
| C12 | 2 -> 4 | Measured 10 Mbp: count 647k q/s (1.15x of best 747k), locate8 990k pat/s (1.15x of rust-bio 1.14M), geometric mean 1.15x; 4 threads 3.10M q/s; file-backed 0.99-1.03 of in-memory on DNA | `bench/results_1b9_rows.tsv:202` (`count_present_qps` 647,176), `:205` (`locate_findall_pps` 990,228), `results_1b9_abapp.txt` (`summarize_1b9.py --gate`) | nothing for the scale (4 is the top); open: the 20-letter text on the file-backed tree is 0.934 of the in-memory tree at four threads (ADR-004 1b.9 open point 1) |
| C13 | 1 -> 2 | Measured: container built at 10^8 symbols (17.9 s, RssAnon 253 MB, 13.25 B/char; 100 Mbp was killed at 600 s before); facade/reader limit 2,147,483,646 symbols, 3x10^9 not run (owner hand-off), 64-bit positions not built | `bench/results_1b9_chain.txt:9`, `src/SuffixTree/Algorithms/SuffixTree.Esa/EsaSuffixTree.cs:74` (`MaxTextLength`), `EsaLimits.cs:12` | 3 = documented working at >= 3 Gbp (SPEC): the owner's 3x10^9 chain (not run), a reader above 2^31 symbols and the 64-bit position type (S2) |
| C14 | 2 -> 3 | Measured: container build with bounded RAM (RssAnon 251-286 MB at 10^8 symbols, 13.25-14.25 B/char on disk, 17.9-26.5 s, 4 threads) and read-only mmap Open (15 ms); builder is internal (Bench CLI only), no cold-cache query timing, 3x10^9 not run. | `bench/results_1b9_chain.txt:9-12`, ADR-004 1b.9 table C, `External/Sesa/SesaBuilder.cs:307` (internal), ADR-004 "1b.9d" (the chain under a page-cache limit: 1.73x slower at an effective cache of 0.71, 5.05x at 0.50) | a public builder API (S8), a cold-cache query measurement, the 3x10^9 run |
| C15 | 0 -> 2 | Measured: container builder runs 4 worker threads by default; 100 Mbp 39.3 / 28.6 / 24.7 s at 1 / 2 / 4 threads (1.6x), 1 Gbp 247 / 151 / 123 s (2.0x), prototype figures; in-memory SA-IS build single-threaded | `analysis/em/e6_threads.tsv` (t1 39.316 s, t2 28.629 s), `analysis/em/e6_main.tsv` (t4, median of 3 = 24.696 s), ADR-004 "E8" (1 Gbp), `External/Sesa/SesaBuilder.cs:47` (`Threads = 4`), test `TheBytes_DoNotDependOnThreadsSizesSeedsOrArenaCaps` | 3 = documented near-linear speed-up on >= 8 threads; the 1 / 2 / 4 thread times are those of the E6 / E8 prototype of the same pipeline (the 1b.9 gate ran the production chain at 4 threads only: 17.88 s at 10^8); the in-memory build is single-threaded |

### 2.2 Functions (21 up)

All evidence sentences are in `analysis/audit/K0_func.tsv` (rows dated 2026-11-23, prefix "Esa:", so the tree-layer-only figure counts them). The public surface is `EsaSuffixTree`
(`SuffixTree.Esa/EsaSuffixTree.cs:28`: `ISuffixTree`, `ISuffixTreeNavigation`, `ISuffixArrayAccess`, `IDisposable`), opened by `EsaSuffixTreeFactory.CreateEsa|Open` and `SuffixIndexFactory.Open`.
Calibration against the other columns: a function that needs only a short loop of caller code over public SA / LCP is 1 (K4, K5 score F43 1 that way); a public accessor is 2 (K4, K9 score F07 / F08 / F59 2).

| id | function | old -> new | source | missing for the next level |
|---|---|---|---|---|
| F07 | Expose SA[i] / inverse SA | 1 -> 2 | `ISuffixArrayAccess.cs:10,13` (`SuffixAt`, `RankOf`) | - |
| F08 | Expose LCP array / PLCP | 1 -> 2 | `ISuffixArrayAccess.cs:16,19` (`LcpAt`, `Lcp`) | - |
| F10 | Expose child table / LCP-interval tree | 1 -> 2 | `ISuffixTreeNavigation.cs:9-21` | the raw child array stays internal (three layouts) |
| F15 | Locate with limit / lazy iterator | 0 -> 2 | `ISuffixIndex.cs:58` (`Locate(locus, skip, Span<long>)`; 2,173,239 patterns/s, `results_1b9_rows.tsv:206`) | `FindAllOccurrences` still materialises the list |
| F17 | Reusable handle for a pattern | 0 -> 2 | `IndexLocus.cs:32` (SA interval `Lo..Hi`, depth `Length`) | - |
| F18 | Extend a match to the left | 0 -> 1 | `ISuffixIndex.cs:50` (`TryExtendLeft`, `Support(ExtendLeft)` = Degraded, O(width of the interval)) | an O(1) backward step (an FM-index) |
| F19 | Extend a match to the right | 0 -> 2 | `ISuffixIndex.cs:41,47` (`TryExtendRight`) | - |
| F20 | Bidirectional extension | 0 -> 1 | `ISuffixIndex.cs:41,50` (both on one `IndexLocus`; the left side Degraded) | a bidirectional index |
| F36 | Matching statistics array | 0 -> 1 | `ISuffixTreeNavigation.cs:24` (`TrySuffixLink`) with `TryExtendRight` and `Count`; user code, as K1 / K4 / K7 are scored | a public `ms[]` routine |
| F43 | Shortest unique substring | 0 -> 1 | `ISuffixArrayAccess.cs:10,16`; a short loop over SA and LCP | a routine |
| F54 | Child by character | 0 -> 2 | `ISuffixTreeNavigation.cs:15` | - |
| F55 | Parent | 0 -> 2 | `ISuffixTreeNavigation.cs:12` | - |
| F56 | Suffix link | 0 -> 2 | `ISuffixTreeNavigation.cs:24` (kernel (c2): 410-415 ns per link at 10 Mbp, PROGRESS-S1 1b.4a) | - |
| F57 | Weiner link | 0 -> 1 | `ISuffixIndex.cs:50` (`TryExtendLeft` on a node locus, when the result is explicit) | an accessor |
| F58 | String depth of a node | 1 -> 2 | `IndexLocus.cs:32` (`Length`), `ISuffixTreeNavigation.cs:27` | - |
| F59 | LCA / LCP(i,j) | 0 -> 2 | `ISuffixArrayAccess.cs:19,22` (`Lcp`, `Lca`) | - |
| F60 | Level ancestor | 0 -> 1 | `ISuffixTreeNavigation.cs:12` (a loop of `TryParent`, O(depth)) | a routine |
| F61 | Leaf rank / leaf label | 1 -> 2 | `ISuffixArrayAccess.cs:10,13` | - |
| F63 | Enumerate children | 1 -> 2 | `ISuffixTreeNavigation.cs:21` (`GetChildren`) | - |
| F65 | Edge / path label access | 1 -> 2 | `ISuffixTreeNavigation.cs:27` (`EdgeSymbol`), `IIndexText.CopyTo` | - |
| F79 | 2-bit packed DNA | 0 -> 1 | `External/PackedTextFormat.cs:7` (`b = ceil(log2 sigma)` bits, 0.25 B/char at sigma 4, `results_1b9_sections.txt`) | the packing is written only by the internal builder and the container holds byte texts; the in-memory text is a byte per symbol |

Tests behind F07-F20, F36-F65 (the public contract members): the conformance suite `Contract/IndexConformanceTests.cs` over eight engines (Ukkonen adapter, three ESA layouts, two facades, file-backed `Checked` and `Trusted`),
`Contract/EsaContractChecker.cs` against `IndexContractOracle` (brute force from the definitions), the 1b.8 differential harness (every interface member, 111,227 cases in 10 minutes, 0 failures, ADR-004 1b.9 "Tests").

## 3. Considered and not moved

Properties (17 unchanged; those I weighed are listed).

| id | score | candidate | why it stays |
|---|---|---|---|
| C1 | 3 | 4 | not a measurable capability; `FindAllOccurrences` is still unsorted `int`; the new lazy `Locate` is `long` but not on `ISuffixTree` |
| C2, C3, C7 | 1 | 2 | nothing was added: the public primitives (`TryChild`, `TryExtendRight`) allow caller backtracking, which is the level-1 definition already used |
| C4 | 3 | 4 | MEM / MUM now run on the ESA (100 kbp query on 10 Mbp: 0.109 s first call, 0.075 s steady, ADR-004 1b.9 A8; MUMmer parity data pass), but there is no head-to-head with `mummer`, forward strand only |
| C5 | 2 | 3 | no public `ms[]` (the brief and the rubric); the public primitives make it caller code, which is what K1 / K2 score 1 for; `FindExactMatchAnchors` is unchanged |
| C6, C9 | 2, 3 | - | no new multi-text index, no new analytics API (the ESA computes the same operations natively; same list) |
| C16 | 3 | 4 | 3 is the highest level any column has; measured again: 3,098,222 count/s at 4 threads = 3.96x of 783,122 at one thread (same instrument), 10,240 thread runs of the fuzz with equal checksums |
| C17, C18 | 0 | 1 | the container is immutable and uncompressed (13.25 B/char); the 2-bit text is a text codec, not a compressed index |
| C19 | 2 | 3 | the 2-bit packed text exists only in the SESA container, whose builder reads **byte texts** (code units < 256); Unicode (UTF-16) only in the in-memory ESA and the v6 tree; no integer alphabets; K2 (<= 256 ranks, compact dna4) scores 2 |
| C20 | 3 | 4 | the SESA container adds CRC-32C per section, validation levels and an atomic, resumable build (all four parts of the level-3 definition were already met by v6); the first `Open` of a fresh process is 13-15 ms against its 10 ms gate (8.5 ms with ReadyToRun), the v6 cold load is 4.8 ms; 4 has no definition beyond best-in-class |
| C21 | 3 | 4 | the new contract (`long` positions, `Span<long>`, `NotSupportedException` naming the operation) is good, but the container builder is internal and the positions of `ISuffixTree` are `int` |
| C22 | 1 | 2 | no PackageId, no tag; the Docker hand-off image and the Bench commands are measurement tools; the MCP tools still use the Ukkonen tree |
| C23 | 3 | 4 | the volume is far above level 3 (`SuffixTree.Esa.Tests` default 4,945 + 37 skipped, `Vectors` 632 against libsais hashes, 111,227 differential cases, over a thousand mutants with a verdict, 1,161 in 1b.9c alone) and `dotnet test Seqeron.sln` in CI runs the default sets, but 4 is "best in class" without a measurement; the Exhaustive and Fuzz categories and the checked build are not in the CI (debt register of S1b, the owner's) |
| C24 | 2 | 3 | still no tag or release (`git tag` is empty), one maintainer |

Functions (60 unchanged; those I weighed).

| id | score | candidate | why it stays |
|---|---|---|---|
| F03 | 1 | 2 | the container builder reads its input through a read-only mapping with bounded RAM, but `SesaBuilder` is `internal` (InternalsVisibleTo the Bench, EmProto and test assemblies only); the only route is the Bench / `s1b.sh` command line, which the scale scores 1 ("CLI-only when an API is expected") |
| F09 | 1 | 2 | BWT, LF and psi are one-line caller expressions over `SuffixAt`, `RankOf`, `Text` (K4 scores 1); the container's PSI section has no public accessor |
| F16 | 2 | - | the 2 is the Genomics-layer evidence; `Locate` and `FindAllOccurrences` of the ESA are in suffix-array order, not ascending |
| F31, F32, F35, F39 | 1 | 2 | the public navigation permits caller backtracking and caller SMEM filtering; no routine was added |
| F44, F53 | 0 | 1 | no routine and, unlike F43, no short loop (K5 with SA + LCP scores F53 0) |
| F50, F52 | 1 | 2 | unchanged hash-based / LZ76 evidence |
| F66-F69 | 2 | - | already 2; the container adds CRC-32C and validation, v6 keeps the SHA-256 export |
| F70 | 1 | 2 | the Bench `sesa` commands and the hand-off scripts (`s1b.sh build`, `verify`, `chain`) are harnesses with fixed protocols, not a documented query CLI |
| F71 | 2 | - | the MCP server is unchanged (Ukkonen tree, rebuild per call, 5,000,000 characters per call by default since S0: a safety limit, configurable, not a missing capability) |
| F72, F73, F74, F75 | 2, 1, 2, 1 | - | no new statistics accessor beyond the old ones, no graph export, thread safety already 2 (the new `ReaderGate` protocol and tests only strengthen the evidence), no batch API |
| F76, F77, F78, F80, F81 | 2, 2, 1, 2, 1 | - | UTF-16 and byte texts as before; the container is byte-only; no integer alphabet |

## 4. Scale arithmetic (properties C10, C11, C12)

Ratios as in `bench_overrides.tsv`: Seqeron over the best measured tool, same machine class. The levels the other columns got: 1.6x (K8) and 2.2x (K6) -> 3, 5.0x (K1) and 5.7x (K7) -> 2, 15x (K10) -> 2, 44.6x (old K0) -> 1,
<= 1.5x -> 4 (K5, K7).

| cell | Seqeron (new) | best | ratio | other columns at that ratio |
|---|---|---|---|---|
| C10 | 1.138 s median [1.118-1.174] of the gate row (3.3x the pooled same-session libsais 0.3442 s, 3.5x the 0.328 s of 2026-11-03); 0.927 s median of the 36 cold builds of the day (2.8x); public `CreateEsa` 0.841 s (2.4x the same-session libsais, 2.6x the 0.328 s); `dna_50m` 5.227 s against libsais 2.087 s (2.5x) | libsais 0.328 s | 2.4-3.5x | between K6 2.2x (3) and K1 5.0x (2); Seqeron is faster than the two members that score 2 (sdsl csa 1.64 s, rust-bio 1.87 s) and builds the child table as well: 3 |
| C11 | 10.00 B/char in RAM, 9.25 (Core) / 13.25 (Full) in the file, 15.4 B/char peak of the process | sdsl csa_wt_rrr 0.27, plain SA 4 (+1 text) | 37x the smallest, 2x a plain SA | plain SA 5 B/char -> 3; K2 about 10 B/char -> 2 (optimistic): 2 |
| C12 | count of present 20-mers 647,176 q/s; locate of 8-mers 990,228 pat/s | libsais 747,415; rust-bio 1,140,278 | 1.155x and 1.152x, geometric mean 1.15x (old: 1.3x and 23.0x, 5.4x) | K7 at 1.5x -> 4: 4 |

Same-session control (`results_1b9_libsais.tsv`, `count_present20_1t` 806,827 / 749,086 / 911,684 and `locate8_1t` 690,854 / 665,996 / 880,281 q/s at the three points): Seqeron's 647k and 990k are 0.80 and 1.18 of the libsais binary search of the same day (ADR-004 1b.9 "Limits"), geometric mean 0.97.

## 5. Differences of the regenerated matrix against the old one

`python3 -I build_docs.py` from `analysis/`; the old file is the committed `COMPETITOR_MATRIX.md` of `50936932` (regenerating it before any change gave a byte-identical file). Comparison by
`analysis/review/1b11_diff_columns.py` (cell by cell over Tables B and A, the group and total rows, the ranking and the gap list):
`git show 50936932:docs/refactoring/suffix-tree-upgrade/COMPETITOR_MATRIX.md > /tmp/stbench/b11/old.md; python3 -I analysis/review/1b11_diff_columns.py /tmp/stbench/b11/old.md COMPETITOR_MATRIX.md`; its output:

```
section 0 (heading, scales), section 1 (competitors), section 6 (evidence): identical
Tables B and A: non-table lines identical; cells changed in the Seqeron column: 7 properties + 21 functions (+ 1 weighted mean, 5 group coverage rows, 1 total row)
cells changed in any other column: 0
ranking: every other tool keeps its two figures; the rank number only changes for r-index 4 -> 5, MUMmer4 5 -> 6, rust-bio 6 -> 7, sdsl 7 -> 8, Python 8 -> 9
Seqeron ranking line: 1.80 / 55.6 % / rank 9 -> 2.28 / 72.8 % / rank 4
coverage lines: "Suffix-tree layer only" 39.5 % -> 56.8 %; "Seqeron scores below the best competitor on 44 of 81" -> 29 of 81
gap list: 44 -> 29 entries; removed F07 F08 F10 F15 F17 F19 F54 F55 F56 F58 F59 F60 F61 F63 F65; changed (Seqeron column only) F18 F20 F36 F43 F57 F79
"Seqeron scores above every competitor on": unchanged (F25, F47, F68, F71)
```

Things in the generated page that the re-score could not change because `build_docs.py` is outside the files this package may edit (owner decision, see section 7): the K0 line of the competitor table still reads
"Seqeron.SuffixTree (this repository, commit 24231842)"; the title and the intro still say 2026-11-03; the sentence "The functional lead is therefore thin and depends on how the library boundary is drawn" is hard-coded
and is now less true (tree layer only 56.8 % against sdsl 52.5 %, whole library 72.8 %); the sentence that lists the overridden cells names "C10-C13, C15, C16" and is still right (C14 is in `K0_audit.tsv`).

## 6. README headline numbers and the rows they come from

`README.md` of this directory, "Headline" table and the paragraph "Measured at 10 Mbp". All figures are the median [min-max] of 5 fresh processes unless stated; 10 Mbp random DNA `dna_10m`.

| figure in the README | value | row |
|---|---|---|
| property mean, function coverage, tree layer only | 2.28, 72.8 %, 56.8 % | regenerated `COMPETITOR_MATRIX.md`, section 2 |
| in-memory build SA + LCP + child table | 1.138 s [1.118-1.174] | `analysis/bench/results_1b9_rows.tsv:9` (`build_total_s`), ADR-004 1b.9 A1 |
| the public `EsaSuffixTreeFactory.CreateEsa` | 0.841 s [0.832-0.907] | `results_1b9_rows.tsv:44` (`entry_build_s`), ADR-004 1b.9 A1 |
| libsais, same session | 0.344 s (pooled median of 15 runs; points 0.3385, 0.3454, 0.3126 s) | `analysis/bench/results_1b9_libsais.tsv` (`build time`), ADR-004 1b.9 "Machine-drift control" |
| Ukkonen tree build (old) | 14.6 s | `BENCHMARK.md` section 1 (unchanged) |
| index size in memory | 10.00 B/char | `results_1b9_rows.tsv:22` |
| container Core / Full | 9.25 / 13.25 B/char (file 13.2515) | `results_1b9_sections.txt:11-13`, ADR-004 1b.9 A2 |
| peak RSS of the building process | 15.38 B/char (VmHWM 153.8 MB) | `results_1b9_rows.tsv:12` |
| count of present 20-mers, 1 thread | 647,176 q/s [611,792-824,474] | `results_1b9_rows.tsv:202`, ADR-004 1b.9 A3 |
| locate of 1,000 8-mers, 1 thread | 990,228 pat/s [876,565-1,627,933]; lazy `Locate` 2,173,239 | `results_1b9_rows.tsv:205-206`, A4 |
| 4 threads on one tree (count present / locate) | 3,098,222 [2,646,911-3,631,238] / 4,841,674 [4,371,764-5,279,213] q/s; its own one-thread figures on the same instrument 783,122 / 1,022,181 | `results_1b9_abapp.txt` through `python3 -I summarize_1b9.py --gate results_1b9_abapp.txt`, ADR-004 1b.9 A3-A4 |
| 100 Mbp (10^8) container, 4 threads | 17.88 s [15.88-19.43], RssAnon 253.2 MB [247.1-256.1], 1,325,019,393 bytes (13.25 B/char); genome-like 20.86 s, 286.3 MB, 13.49 B/char | `results_1b9_chain.txt:9-10`, ADR-004 1b.9 table C |
| not met | first `Open` of a fresh process 15.24 ms [13.09-15.37] against the 10 ms of 1b.7 (8.5 ms with ReadyToRun); protein count on the file-backed tree 0.934 of the in-memory tree at four threads | ADR-004 1b.9 A9 (`results_1b9_rows.tsv:527`, `warm_open` = the page-cache-warm first call, 0.015237 s) and the `abapp` gate table, `results_1b9_open_r2r.txt` |

The row "first place over 2000 random weightings" of the headline table is a figure of `plan_model.py`, which this scorer does not read; it was marked "not recomputed" instead of being carried as a current number; **the orchestrator recomputed it afterwards (2.9 % of 2000 weightings, same method as `plan_model.py`: seed 7, weights 1-3 per property, run on the regenerated matrix; worst gap -0.76, median -0.26, best +0.23) and the README carries that figure.**

## 7. Defects and questions (reported, not fixed; outside the files of this package)

1. **`analysis/merge.py` reads `results/K*.tsv`, and `results/K0.tsv` is a copy of the old `audit/K0_audit.tsv`** (identical before this package). After the re-score `merge.py` shows Seqeron 2.20 (it sees the five new rows of
   `bench_overrides.tsv` but not C8 and C14 of `K0_audit.tsv`); the matrix (`build_docs.py`, which reads `audit/K0_audit.tsv`) shows 2.28. `merge_func.py` agrees with the matrix (72.8 %). Fix: read `audit/K0_audit.tsv`
   in `merge.py`, or let the owner allow the new rows in `results/K0.tsv`.
2. **`build_docs.py` hard-codes** the K0 version line (commit 24231842, version 1.0.0), the date 2026-11-03 and the "thin lead" sentence (section 5 above). The Seqeron column now scores the repository at the head of 2026-11-23.
   Question for the owner: may `build_docs.py` get a K0 line "this repository, head 50936932 (production code 12495c3e)" and the lead sentence computed from the numbers?
3. **`SesaBuilder` and `SesaBuildOptions` are `internal`**: the external-memory builder, the one thing S1b was built for, has no public entry. `EsaSuffixTreeFactory.Open` reads a container but nothing public writes one; the
   builder is reachable from the Bench app, EmProto and the tests by `InternalsVisibleTo`. This is why F03, F70 and F79 stay at 1 and why C14 / C15 carry "builder is internal" in the evidence. A public `Build(input, target, options)` is the item that would lift them.
4. **The container holds byte texts only** (`SesaBuilder.Build`: "one byte per code unit"; the conformance factory `esa-file` accepts texts below 256 only); the v6 tree keeps UTF-16. Not stated in a user-facing document.
5. **There is no user documentation of the ESA backend**: `docs/algorithms/Pattern_Matching/Suffix_Tree.md`, `SuffixTree.Persistent/README.md` and the wiki describe the Ukkonen and v6 trees only; the container format lives in ADR-003, the usage in the Bench README and `S1B-RUNBOOK.md`. `analysis/audit/K0_func_extra.tsv` was not extended
   (not in the list of files of this package) and lacks `IndexLocus`, `ISuffixIndex`, `ISuffixTreeNavigation`, `ISuffixArrayAccess`, `EsaOpenOptions` (validation levels, access modes) and `SuffixIndexFactory.Open`.
6. The file-backed tree reads at 0.934 of the in-memory tree on the 20-letter text at four threads, and the first `Open` of a fresh process is 13-15 ms: both already open points of ADR-004 1b.9 (the orchestrator's decisions are there); no cell depends on them.
7. Not scored because only the owner's run can prove it: the production chain at 3x10^9 and 2x10^9, arm64, macOS, a 64-bit position type, S2 and later stages, the FM-index, both strands, matching statistics as a public API.

## 8. Cells I am unsure about

- **C10 = 3** (ratio 2.4-3.5x lies between the 3 of K6 at 2.2x and the 2 of K1 at 5.0x; the median of the gate row alone, 3.5x, is on the 2 side of the logarithmic midpoint 3.3x, the pooled and the user-facing figures are on the 3 side).
- **C14 = 3** (the capability is measured at 10^8 with bounded RAM, but the builder is not public API and the cold-cache query cost is not measured).
- **C15 = 2** (the 1 / 2 / 4 thread times are those of the prototype of the same pipeline; the production default is 4 threads and its 4-thread time is measured, its 1-thread time is not).
- **C12 = 4** (the best-tool rows are of 2026-11-03 and the Seqeron rows of 2026-11-23; the same-session libsais control shows no drift, +-5 %, and the margin to the 1.5x line is 30 %).
- **F18 / F20 = 1** (`TryExtendLeft` is correct and public but flagged Degraded by the library itself; an argument for 2 exists).
- **C19 = 2, C20 = 3 and C23 = 3 kept** (the candidates for 3 / 4 weighed in section 3; C23 is the one where the volume of tests is far above level 3 and the level-4 definition has no measurement).


## 9. Corrections after the review of 1b.10 (2026-11-23, by the orchestrator)

The independent review of S1 (`1b10_review.md`, finding R05 and its section 4) found that two cells of this ledger contradict the ledger's own F03 (build from a file without loading the text, held at 1 because the builder is `internal` and the only route is the Bench / `s1b.sh` command line), and that the base of C10 was not stated. The orchestrator accepted the reasoning (the same doubt is written in section 8 above) and corrected the evidence files by appending new rows; the 1b.11 rows stay in the files and the later row wins, as in the other evidence files.

| cell | after 1b.11 | after the review | the row that now decides |
|---|---|---|---|
| C14 external / larger-than-RAM | 3 | **2** | the last `K0 C14` row of `audit/K0_audit.tsv`: build and query with bounded RAM are measured (RssAnon 251-286 MB at 10^8), but the builder is internal; level 3 needs a public builder |
| C15 parallel construction | 2 | **1** | the last `K0 C15` row of `audit/bench_overrides.tsv`: four worker threads by default, but only in the internal builder (Bench `--threads`); the 1 / 2 / 4 thread figures are the E6 prototype's; the public `CreateEsa` is single-threaded |
| C10 construction speed | 3 (V) | 3 (**U**) | the last `K0 C10` row of `audit/bench_overrides.tsv`: the sentence now names the base, 3.5x of libsais SA only (0.328 s) and 2.0x of libsais SA+LCP (0.564 s, ADR-002 amendments); the comparison base is not like for like, hence `U` |

Effect: weighted property mean **2.28 -> 2.20**, rank **4 of 13 unchanged** (ropebwt3 2.48, GenomeTools 2.30, BWA/Bowtie2 2.30, Seqeron 2.20, r-index 2.13); function coverage unchanged (72.8 %, 56.8 % on the tree layer). First place over 2000 random weightings (the method of `plan_model.py`, seed 7, weights 1-3): **2.9 % -> 0.4 %** (worst gap -0.88, median -0.34, best +0.15). Diff of the regenerated matrix against the one before the correction: 5 line pairs, all in the Seqeron column (C10 flag, C14, C15, the weighted-mean cell) and the ranking line; no competitor cell changed.

**What raises them again:** a public builder entry (or the Bench command shipped and documented as the tool, as the competitors' command lines are scored): C14 -> 3, C15 -> 2 and F03 -> 2 move together, because they stand on the same builder. C10 stays 3 until an SA+LCP build within 1.5x of the best tool is measured on the public path.
