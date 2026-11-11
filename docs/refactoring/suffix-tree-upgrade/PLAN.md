# Suffix-tree upgrade plan (ST-UPGRADE) - proposed, NOT approved

Status: **proposal, 2026-11-03.** Nothing in this plan has been started. Do not start any stage without an explicit yes from the
owner (rule in `docs/Validation/review-2026-09/QUEUE.md`). This document replaces the one-line mention "ST-UPGRADE, stages 0-8"
there; the original stage text was not found in the repository or its history, so the stages below are derived from the
competitor analysis and their numbering is new.
(The earlier `suffix-tree-controlled-refactoring-plan.md` is an unrelated, finished internal refactoring.)

Evidence: [COMPETITOR_MATRIX.md](COMPETITOR_MATRIX.md) (13 tools x 24 properties and x 81 functions, sourced per cell),
[BENCHMARK.md](BENCHMARK.md) (measurements on one machine), `analysis/` (data, scripts, hand-off package).

## 1. Goal and success criterion

Owner's goal: the Seqeron suffix tree should be, on average, better than every competitor.

Measurable form, computed by `analysis/plan_model.py` from the matrices:

| criterion | now | target | competitor to beat |
|---|---|---|---|
| property mean, 24 properties, weighted 0-4 | 1.80 | >= best + 0.4 (>= 2.88) | ropebwt3 2.48 |
| function coverage, 81 functions | 55.6% (39.5% tree layer only) | >= best + 15 pp (>= 67.5%) | sdsl 52.5% |
| share of 2000 random weightings in which Seqeron is first | 0% | >= 95% | - |

Honest reading of "now": Seqeron is first only in function coverage, and that lead depends on Seqeron.Genomics scans
(motif, approximate, repeat finders) that are not on the tree. On the properties it is 9th of 13. The measured gaps
are build speed (45x slower than libsais at 10 Mbp, 186x at 50 Mbp), memory (169 B/char in memory, 56 B/char persistent file,
against 4 B/char for a suffix array), locate throughput (23x slower than rust-bio), scale (50 Mbp finished, 100 Mbp did not in
600 s, text length is `int`), no parallel build, no compression, no both-strand / approximate / IUPAC search on the index,
no public navigation API. Query count latency is on par (1.3x of the best in memory).

## 2. Stages

Cost is relative effort, not calendar time: S = 1, M = 3, L = 8. `dep` = needs the listed stage first.

| stage | name | cost | dep | main cells moved (before -> after, if its gate passes) |
|---|---|---|---|---|
| S0 | Safety, honest baseline, perf harness | 1 | - | none by design (removes defects, makes later claims measurable) |
| S1 | Enhanced-suffix-array (ESA) backend behind `ISuffixTree` | 8 | - | C8 2->3, C10 1->2, C11 1->2, C12 2->3; about 20 functions (SA/LCP/child table, navigation, handles, lazy locate) |
| S2 | Genome scale | 3 | S1 | C13 1->3 (needs external run), C19 2->3; 2-bit DNA, integer alphabet, file/stream build |
| S3 | Both strands, matching statistics, seeds, collections, batch | 3 | S1 | C5 2->3, C6 2->3, C7 1->3; ms[], SMEM, chaining, document mapping, batch queries |
| S4 | Approximate and degenerate search on the index | 3 | S1 | C2 1->3, C3 1->3; mismatch, edit, IUPAC, N, bidirectional via reversed-text index |
| S5 | Parallel construction | 3 | S1 | C15 0->2, C10 2->3 |
| S6 | Online append, merge, segment delete | 1 | - | C17 0->2 |
| S7 | Compressed FM-index | 8 | S1 | C18 0->3, C11 2->4; backward step, bidirectional, sampled locate |
| S8 | Packaging and proof | 1 | - | C22 1->3, C24 2->3, C23 3->4; NuGet, dotnet tool, MCP sessions, published benchmark |

Total 31 units; the recommended core (everything except S7) is 23.

### S0 - Safety, honest baseline, perf harness (cost 1)

Scope (all defects found by the audit of the current code, `analysis/audit/K0_gaps.md`):
- Correct the documentation claim "20 bytes per character" (`SuffixTree.cs:76`); measured 169 B/char in memory, 56 B/char persistent file.
- `PersistentSuffixTree.Load` validates only the header; query loops use raw pointers with `Debug.Assert`-only bounds. Add a
  structural validation pass (opt-in strict mode) and a documented trust boundary; hash/verify on request for v6 files.
- Guard int32 limits with explicit exceptions: `_nodeCount++` unchecked, `NodeCount` overflow above ~1.07 Gchar with a custom
  `ITextSource`, `HeapStorageProvider` 2 GB cap, `LoadStringInternal` above `int.MaxValue`, `PrintTree` capacity `estimatedNodes * 50` overflow above ~21 M nodes.
- `GenomicAnalyzer.FindRepeats` uses the O(n^2) `GetAllSuffixes`; switch to a tree/LCP-based method.
- MCP `suffix_tree_*` tools rebuild the tree on every call and have no size limit: add a documented limit now, sessions in S8.
- Commit the benchmark harness (`analysis/bench`) as a repeatable perf project (not run in CI by default) with a regression
  threshold of +-15% against the numbers in `BENCHMARK.md`.

Oracle and gate: existing suffix-tree test suites unchanged and green (the three `dotnet test` gates of the earlier refactoring plan);
each defect gets a failing test first; harness reproduces the 10 Mbp row of `BENCHMARK.md` within 15%.

### S1 - Enhanced-suffix-array backend (cost 8): the key stage

Why: it is the largest single gain (+18.5 pp functions, +0.20 property mean) and every other large stage depends on it. A suffix
array plus LCP plus child table is a suffix tree in the sense of Abouelhoda, Kurtz, Ohlebusch (2004); nodes are LCP intervals.

Scope:
- Linear-time SA construction (SA-IS) in managed code, Kasai/Phi LCP, child table (up/down/next index), optional bucket table for the first k characters.
- New backend implementing the existing `ISuffixTree` and the existing generic algorithms `SuffixTreeAlgorithms.*<TNode, TNav>` through an LCP-interval navigator, so
  `LongestRepeatedSubstring`, LCS, MEM, MUM, anchors, repeats are reused rather than rewritten.
- Public navigation API (child by character, parent, string depth, leaf rank/label, subtree size, suffix link, LCA/LCP(i,j), edge label), public pattern handle (interval) and
  lazy / limited locate, `SA[i]`, `LCP[i]`, inverse SA accessors.
- Factory selecting the backend; Ukkonen tree stays (online mode, differential oracle). Persistent container for the ESA arrays (mmap).

Oracle (oracle-first, as in the validation campaign): the existing Ukkonen tree on every operation; brute force for n <= 2000; SA equality against libsais on generated vectors
(`analysis/bench` inputs, hashes in `MANIFEST.sha256`); adversarial strings (runs, Fibonacci, Thue-Morse, de Bruijn, all-equal, binary and DNA/protein alphabets).

Acceptance gate (same harness and machine class as `BENCHMARK.md`, 10 Mbp random DNA, 1 thread). If it fails, stop: the projection in section 4 is void.
- Build (SA + LCP + child table) <= 3.3 s (<= 10x libsais 0.328 s). Expected 1-2 s (Go's SA-IS alone: 0.73 s).
- Index size <= 16 B/char including text (expected ~10: SA 4 + LCP ~1.25 + child 4 + text 1).
- Count 20-mers >= 450k q/s; locate of 1000 random 8-mers >= 250k patterns/s (geometric-mean ratio to the best measured tools <= 3x, which is what scores 3 on C12).
- All Ukkonen-oracle differential tests identical.

Risk: managed SA-IS is typically 2-4x slower than libsais; the gate allows 10x. Locate returns a contiguous SA slice, which removes the current subtree walk (7.6 M positions/s vs 71 M/s for Go).

### S2 - Genome scale (cost 3)

Scope: 64-bit position type on the ESA API (add `long` overloads, keep `int` API; breaking-change decision for the owner, recommended: additive), 2-bit packed DNA text and sigma=4 fast path,
integer alphabet, format v7 (v6 stays readable), build from FASTA/stream into a memory-mapped 2-bit text without a managed string.
Gate here (15 GB RAM, 4 vCPU): 100 Mbp builds in <= 60 s and <= 12 B/char; 500 Mbp builds and answers queries correctly (checksums against Go/libsais SA on the same input).
**Hand-off package required:** C13 level 3 means documented working at >= 3 Gbp (human genome), needs >= 64 GB RAM or a semi-external build; per the campaign rule this becomes a
package with scripts and expected checksums for another machine, result returned to this repository.

### S3 - Both strands, matching statistics, seeds, collections, batch (cost 3)

Scope: public `ms[]` (matching statistics) API; SMEM and seed generation; reverse-complement search done on the query side (no 2x index memory) for count, locate, MEM, MUM, seeds, with
coordinate mapping; generalized index over many texts with document id/offset mapping, document listing/counting, k-common substrings; public nucmer-style chaining of matches;
batch query API (many patterns per call, parallel).
Oracle: definitions by brute force (`...ByDefinition` implementations already exist for MEM/MUM/LCS), MUMmer 3.23 parity already established for forward strand, extended to reverse strand
(`mummer -b`/`-r`). Gate: all oracle tests identical; batch of 10^5 queries scales >= 3.5x on 4 threads.

### S4 - Approximate and degenerate search on the index (cost 3)

Scope: k-mismatch and k-edit search by pruned descent over the ESA intervals; IUPAC and N in the pattern as per-position alternation sets; policy for N in the indexed text;
bidirectional extension and search schemes (pigeonhole / optimal schemes) using an ESA of the reversed text.
Oracle: the existing `ApproximateMatcher.FindWithMismatches` / `FindWithEdits` and `MotifFinder.FindDegenerateMotif` scans (differential: identical hit sets on random and adversarial inputs, k <= 3).
Gate: identical to oracle; on 1 Mbp DNA, 10^4 random 30-mers, k = 2 faster than the O(n*m) scan by >= 100x and compared on the same machine with SeqAn3 (open source, header-only; scoring leader on C3).

### S5 - Parallel construction (cost 3)

Scope: parallel bucket partition + per-bucket sorting for SA, parallel LCP and child table.
Gate here (4 vCPU): 4-thread build of the ESA at 10 Mbp <= 1.0 s (3x libsais serial 0.328 s) and >= 2x faster than 1 thread; results byte-identical to the serial build.
Honest limit: libsais itself reaches only 1.8x on 4 threads at 10 Mbp. **Hand-off:** level 3 on C15 needs a documented near-linear speed-up on >= 8 threads; not verifiable on this machine.

### S6 - Online append, merge, segment delete (cost 1)

Scope: keep the Ukkonen builder as the online tail segment, immutable ESA segments behind it; queries cover all segments plus the m-1 boundary overlap; compaction rebuilds; tombstones by
document id; merge of two indexes as rebuild of the concatenation.
Oracle: rebuild-from-scratch equality of all query results after random append/delete/merge sequences.

### S7 - Compressed FM-index (cost 8)

Scope: BWT from the ESA, wavelet matrix / RRR rank-select (2-level for DNA), sampled suffix array with rate s and sampled inverse for extract, backward search (left step), locate with sampling,
bidirectional via a reversed-text index, persistent format.
Gate (10 Mbp random DNA, the hardest case for compression): <= 1.0 B/char including samples at the default rate (sdsl: 0.27-0.56); count >= 220k q/s (1/3 of sdsl `csa_wt_huff` 660k);
results identical to the ESA backend; locate speed documented per sampling rate.
This is the only stage that closes the compression/memory gap to ropebwt3, r-index and sdsl; without it C18 stays 0 and C11 stays 2.

### S8 - Packaging and proof (cost 1)

Scope: NuGet packages (`SuffixTree`, `SuffixTree.Persistent`), `dotnet tool` command line (`build`, `query`, `mem`, `stats`), MCP server sessions (build once, query many, persisted index handle), benchmark
results published in the README and regenerated by the harness (the re-score of the Seqeron column of the matrices is not here: it is work package 1b.11 at the end of S1, owner decision 2026-11-11), fuzzing job comparing backends, `CHANGELOG` and tags (today: no tags, version 1.0.0).
Gate: `dotnet pack` output installs and runs the sample; CI job runs the differential fuzz for a fixed budget.

## 3. Recommended order and decision points

S0 -> S1 -> **go/no-go on the S1 gate** -> S3 -> S4 -> S5 -> S2 -> S6 -> S8 -> decide S7 with the then-current numbers.

- If the S1 gate fails (build or memory off by more than 2x), stop and revisit: the projection assumes S1 succeeds.
- S7 is last because it is the most expensive (8) and the only stage the core result needs only for the +0.4 margin (see below).

## 4. Projection (`analysis/plan_model.py`, output in `analysis/plan_model_output.txt`)

The model sets Seqeron's cells to the score a strict auditor would grant if the stage gate passes (never lowering a cell) and recomputes both matrices with competitors unchanged.

| stages done | property mean | function coverage |
|---|---|---|
| now | 1.80 | 55.6% |
| S0+S1 | 2.00 | 74.1% |
| S0+S1+S3+S4 | 2.37 | 87.7% |
| S0+S1+S3+S4+S5+S8 | 2.61 | 88.3% |
| all except S7 (recommended core, cost 23) | 2.80 | 93.8% |
| all (cost 31) | 3.02 | 95.1% |
| best competitor | 2.48 (ropebwt3) | 52.5% (sdsl) |

Robustness to the weights (2000 random weightings, each weight 1-3): first place now 0%; S0+S1+S3+S4+S5+S8 77.7%; all except S7 99.5%; all 100%.
So "on average better than all" holds with the core set (+0.32 and +41 pp); the stricter +0.4 margin needs S7.

## 5. What the plan does not achieve

- Not best in every cell: after all stages Seqeron stays one level under the best tool on C3 (SeqAn3), C4 (MUMmer4), C5 (r-index), C7 (ropebwt3), C8 (sdsl), C10 (libsais), C12, C13 (ropebwt3), C18, and below on C14, C15, C17.
  The claim is "better on average", not "best at everything".
- A suffix tree, a suffix array and an FM-index do different jobs; build/size/count/locate are comparable, tree functions (LRS, MEM/MUM, navigation) are not.
- The benchmark is one 4-vCPU VM, random and lightly mutated DNA, warm cache, no out-of-core run; real genomes behave differently.
- 47 of 288 competitor property cells and 11 of 972 function cells are unverified (scored optimistically for the competitor, so the gap can only be larger than shown).
- Competitors move: Movi 2, Move-r, Columba/b-move, suffixient arrays and grlBWT were found by scouting and are not scored (`analysis/results/scout.md`).
- Targets are conditional on the gates; none is a promise.

## 6. Hand-off items (cannot be verified on the current 4-vCPU / 15 GB machine)

1. S2: >= 3 Gbp build/query (C13 level 3).
2. S5: >= 8 thread speed-up evidence (C15 level 3).
3. libsais64 / libdivsufsort at 1 Gbp, sdsl `cst_sct3` at >= 20 Mbp, BWA `-a is`, ropebwt3 / r-index / grlBWT on a repeat-rich multi-genome input (`analysis/bench/HANDOFF.md` section 3).
4. Out-of-core (index larger than RAM) run of the persistent tree.
Scripts: `analysis/bench/handoff_package.tar.gz`; results must not be merged into `bench.tsv` without a different host tag.

## 7. Decisions needed from the owner

1. Approve S0 and S1 (and thereby the go/no-go gate), or decline.
2. API policy for 64-bit positions (S2): additive `long` overloads (recommended) or breaking change.
3. Whether S7 (largest cost, needed only for the +0.4 margin) is in scope.
4. Where the hand-off runs happen (S2, S5), by whom.

## 8. Consumer requirements from the neoantigen mRNA vaccine pipeline (added 2026-11-10)

On 2026-11-10 the owner put the implementation of the missing vaccine-pipeline algorithms into the library plan
(`ALGORITHMS_ROADMAP.md`, review batches B27–B33; list in `docs/NEOANTIGEN_VACCINE_GAP_ANALYSIS.md`). Some of those units are
consumers of the index this plan builds. This section records what they need from it; **it makes no claim that the current
contract or the ESA meets any of it**, each row is checked at the named point. A row that fails is a finding for the stage, not
a reason to change the other plan.

| # | Need | Consumer unit | What to check | Checked at |
|---|------|---------------|---------------|-----------|
| R1 | Protein alphabet (about 25 symbols) and many short strings (a proteome is about 2×10^7 residues in tens of thousands of isoforms): exact lookup of 8–15-mers | NEO-SELF-001 | Does the contract (ADR-001) cover several strings with separators and σ above 4 in the public API, not only in unit tests? Measure build and lookup on a synthetic 20 Mresidue proteome | Before 1b.5 (contract reading), measured at gate 1b.9. **Contract reading done 2026-11-11 (documents only, nothing run):** symbols are UTF-16 code units 0..65535 (ADR-001 section 6, rule 1), so a 25-symbol protein alphabet is inside the contract; several strings are not a contract feature, the consumer joins the records with one separator symbol outside the alphabet (no peptide query contains it, so no match crosses a record) and maps positions back to records itself; no ADR change needed. Build and lookup on a 20 Mresidue proteome remain to be measured at 1b.9 |
| R2 | Backward search (BWT, rank checkpoints) for read seeds | IDX-FM-001, MAP-SEEDEXTEND-001 | Can the SESA container of ADR-003 carry or derive the BWT without a second build? Decide whether to reserve the field | Before 1b.7 (container), ADR-003 review |
| R3 | Seeds for short-read mapping: MEMs of a read against a 3 Gbp reference at read-batch throughput | MAP-SEEDEXTEND-001, HLA-TYPE-I-001 | The generic MEM path pays two suffix links and two rescans per query start; the one-locus path (candidate of S3) halves that. Measure reads per second on a large reference | S3 decision point |
| R4 | Approximate (k-mismatch) lookup of peptides and short reads on the index | NEO-SELF-001, MAP-SEEDEXTEND-001 | Which part is served by the index and which by `ApproximateMatcher` (B05) | 1b.6 facade review |
| R5 | A persistent index reused across patients (the proteome and the reference do not change) with fast load | NEO-SELF-001, MAP-* | Load time and resident size of the persistent form; builder must stream (ADR-004) | 1b.7 and E4 |
| R6 | Many concurrent readers over one index | every consumer | The navigation pack is built once and thread-safe (tested in 1b.4); the same for the facade | 1b.6 |
| R7 | Bounded-memory streaming of the sequencing files around the index | INFRA-STREAM-001 | Not an index property; the index must not require the whole genome string resident beyond what ADR-004 states | ADR-004 review |

Decision for the owner: whether the stage order of section 3 changes because of R1–R3 (for example, reading the contract for R1
before 1b.5). Until the owner says so, nothing in sections 2–3 changes.
