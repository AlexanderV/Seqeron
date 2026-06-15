# Findings Register — full disposition of every note across the validation campaign

**Date:** 2026-06-12   **Library:** Seqeron.Genomics (mission-critical)
**Source:** every note/limitation/follow-up across all 86 reports in `docs/Validation/reports/`.

Each finding is placed in exactly one category:

| Category | Meaning |
|----------|---------|
| **FIXED-NOW** | Non-radical (doc/comment/spec text or small safe code change). Fixed in this pass. |
| **FEASIBLE → IMPLEMENT** | Needs real code change but is implementable safely with strict tests, in its own context. |
| **NOT-POSSIBLE (radical)** | Requires a redesign / public-API change / new model. Documented here, not changed. |
| **BY-DESIGN** | Note documents *correct, sourced* intended behaviour; "fixing" it would be wrong. No change. |

---

## A. FIXED-NOW — non-radical (doc/comment/spec) + small safe code fixes

| # | Unit | Finding | Fix |
|---|------|---------|-----|
| A1 | DISORDER-PRED-001 | Stale residue-ranking string `…Q,K,S,E,P` in Evidence doc & TestSpec (S<K, so order is `Q,S,K,E,P`). Code comment already fixed mid-campaign. | Correct the string in Evidence + TestSpec. |
| A2 | MIRNA-PRECURSOR-001 | Evidence "Dataset 3": hsa-mir-21 "8 consecutive pairs" → actually 16. | Correct Evidence text. |
| A3 | MIRNA-PRECURSOR-001 | INV-9 "all else equal" wording imprecise vs M11 (varies stem+loop). | Tighten spec wording. |
| A4 | MIRNA-PRECURSOR-001 | Stale checklist method names `FindPreMiRnas`/`ValidateHairpin` → real `FindPreMiRnaHairpins`. | Fix `ALGORITHMS_CHECKLIST_V2.md`. |
| A5 | SEQ-GCSKEW-001 | Stale method name `FindOriginOfReplication` → `PredictReplicationOrigin` in checklist line ~1543. | Fix `ALGORITHMS_CHECKLIST_V2.md`. |
| A6 | CHROM-KARYO-001 | `GetChromosomeBaseName` comment mentions "chr1a, chr1b" but code handles only `_N` integer suffixes. | Correct the code comment. |
| A7 | PROBE-DESIGN-001 | Stale inline comment "50-70 bp" in M11 test (params are 50-60); assertion already uses live params. | Correct the test comment. |
| A8 | PHYLO-DIST-001 | Evidence §5.3 prose rounds pure-transversion K2P to 0.31726; exact is 0.317162 (test tolerance fine). | Correct Evidence prose. |
| A9 | CRISPR-PAM-001 | Reverse-strand `PamSite.TargetStart` is a revComp-string index while `Position` is forward — internally consistent but easy to misread. | Add an XML `<remarks>` caveat on the field (no behaviour change). |
| A10 | META-CLASS-001 | Docs/spec conflate the implemented flat best-hit with Kraken/LCA (overclaim). | De-overclaim the description/XML docs/spec so it honestly states "best-hit, no LCA" (the LCA *implementation* is NOT-POSSIBLE → §C1). |
| A11 | SV-CNV-001 | `LogRatioToCopyNumber` rounded the absolute copy number `round(2·2^log2)` with `MidpointRounding.AwayFromZero`, diverging from the CNVkit reference (`do_call` uses NumPy `ndarray.round()` = round-half-to-even) at every exact half-integer copy number, e.g. copies=0.5→code CN1 vs CNVkit CN0; copies=2.5→3 vs 2. | Switch to `MidpointRounding.ToEven` to match NumPy/CNVkit (also aligns with sibling `CreateSegment`). Strengthened the S1 property test to the exact sourced sequence `[0,1,2,4,8]` (which exposed the bug) and added a half-integer round-half-to-even `[TestCase]`. Full suite 6493 green. |
| A11 | EPIGEN-METHYL-001 (Phase-2) | Test gaps only (no algorithm defect): CHH invalid-third-base branch (`CAN`→null) untested; `FindMethylationSites` positions test used a vacuously-true `All(...)`. | Added `GetMethylationContext_NonAcgtThirdBase_ReturnsNull`; guarded the positions test with `Is.Not.Empty`. Suite 6478 green. |
| A12 | EPIGEN-DMR-001 (Phase-2) | Test-quality only (no algorithm defect): M8 PValue test asserted only `InRange(0,1)` (a deliberately-wrong impl would pass); `FisherExactProbability` error/`n==0`/non-extreme branches and `AnnotateDMRs` null-input were untested. Algorithm itself matches scipy two-sided Fisher exactly across 7 tables and Wikipedia single-table 0.001346076. | Locked M8 + M1 PValues to scipy values (2.475428262210228e-39, 2.070073888186964e-35); added FisherExactProbability negative-cell throw, `n==0`→1.0, symmetric (5,5,5,5)→0.34371820130334063, and AnnotateDMRs null-input tests. Suite 6482 green. |
| A13 | VARIANT-ANNOT-001 (Phase-2) | Test-quality only (no algorithm defect): S1 (ambiguous codon) asserted only `Is.Not.EqualTo(SynonymousVariant)` — a deliberately-wrong impl could pass; public `GetConsequenceRank` had no direct exact-value coverage (only implicit via M13). Algorithm itself matches Constants.pm ranks/impacts and VariationEffect.pm predicates + NCBI transl_table 1 exactly (verified by hand-decode). | Strengthened S1 to assert exact `CodingSequenceVariant`/MODIFIER; added `GetConsequenceRank_MatchesConstantsPmRanks` (exact Constants.pm ranks) and `GetConsequenceRank_OrderingConsistentWithImpactTiers` (INV-5). Suite 6484 green. |
| — | CRISPR-OFF-001, SEQ-ENTROPY-001, PROTMOTIF-FIND-001, KMER-FIND-001 | Backwards comment / rounding typo / stale citation / weak asserts. | **Already fixed mid-campaign** (here for completeness). |
| A14 | SV-DETECT-001 (Phase-2) | **Real algorithm defect.** `IsConcordantOrientation` treated RF (reverse-forward, outward-facing / "everted") pairs as concordant, so an RF pair was never flagged discordant; `ClassifySV` had no Duplication branch. Per DELLY (Rausch et al. 2012), LUMPY, Manta and SVXplorer (Kumar et al. 2020) — and the unit's own cureffi/BWA source ("RF, FF or RR … that's a problem") — RF is the basic **tandem-duplication** signature for a short-insert FR library, hence discordant; only FR is the proper-pair (FLAG 0x02) orientation. Test S4 was green-washed (asserted RF "not discordant"). Description (SV_Detection.md), Evidence and TestSpec all repeated the wrong "FR/RF both concordant" claim. | Fixed `IsConcordantOrientation` → FR only; added RF→`SVType.Duplication` branch to `ClassifySV`; rewrote S4 to assert RF discordant→Duplication; added M9 ClassifySV Duplication test; corrected the description/Evidence/TestSpec and added DELLY + SVXplorer refs. Full suite 6485 green. |

## B. FEASIBLE → IMPLEMENT (own context, strict tests first, project patterns)

| # | Unit | Finding | Planned implementation |
|---|------|---------|------------------------|
| B1 | PROTMOTIF-PROSITE-001 | Extended ScanProsite `*` (Kleene-star) query metachar silently dropped → a pattern using it is silently mis-parsed. | **Reject** unsupported PROSITE metacharacters with a clear exception (mirror the Newick "throw, don't silently drop" fix). Strict tests. |
| B2 | POP-FST-001 | `CalculateFst` silently truncates to `min(count1,count2)` on mismatched locus counts. | Validate equal locus counts → throw `ArgumentException` on mismatch (defensive contract). Strict tests. |
| B3 | PARSE-EMBL-001 | INSDC doubled-quote escaping `""`→`"` not unescaped in qualifier values; plus a dead `QualifierRegex` path. | Unescape `""`→`"` in qualifier parsing; remove dead code. Strict tests with the doubled-quote case. |
| B4 | SPLICE-PREDICT-001 (LIMITED) | Region shorter than `minExonLength` dropped by `DeriveExons` but kept by `GenerateSplicedSequence` → coverage/spliced-length invariants (INV-3/4/11) inconsistent. | Make the two consistent (spliced sequence derives from the same reported-exon set), restoring all invariants while preserving `minExonLength` semantics and the `start<end` property test. Strict tests. |

## C. NOT-POSSIBLE (radical — documented, not changed in this pass)

| # | Unit | Gap | Why radical / what a real fix needs |
|---|------|-----|--------------------------------------|
| C1 | META-CLASS-001 | True taxonomic classification (Kraken weighted root-to-leaf / LCA). | Needs a taxonomy DAG data model, LCA-at-DB-build, per-read classification tree — changes the public API and the semantics of 27 locked tests. (Overclaim text is corrected in A10.) |
| C2 | Phylogenetics (PHYLO-NEWICK / TREE / COMP) | Full N-ary (multifurcating) trees. | `PhyloNode` is a binary `Left`/`Right` model; N-ary requires refactoring the node type and every consumer (UPGMA, NJ, RF, MRCA, Newick). Mitigated: Newick now **throws** instead of silently truncating (done mid-campaign). |
| C3 | PHYLO-COMP-001 | Unrooted-bipartition Robinson-Foulds as an alternative to the current rooted-clade metric. | Different metric; additive feature + new tests. Current rooted-clade RF is correct and documented. |
| C4 | ALIGN-MULTI-001 | Guide-tree progressive MSA. | Current star alignment is honestly declared; progressive MSA is a different algorithm (guide tree + profile-profile). |
| C5 | RESTR-DIGEST-001 | Circular-molecule digest. | Needs a topology parameter + wrap-around fragment joining; current linear digest is correct and documented. |
| C6 | ANNOT-GENE-001 | Reverse-strand ribosome-binding-site (SD) reporting in the RBS helper. | Feature addition (PredictGenes already covers both strands). |
| C7 | CRISPR-GUIDE/OFF-001, PRIMER-TM-001 | Doench/Azimuth on-target, MIT/Hsu off-target weight model, SantaLucia NN Tm. | New scientific models; current methods are honestly-declared heuristics, not defects. |
| C8 | PARSE-EMBL-001 | Remote references, deprecated single-dot ranges, dedicated site-between flag. | Rare/edge INSDC location features; out of standard-record scope. |

## D. BY-DESIGN — correct, sourced behaviour; intentionally unchanged

| Unit | Note | Why no change |
|------|------|---------------|
| PARSE-FASTA-001 | Orphan sequence (no header) / header with no sequence silently skipped. | Documented, matches Biopython lenient convention; "fixing" would break the sourced contract. |
| PARSE-FASTQ-001 | Negative-Q (Solexa) unsupported, clamps to Q0. | Spec scopes Phred+33/+64 only; Solexa obsolete. |
| PAT-EXACT-001 | Empty pattern → all positions `[0..n-1]`. | Defined via formal-language ε-substring convention. |
| SEQ-VALID-001 | Empty sequence → valid (true). | Only choice consistent with Biopython + universal-quantifier semantics. |
| ALIGN-SEMI-001 | Fitting-alignment variant (query fully aligned, reference ends free). | Sourced (Rosalind SIMS); correct semi-global variant. |
| ALIGN-MULTI-001 / RF / heuristics | Star MSA, rooted-clade RF, declared heuristics. | Honestly scoped; correct for their stated contract (see C3-C7 for optional upgrades). |
| RNA-STRUCT-001 | `PredictStructure` greedy dot-bracket vs Nussinov/Zuker DP traceback. | Documented approximation; the DP score values themselves are correct. |
| ASSEMBLY-STATS-001 | No defect. N50/L50/Nx/Lx inclusive boundary (`cumulative*100 >= total*threshold`) is algebraically identical to QUAST `s <= limit`; auN = Σl²/Σl; gaps = 0-based inclusive maximal N-runs. All values re-derived from live QUAST `N50.py` + Wikipedia worked examples (A→70/2, B→50/3) + lh3 auN blog. | Implementation and tests both correct; tests assert exact sourced values (no green-washing). Full suite 6497 green, 0 code/test changes. |
| SV-BREAKPOINT-001 | `FindBreakpoints` uses single-linkage on ADJACENT sorted junctions; INV-03 ("reported position within clusterTolerance of EVERY member") can fail on long chains where reported coord = rounded mean. | Single-linkage matches ClipCrop "sorted and clustered within 5-base differences" (adjacent diffs); reported coord is the explicit ASSUMPTION A1 (mean). INV-03 wording overstates the guarantee but the clustering+support logic is sourced and correct. Added chaining test to lock the actual behaviour. |
| ASSEMBLY-OLC-001 | `AssembleOLC` layout is greedy best-successor chaining (not transitive-reduction + Hamiltonian-path) and consensus is realized as chain concatenation (not per-column majority vote). | Exact OLC layout is NP-complete (Compeau et al. 2011); both simplifications are documented and sourced (Langmead OLC/SCS), and for exact-overlap chains concatenation == majority vote. Overlap detection (`FindOverlap`/`FindAllOverlaps`) is exactly correct — 12-edge GTACGTACGAT graph and all overlap/merge values re-derived by hand and matched. Test gap (MinContigLength discard branch had no test) fixed in session: added M5b with sourced exact values. Suite 6494 green. |
| ASSEMBLY-DBG-001 | `BuildDeBruijnGraph`/`AssembleDeBruijn` correct: (k-1)-mer multigraph, prefix→suffix edges, Hierholzer Eulerian walk, `p₀ + last-char(pᵢ)` spelling. No implementation defect — every published example (AAABBBA, a_long_long_long_time k=5, to_every… k=4, ATGGCGTGCA) re-derived from the **fetched** Langmead PDF + an independent Hierholzer reference and matched verbatim. Walk-selection among multiple Eulerian walks is a deterministic single choice (ASSUMPTION-1), so `to_every…` k=3 recovers the input under this walk (Langmead's printed mis-order is his arbitrary walk; C1 asserts only the branch node, correctly). | Test-coverage gap (3 documented branches untested) fixed in session: added `AssembleDeBruijn_DisconnectedGraph_OneContigPerComponent` (two disjoint reads → 2 contigs, Langmead p.24), `AssembleDeBruijn_MinContigLength_FiltersShortContigs` (length-10 contig dropped at MinContigLength=11), `BuildDeBruijnGraph_NullReads_ProducesEmptyGraph`. All values sourced/hand-derived, not code echoes. Fixture 17→20; full suite 6497 green. |

---

## Status

- A (FIXED-NOW): see per-row commits / this pass.
- B (FEASIBLE): implemented in separate contexts with strict tests — status tracked below.
- C (NOT-POSSIBLE): documented above; carried in the ledger's "Deferred BIG fixes" backlog.
- D (BY-DESIGN): no action.

| Item | State |
|------|-------|
| A1–A10 doc/comment fixes | ✅ DONE — all text/comment-only, 0 behaviour change, suite 4503 green. (A1 TestSpec had no such string — only Evidence; A8 corrected to exact 0.317128.) META-CLASS overclaim de-claimed in code XML docs + spec + Evidence. |
| B1 PROSITE reject `*` | ✅ DONE — throws FormatException on unsupported tokens, +6 tests, suite 4495 green, operators byte-for-byte unchanged |
| B2 FST mismatch throw | ✅ DONE — ArgumentException on mismatched locus counts (PairwiseFst inherits), +3 tests, suite 4498 green |
| B3 EMBL `""` unescape | ✅ DONE — INSDC `""`→`"` in EMBL **and** GenBank (old `.Trim('"')` was doubly wrong), dead QualifierRegex removed, +5 tests, suite 4503 green |
| B4 SPLICE-PREDICT invariants | ✅ DONE — consistent-filter fix, +3 strict tests, suite 4489 green, unit now CLEAN |
