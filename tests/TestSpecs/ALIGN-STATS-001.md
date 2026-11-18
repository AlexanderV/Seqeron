# Test Specification: ALIGN-STATS-001

**Test Unit ID:** ALIGN-STATS-001
**Area:** Alignment
**Algorithm:** Pairwise Alignment Statistics (Identity / Similarity / Gaps) and Alignment Formatting
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-18 (B09 re-validation)

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Rice, Longden & Bleasby (2000), EMBOSS | 1 | https://doi.org/10.1016/S0168-9525(00)02024-2 | 2026-06-13 |
| 2 | EMBOSS needle docs (rel 6.6) | 3 | https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html | 2026-06-13 |
| 3 | EMBOSS Alignment Formats (markup legend) | 3 | https://emboss.sourceforge.net/docs/themes/AlignFormats.html | 2026-06-13 |
| 4 | NCBI BLAST QuickStart | 3 | https://www.ncbi.nlm.nih.gov/books/NBK1734/ | 2026-06-13 |
| 5 | pseqsid reference implementation | 3 | https://github.com/amaurypm/pseqsid | 2026-06-13 |
| 6 | EMBOSS 6.6 source `ajax/core/ajalign.c` (`alignConsStats`, `alignSim`) | 1 (reference implementation source) | https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/ajax/core/ajalign.c | 2026-11-18 |
| 7 | EMBOSS `needle` / `water` 6.6.0 binaries + `EBLOSUM62` / `EPAM250` / `EDNAFULL` data files (Debian `emboss`) | 1 | run locally | 2026-11-18 |
| 8 | Biopython 1.88 `Bio.Align.Alignment.counts`, `substitution_matrices` | 1 | PyPI `biopython==1.88` | 2026-11-18 |

### 1.2 Key Evidence Points

1. Identity = identical columns / alignment Length × 100; Length **includes gap columns**. — EMBOSS needle docs (source 2).
2. Similarity = (identical + similar) columns / Length × 100; a non-identical column is "similar" iff its substitution score is **positive**. — EMBOSS needle (source 2), BLAST "positives" (source 4), pseqsid (source 5).
3. Gaps = gap columns / Length × 100; a gap column is one where either aligned character is `-`. — EMBOSS needle (source 2), BLAST gap definition (source 4).
4. Published worked example (HBA vs HBB, EBLOSUM62): Length 149, Identity 65/149 = 43.6%, Similarity 90/149 = 60.4%, Gaps 9/149 = 6.0%. — EMBOSS needle (source 2).
5. srspair markup (needle/water 6.6.0 binary + `alignSim`, sources 6-7): `|` identity, `:` similar (positive score), `.` aligned but not similar, space gap column. The AlignFormats page (source 3) differs (`:` > 1.0, `.` small positive, space mismatch); binary and source win.
6. Column similarity per `alignConsStats` (source 6): `M(x,y) > 0` and a positive self score of x or y; an identical pair with a non-positive own score (X/X in BLOSUM62, N/N in EDNAFULL) is an identity, not similar, markup space. First-row ambiguity residues give the markup quirk (R/A `.`, A/R `:`; X/P space, P/X `.`).
7. Gap openings = maximal `-` runs per row (BLAST tabular `gapopen`; Biopython `open_gaps`) — source 8.

### 1.3 Documented Corner Cases

- Gap columns are counted in the percentage denominator (Length). — source 2.
- Similarity% ≥ Identity% always (similar set ⊇ identical set). — sources 2, 5.
- DNA simple match/mismatch model (Mismatch < 0): no non-identical column is similar ⇒ Similarity equals Identity. — sources 2, 4 + model definition.

### 1.4 Known Failure Modes / Pitfalls

1. Using `(matches + mismatches)/Length` as Similarity is incorrect — that is the non-gap fraction, not the EMBOSS/BLAST similarity (which requires a positive substitution score). — sources 2, 4.
2. Excluding gap columns from the denominator inflates all percentages. — source 2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateStatistics(AlignmentResult, ScoringMatrix?)` | SequenceAligner | Canonical | scalar model; Identity/Similarity/Gaps/Positives/GapOpens per EMBOSS needle |
| `CalculateStatistics(AlignmentResult, SubstitutionMatrix)` | SequenceAligner | Canonical | residue-pair model (BLOSUM62/PAM250/EDNAFULL) |
| `FormatAlignment(AlignmentResult, int, ScoringMatrix?)` | SequenceAligner | Canonical | needle/water three-line markup, scalar model |
| `FormatAlignment(AlignmentResult, int, SubstitutionMatrix)` | SequenceAligner | Canonical | markup per EMBOSS alignConsStats |
| `SubstitutionMatrix.Blosum62/Pam250/EdnaFull` | AlignmentTypes | Data | cell-for-cell vs EMBOSS data files + Biopython |
| MCP `alignment_statistics`, `format_alignment` | AlignmentTools | Wrapper | delegate; `matrix`/`match`/`mismatch` parameters |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Matches + Mismatches + Gaps = AlignmentLength | Yes | Column classification partitions the alignment (sources 2, 4) |
| INV-2 | Identity% ≤ Similarity% ≤ 100% | Yes | Similar set ⊇ identical set (sources 2, 5) |
| INV-3 | Identity% + (Similarity−Identity)% + (mismatch-only)% + Gap% = 100% | Yes | Length-denominator partition (source 2) |
| INV-4 | Markup line length per block = sequence block length; chars ∈ {`|`,`:`,`.`,space} | Yes | srspair (sources 6-7) |
| INV-5 | GapOpens ≤ Gaps; = 0 iff Gaps = 0 | Yes | run definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | EMBOSS formula/denominator | Counts 65 id, 9 gaps over Length 149 (similar=90) reproduced as percentages | Identity 43.6%, Similarity 60.4%, Gap 6.0% (Within 1e-3) | EMBOSS needle worked example (source 2) |
| M2 | DNA SimpleDna similarity = identity | Aligned DNA, default SimpleDna scoring | Similarity == Identity; similar count == match count | sources 2,4 + Mismatch<0 |
| M3 | Positive-scoring mismatch ⇒ similar | Mismatch=+1 scoring, one non-identical column | Identity 75%, Similarity 100% | EMBOSS/BLAST/pseqsid similarity rule |
| M4 | Gap counting | Column with `-` counts as gap, not identity/mismatch | Gaps counted, excluded from matches/mismatches | source 2,4 |
| M5 | Exact counts hand alignment | 9-col DNA alignment (6 id,1 mm,2 gap) | Matches 6, Mismatches 1, Gaps 2, Length 9, Id 66.6̅%, Sim 66.6̅%, Gap 22.2̅% | Direct EMBOSS-formula derivation |
| M6 | FormatAlignment legend (identity/gap/mismatch) | SimpleDna DNA alignment with mismatch and gap | markup `|` for id, `.` for mismatch, space for gap (`ACG-T`/`ACATT` → `||. |`) | needle 6.6.0 (sources 6-7) |
| M7 | FormatAlignment similarity mark | Mismatch=+1 scoring, similar column | markup `:` at similar column | EMBOSS AlignFormats (source 3) |

| M8 | Published EMBOSS example, real BLOSUM62 | HBA/HBB needle rows | Length 149, Matches 65, Positives 90, Gaps 9, GapOpens 4; markup lines of blocks 1-3 equal needle's | EMBOSS needle docs + binary (sources 2, 7) |
| M9 | 301-alignment oracle | needle/water 6.6.0 header + markup, Biopython 1.88 counts, DNA / IUPAC DNA (EDNAFULL) / protein (BLOSUM62, PAM250 incl. B/Z/X) | Length, Identity, Similarity, Gaps, mismatches, GapOpens and markup line of every alignment equal | sources 6-8 (`TestData/Alignment/needle_water_stats.tsv`) |
| M10 | Identical pair with non-positive score | `AXK`/`AXR` BLOSUM62; Match=0 scalar | Positives 2 (X/X not similar); markup `\| :`; scalar Positives 0 | source 6 |
| M11 | First-row ambiguity quirk | EDNAFULL R/A, A/R, Y/C, C/Y, R/C, N/A, A/N; BLOSUM62 X/P, P/X, X/X | markup `.` `:` `.` `:` space space `.`; space `.` space | needle 6.6.0 binary |
| M12 | Case-insensitive; unequal rows | `acgt` = `ACGT`; rows of unequal length | same stats; `ArgumentException` | EMBOSS; contract |
| M13 | GapOpens | `AC--GT`/`AGTTGT` = 1; two-row runs; affine alignment with one long gap = 1 | counts of maximal runs | Biopython `open_gaps`, BLAST gapopen |
| M14 | Matrix tables | spot values + symmetry + alphabets; scalar presets (`BlastDna` 2/−3/−5/−2 = blastn defaults, `HighIdentityDna` match/mismatch = EDNAFULL ACGT) | published values | EMBOSS data files, Biopython, BLAST+ manual |
| M15 | MCP binding | `alignment_statistics(..., "blosum62")` = published example; scalar `mismatch = 1`; unknown matrix rejected; `format_alignment` with matrix `::.: |` | as listed | tool contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty alignment stats | AlignmentResult.Empty | AlignmentStatistics.Empty (all 0) | documented empty handling |
| S2 | Empty alignment format | AlignmentResult.Empty | "" | documented empty handling |
| S3 | Null alignment | CalculateStatistics/FormatAlignment(null) | ArgumentNullException | input validation |
| S4 | Non-positive lineWidth | FormatAlignment(.., 0) and (.., -1) | ArgumentOutOfRangeException | input validation |
| S5 | Perfect identity | identical aligned sequences | Identity 100%, Similarity 100%, Gap 0% | trivial boundary |
| S6 | All gaps | every column has a `-` | Gaps == Length, Identity 0% | gap boundary |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Line wrapping | lineWidth smaller than alignment ⇒ multiple blocks | block count = ceil(Length/lineWidth); each markup line correct length | INV-4 |
| C2 | Invariant property | Random-but-deterministic alignment columns | INV-1, INV-2, INV-3 hold | property-based |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `SequenceAligner*`: existing files cover GlobalAlign, LocalAlign, SemiGlobalAlign, MultipleAlign and a legacy `SequenceAlignerTests.cs`. No dedicated test file for `CalculateStatistics` / `FormatAlignment`.
- `grep` for `CalculateStatistics`/`FormatAlignment` in test tree: only incidental usage in the legacy file (no exact evidence-based assertions on Similarity semantics).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new |
| M2 | ❌ Missing | new |
| M3 | ❌ Missing | new |
| M4 | ❌ Missing | new |
| M5 | ❌ Missing | new |
| M6 | ❌ Missing | new |
| M7 | ❌ Missing | new |
| S1 | ❌ Missing | new |
| S2 | ❌ Missing | new |
| S3 | ❌ Missing | new |
| S4 | ❌ Missing | new |
| S5 | ❌ Missing | new |
| S6 | ❌ Missing | new |
| C1 | ❌ Missing | new |
| C2 | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAligner_CalculateStatistics_Tests.cs` — all ALIGN-STATS-001 cases (both methods, one canonical file per unit).
- **Remove:** nothing; legacy `SequenceAlignerTests.cs` is out of scope for this unit and untouched.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `Unit/Alignment/SequenceAligner_CalculateStatistics_Tests.cs` | Canonical (ALIGN-STATS-001) | 69 (incl. TestCase rows) |
| `Seqeron.Mcp.Alignment.Tests/AlignmentStatisticsTests.cs`, `FormatAlignmentTests.cs` | MCP binding | see file |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented | ✅ Done |
| 11 | S4 | ❌ Missing | Implemented | ✅ Done |
| 12 | S5 | ❌ Missing | Implemented | ✅ Done |
| 13 | S6 | ❌ Missing | Implemented | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented | ✅ Done |
| 15 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | CalculateStatistics_EmbossExample_MatchesPublishedPercentages |
| M2 | ✅ Covered | CalculateStatistics_SimpleDna_SimilarityEqualsIdentity |
| M3 | ✅ Covered | CalculateStatistics_PositiveScoringMismatch_SimilarityExceedsIdentity |
| M4 | ✅ Covered | CalculateStatistics_GapColumns_CountedAsGaps |
| M5 | ✅ Covered | CalculateStatistics_HandAlignment_ExactCountsAndPercentages |
| M6 | ✅ Covered | FormatAlignment_SimpleDna_UsesIdentityAndGapLegend |
| M7 | ✅ Covered | FormatAlignment_PositiveScoringMismatch_UsesSimilarityMark |
| S1 | ✅ Covered | CalculateStatistics_EmptyAlignment_ReturnsEmpty |
| S2 | ✅ Covered | FormatAlignment_EmptyAlignment_ReturnsEmptyString |
| S3 | ✅ Covered | CalculateStatistics_Null_Throws / FormatAlignment_Null_Throws |
| S4 | ✅ Covered | FormatAlignment_NonPositiveLineWidth_Throws |
| S5 | ✅ Covered | CalculateStatistics_PerfectIdentity_HundredPercent |
| S6 | ✅ Covered | CalculateStatistics_AllGaps_GapsEqualLength |
| C1 | ✅ Covered | FormatAlignment_NarrowLineWidth_WrapsIntoBlocks |
| C2 | ✅ Covered | CalculateStatistics_Invariants_Hold |
| M8 | ✅ Covered | CalculateStatistics_EmbossHbaHbb_Blosum62_ReproducesPublishedExample / FormatAlignment_EmbossHbaHbb_Blosum62_MatchesNeedleMarkupLine |
| M9 | ✅ Covered | CalculateStatistics_And_FormatAlignment_MatchEmbossAndBiopython_OnOracleTable |
| M10 | ✅ Covered | CalculateStatistics_IdenticalPairWithNonPositiveScore_IsNotSimilar / FormatAlignment_IdenticalPairWithNonPositiveScore_IsBlank |
| M11 | ✅ Covered | FormatAlignment_EdnaFull_FirstRowAmbiguityCode_FollowsEmbossQuirk / FormatAlignment_Blosum62_FirstRowX_FollowsEmbossQuirk |
| M12 | ✅ Covered | CalculateStatistics_LowerCaseRows_AreCaseInsensitive / CalculateStatistics_And_FormatAlignment_UnequalRows_Throw |
| M13 | ✅ Covered | gap-open tests incl. CalculateStatistics_OnAffineGlobalAlignment_SingleLongGap_OneOpening |
| M14 | ✅ Covered | Blosum62/Pam250/EdnaFull_PublishedValues, SubstitutionMatrix_Tables_AreSymmetric_AndRejectUnknownSymbols, ScalarPresets_MatchPublishedMatchMismatch |
| M15 | ✅ Covered | AlignmentStatisticsTests (Blosum62 example, scalar, case/unknown matrix), FormatAlignmentTests (matrix marks) |

---

## 6. Assumption Register

**Total assumptions:** 0 (the former assumption "srspair `.` tier unreachable / collapses to `:`" was wrong: needle 6.6.0 writes `.` for every aligned non-similar column and the identity/similar marks depend on self scores; implemented per `alignConsStats`).

---

## 7. Open Questions / Decisions

1. The pre-existing implementation computed Similarity as `(matches+mismatches)/Length` (corrected 2026-06 to the positive-score rule).
2. B09 re-validation (2026-11-18): oracle-confirmed defects fixed — case-sensitive comparison (F11), unequal rows accepted (F12), no gap-opening / positives output (F13), no residue-pair similarity and a markup that differed from needle (F14). Decision: markup follows the needle 6.6.0 binary/source rather than the AlignFormats web page where they differ; the first-row-ambiguity quirk is reproduced (documented, tested) because the unit's contract is EMBOSS parity.

## 2026-11 audit round 1 addition (A1-9, B09 F42/F43): BLAST integer percentages and multiple-alignment statistics

- `SequenceAligner_BlastStatistics_Tests` (BLAST+ 2.12.0 blastn/blastp, `TestData/Alignment/blast_stats_oracle.tsv`, 320 best HSPs: counts, `BlastPositives`, integer percentages, `pident` text, report line; `GetPercentMatch` table incl. the 99-not-100 cap).
- `MultipleAlignmentStatistics_Tests` (`TestData/Alignment/msa_stats_oracle.tsv`: EMBOSS infoalign 6.6.0 3 268 per-row records with explicit/consensus reference, distmat 6.6.0 600 matrices (uncorrected and Jukes-Cantor, gap weights 0 / 0.5), ClustalW 2.1 conservation lines (151 DNA/protein MSAs), Biopython `DistanceCalculator('identity')` 77 matrices) and `Seqeron.Mcp.Analysis.Tests.MultipleAlignmentStatisticsToolTests`; generator `scripts/alignment/generate_msa_stats_oracle.py`.
