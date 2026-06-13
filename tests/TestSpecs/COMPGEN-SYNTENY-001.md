# Test Specification: COMPGEN-SYNTENY-001

**Test Unit ID:** COMPGEN-SYNTENY-001
**Area:** Comparative Genomics
**Algorithm:** Synteny / Collinearity Block Detection (MCScanX collinearity model)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wang et al. (2012) MCScanX, Nucleic Acids Res. 40(7):e49 (PMC full text) | 1 / 3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC3326336 | 2026-06-13 |
| 2 | Wang et al. (2012) MCScanX (Oxford Academic) | 1 | https://academic.oup.com/nar/article/40/7/e49/1202057 | 2026-06-13 |
| 3 | Wikipedia — Synteny (definitions + primaries) | 4 | https://en.wikipedia.org/wiki/Synteny | 2026-06-13 |

### 1.2 Key Evidence Points

1. DP recurrence `Score(v) = max(MatchScore(v), max(Score(u) + MatchScore(v) + GapPenalty × NumberofGaps(u,v)))` — Wang et al. (2012), PMC3326336.
2. `MatchScore = 50` per gene pair; `GapPenalty = −1`; `NumberofGaps < 25` — PMC3326336.
3. Reported blocks have score > 250, i.e. **at least 5 collinear gene pairs** — PMC3326336.
4. Both transcriptional directions → forward and inverted (reverse) blocks — PMC3326336.
5. "Collinearity, a more specific form of synteny, requires conserved gene order." — Oxford Academic HTML.
6. Anchors are homologous gene pairs (orthologs/paralogs); E-value cutoff 1e-5 for anchor generation (out of scope here) — PMC3326336.

### 1.3 Documented Corner Cases

- Chains scoring ≤ 250 (< 5 pairs) are not reported.
- An anchor ≥ 25 intervening genes from the previous anchor cannot extend the chain.
- Only non-overlapping chains are reported.
- Decreasing target order → inverted block (valid).
- No anchors → no blocks.

### 1.4 Known Failure Modes / Pitfalls

1. Mixing increasing and decreasing target order in one chain (direction must be consistent) — PMC3326336 (separate sorting per direction).
2. Reporting blocks below the 5-pair / score-250 minimum — PMC3326336.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindSyntenicBlocks(genome1Genes, genome2Genes, orthologMap, minAnchors, maxGap)` | ComparativeGenomics | Canonical | DP collinearity chaining with MCScanX scoring |
| `VisualizeSynteny(blocks)` | ComparativeGenomics | Delegate | Text rendering of blocks (smoke only) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported block has `GeneCount ≥ MinAnchors` (≥ 5 with defaults). | Yes | PMC3326336 (≥5 pairs / score 250) |
| INV-2 | Every reported block has a DP chain score ≥ `MinChainScore` (250 default). | Yes | PMC3326336 (scores over 250) |
| INV-3 | Within a block, consecutive anchors differ in genome2 position by `1 ≤ |Δ| ≤ maxGap+1` and direction is consistent (all increasing OR all decreasing). | Yes | PMC3326336 (NumberofGaps < 25; both directions) |
| INV-4 | A block is `IsInverted = true` iff the genome2 order across its anchors is decreasing. | Yes | PMC3326336 (both transcriptional directions) |
| INV-5 | Block coordinates satisfy `Start1 ≤ End1`, `Start2 ≤ End2`, and lie within the parent genes' coordinate spans. | Yes | Coordinate definition (implementation contract) |
| INV-6 | Reported chains are non-overlapping: each anchor (genome1 gene) belongs to at most one block. | Yes | PMC3326336 (non-overlapping chains) |
| INV-7 | Empty genome list or no anchors → empty result (no exception). | Yes | Definition requires anchors |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Forward 5-anchor chain | 5 adjacent orthologs in same order both genomes | 1 block, GeneCount=5, IsInverted=false, score 250 | PMC3326336 (MatchScore 50, ≥5 pairs) |
| M2 | Reverse 5-anchor chain | 5 anchors, genome2 order reversed | 1 block, GeneCount=5, IsInverted=true | PMC3326336 (both directions) |
| M3 | Sub-threshold 4-anchor chain | 4 adjacent anchors (score 200) | 0 blocks | PMC3326336 (min 5 pairs / 250) |
| M4 | Gap exceeds cutoff | two 3-anchor runs separated by ≥ maxGap intervening genes in genome2 | chain breaks; no 5-pair block → 0 blocks | PMC3326336 (NumberofGaps < 25) |
| M5 | Empty genome | genome1 empty | 0 blocks, no exception | INV-7 |
| M6 | No orthologs | non-empty genomes, empty ortholog map | 0 blocks | Definition requires anchors |
| M7 | Two separated forward blocks | two 5-anchor runs, distinct genome regions | 2 non-overlapping forward blocks | PMC3326336 (non-overlapping chains) |
| M8 | Coordinate bounds | forward block coordinates | Start1=min gene Start, End1=max gene End, Start2≤End2 | INV-5 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Gap-penalised long chain | 6 anchors with one 1-gene gap (NumberofGaps=1) | 1 block, GeneCount=6 (score 6×50−1=299≥250) | GapPenalty=−1 reduces but keeps chain |
| S2 | Null genome1 argument | genome1 = null | ArgumentNullException | Input validation |
| S3 | Ortholog points to absent gene | orthologMap target id not in genome2 | anchor skipped; no crash | Robustness |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | VisualizeSynteny smoke | render 1 forward block | non-empty string containing genome ids and gene count | Delegate |
| C2 | Property: all blocks valid | random-free constructed input | every block GeneCount≥5 and coords within parent bounds | O(n²) invariant property |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomicsTests.cs` — contains a `#region FindSyntenicBlocks Tests` with 4 tests (`FindSyntenicBlocks_CollinearGenes_ReturnsSyntenicBlock`, `_EmptyGenome_ReturnsEmpty`, `_NoOrthologs_ReturnsEmpty`, `_InvertedBlock_MarksAsInverted`). These use `minBlockSize:3`, assert only `Has.Count`/existence and a single boolean, lack exact GeneCount / score / coordinate assertions, and lack assertion messages — ⚠ Weak. They also encode the old non-conformant default (`minBlockSize=3`).
- A separate unit `ChromosomeAnalyzer_Synteny_Tests.cs` (CHROM-SYNT-001) tests a different class (`ChromosomeAnalyzer.FindSyntenyBlocks`); out of scope here.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 Forward 5-anchor | ❌ Missing | old test used 3 collinear genes, weak asserts |
| M2 Reverse 5-anchor | ⚠ Weak | `_InvertedBlock_MarksAsInverted` checks bool only, 3 genes |
| M3 Sub-threshold | ❌ Missing | no 4-anchor negative case |
| M4 Gap cutoff | ❌ Missing | none |
| M5 Empty genome | ⚠ Weak | `_EmptyGenome_ReturnsEmpty` exists, no message |
| M6 No orthologs | ⚠ Weak | `_NoOrthologs_ReturnsEmpty` exists, no message |
| M7 Two blocks | ❌ Missing | none |
| M8 Coordinate bounds | ❌ Missing | none |
| S1 Gap-penalised chain | ❌ Missing | none |
| S2 Null genome1 | ❌ Missing | none |
| S3 Ortholog→absent gene | ❌ Missing | none |
| C1 VisualizeSynteny | ❌ Missing | method did not exist |
| C2 Property | ❌ Missing | none |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_FindSyntenicBlocks_Tests.cs` — all COMPGEN-SYNTENY-001 cases (M/S/C) with exact evidence-based values.
- **Remove:** the `#region FindSyntenicBlocks Tests` block in `ComparativeGenomicsTests.cs` (4 weak tests superseded by the canonical file).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ComparativeGenomics_FindSyntenicBlocks_Tests.cs` | Canonical for COMPGEN-SYNTENY-001 | 13 |
| `ComparativeGenomicsTests.cs` | Other ComparativeGenomics methods (synteny region removed) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented forward 5-anchor test | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote inverted block test with exact values | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented 4-anchor sub-threshold test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented gap-cutoff test | ✅ Done |
| 5 | M5 | ⚠ Weak | Rewrote empty-genome test with message | ✅ Done |
| 6 | M6 | ⚠ Weak | Rewrote no-ortholog test with message | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented two-block test | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented coordinate-bounds test | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented gap-penalised chain test | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented null-argument test | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented ortholog→absent-gene test | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented VisualizeSynteny smoke test | ✅ Done |
| 13 | C2 | ❌ Missing | Implemented property test | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Forward 5-anchor → 1 block, GeneCount 5 |
| M2 | ✅ | Reverse 5-anchor → IsInverted true |
| M3 | ✅ | 4 anchors → 0 blocks |
| M4 | ✅ | Gap ≥ cutoff → 0 blocks |
| M5 | ✅ | Empty genome → empty |
| M6 | ✅ | No orthologs → empty |
| M7 | ✅ | Two runs → 2 blocks |
| M8 | ✅ | Coordinates within parent bounds |
| S1 | ✅ | 6-anchor gapped chain → 1 block of 6 |
| S2 | ✅ | Null → ArgumentNullException |
| S3 | ✅ | Dangling ortholog skipped |
| C1 | ✅ | VisualizeSynteny smoke |
| C2 | ✅ | Property holds for all blocks |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Report rule operationalised as score ≥ 250 AND anchors ≥ 5 (resolves "over 250" vs "≥5 pairs" wording; source-backed). | M1, M3, INV-1, INV-2 |
| 2 | Anchors supplied via `orthologMap` (anchor generation = COMPGEN-ORTHO-001, out of scope); collinearity algorithm unaffected. | All MUST tests |

---

## 7. Open Questions / Decisions

1. None — MCScanX provides the scoring scheme and all numeric defaults; both are source-backed.
