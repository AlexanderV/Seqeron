# Test Specification: ANNOT-GENE-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | ANNOT-GENE-001 |
| **Area** | Annotation |
| **Algorithm** | Gene Prediction |
| **Canonical Methods** | `GenomeAnnotator.PredictGenes`, `GenomeAnnotator.FindRibosomeBindingSites` |
| **Complexity** | O(n) |
| **Status** | ☑ Complete |

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Gene prediction | Encyclopedia | ORF-based gene finding in prokaryotes, signal detection |
| Wikipedia: Shine-Dalgarno sequence | Encyclopedia | SD consensus AGGAGG, distance 4-15bp upstream of start |
| Wikipedia: Ribosome-binding site | Encyclopedia | RBS role in translation initiation |
| Shine & Dalgarno (1975) | Primary | Original SD sequence identification |
| Chen et al. (1994) | Primary | Optimal aligned spacing: 5 nt (E. coli); functional range 4-15 bp |
| Laursen et al. (2005) | Review | Bacterial translation initiation mechanisms |

## Method Specifications

### 1. PredictGenes(dna, minOrfLength, prefix)

**Description:** Predicts genes using ORF-based approach on both strands.

**Invariants:**
- All returned genes have Type = "CDS"
- All genes have strand '+' or '-'
- All genes have Start < End
- Gene IDs follow pattern "{prefix}_{number:D4}"
- Protein length in attributes matches (End-Start)/3 - 1 (excludes stop codon)
- At most one gene per (strand, stop codon), started at the first in-frame start after the previous stop (EMBOSS getorf -find 1; Prodigal one start per stop) — review 2026-09

### 2. FindRibosomeBindingSites(dna, upstreamWindow, minDistance, maxDistance)

**Description:** Finds Shine-Dalgarno sequences upstream of ORFs.

**Invariants:**
- Position is within upstream window of a valid ORF
- Sequence matches one of: AGGAGG, GGAGG, AGGAG, GGAG, GAGG, AGGA (all contiguous ≥4-nt substrings of AGGAGG; Prodigal shine_dalgarno_exact)
- Only maximal motifs per start (no hit contained in another in-range hit); each (position, motif) reported once
- Distance to start codon (aligned spacing) is within [minDistance, maxDistance]
- Score is normalized (motif.Length / 6.0) — declared heuristic: published RBS scores need genome-trained weights (Prodigal) or a thermodynamic model (RBS Calculator)

## Test Cases

### Must Tests (Evidence-Based)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | PredictGenes returns CDS type for all genes | Gene annotation standard | Implementation |
| M2 | PredictGenes assigns sequential gene IDs with prefix | Implementation contract | Implementation |
| M3 | PredictGenes includes strand info (+ or -) | Biological requirement | Wikipedia |
| M4 | PredictGenes filters by minOrfLength | Parameter contract | Implementation |
| M5 | PredictGenes finds genes on both strands | Both strands can encode genes | Wikipedia |
| M6 | FindRibosomeBindingSites detects AGGAGG consensus | SD consensus sequence | Shine & Dalgarno 1975 |
| M7 | FindRibosomeBindingSites validates distance constraints | Functional range 4-15 bp; optimal 5 nt | Chen et al. 1994, Wikipedia |
| M8 | Empty sequence returns empty result | Edge case | Implementation |
| M9 | Sequence without ORFs returns empty | No start/stop = no gene | Definition |
| M10 | PredictGenes protein length = (End-Start)/3 - 1 | Data integrity — excludes stop codon | Implementation, Biology |

### Should Tests (Recommended)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| S1 | Multiple genes in sequence all detected | Multi-gene operons common |
| S2 | Overlapping genes both reported | Can occur in different frames |
| S3 | Alternative start codons (GTG, TTG) recognized | Prokaryotic start codons |
| S4 | RBS shorter motifs (GAGG, AGGA) detected | Variant SD sequences |
| S5 | RBS score reflects motif length | Quality metric |

### Could Tests (Optional)

| ID | Test Case | Rationale |
|----|-----------|-----------|
| C1 | Very long ORF (>1000 aa) handled | Stress test |
| C2 | Case-insensitive sequence handling | Robustness |
| C3 | Multiple RBS sites upstream of same ORF | Edge case |

## Edge Cases & Boundaries

| Case | Input | Expected | Source |
|------|-------|----------|--------|
| Empty | "" | Empty result | Definition |
| Too short | "ATG" | Empty (no stop) | Definition |
| No start codon | "TAATAG" | Empty | Definition |
| Minimal ORF | ATG + 99aa + stop | Filtered if minOrfLength=100 | Parameter |
| SD at exact minDistance (4 bp) | padding+AGGAGG+4bp+ORF | Detected | Boundary, Wikipedia |
| SD at exact maxDistance (15 bp) | padding+AGGAGG+15bp+ORF | Detected | Boundary, Wikipedia |
| SD beyond maxDistance (16 bp) | padding+AGGAGG+16bp+ORF | Not detected | Boundary, Wikipedia |
| SD at optimal (5 nt) | padding+AGGAGG+5bp+ORF | Detected | Chen et al. 1994 |

## Test Pool Consolidation

| File | Status |
|------|--------|
| `GenomeAnnotator_Gene_Tests.cs` | Canonical: 32 tests (26 methods, 6 parametric cases) |
| `GenomeAnnotatorTests.cs` | Retained: GFF3/promoter/other tests |

## Coverage Classification

**Total: 32 tests in canonical file (was 24 — 1 duplicate removed, 5 missing implemented, 4 parametric cases added)**

### Summary

| Category | Count |
|----------|-------|
| ❌ Missing → Implemented | 5 (S2, S3-TTG, C1, C3, Edge-NoStop) |
| ⚠ Weak → Strengthened | 8 (M2, M5, M6, M10, S1, S3→parametric, S4, FrameAttribute) |
| 🔁 Duplicate → Removed | 1 (M7 RespectsDistanceConstraints) |
| ✅ Covered | All remaining |

### ❌ Missing (Implemented)

| ID | Test | Action |
|----|------|--------|
| S2 | Overlapping genes both reported | Implemented `PredictGenes_OverlappingGenes_BothReported` — two ORFs in different frames with verified overlap |
| S3 | TTG start codon not tested | Added `[TestCase("TTG")]` to parametric `AlternativeStartCodons_Recognized` |
| C1 | Very long ORF (>1000 aa) | Implemented `PredictGenes_VeryLongOrf_Handled` — 1500 aa ORF, verifies >3000 nt span |
| C3 | Multiple RBS upstream of same ORF | Implemented `FindRibosomeBindingSites_MultipleUpstreamOfSameOrf` — AGGAGG at spacing 13 + GAGG at spacing 5 |
| Edge | ATG without stop codon | Implemented `PredictGenes_StartCodonOnly_NoStop_ReturnsEmpty` |

### ⚠ Weak (Strengthened)

| ID | Test | Before | After |
|----|------|--------|-------|
| M2 | AssignsSequentialGeneIds | `genes[0]` only, `Does.Match` | All genes verified with exact `$"test_{i+1:D4}"` |
| M5 | FindsGenesOnBothStrands | Only forward strand asserted | Both `+` and `-` strands asserted via revcomp construction |
| M6 | DetectsConsensusAggagg | `Contains("AGGAGG") \|\| Contains("GGAGG")` | `s.sequence == "AGGAGG"` exact match |
| M10 | ProteinLengthAttribute | `LessThanOrEqualTo(nucleotideLength/3)` | Implementation bug fixed: `TrimEnd('*')`; test asserts exact `(End-Start)/3 - 1` |
| S1 | MultipleGenes_AllDetected | `GreaterThanOrEqualTo(2)` | Exact coordinate assertions: `[0,303)`, `[323,626)` |
| S3 | AlternativeStartCodons | Only GTG tested | Parametric: `[TestCase("GTG")][TestCase("TTG")]` |
| S4 | ShorterMotifs_Detected | `s.sequence.Contains(sdMotif)` | `s.sequence == sdMotif` exact match |
| — | FrameAttributeIsValid | `TryGetValue` (no assert on key existence) | `ContainsKey` assertion + exact range check |

### 🔁 Duplicate (Removed)

| Test | Reason |
|------|--------|
| `FindRibosomeBindingSites_RespectsDistanceConstraints` | Fully subsumed by `AtMinDistance_Detected` (4bp), `AtMaxDistance_Detected` (15bp), `BeyondMaxDistance_NotDetected` (16bp), `OptimalSpacing_5nt_ChenEtAl1994` (5nt) |

### Canonical File (`GenomeAnnotator_Gene_Tests.cs`) — 32 tests

| # | Test Method | Spec ID | Status |
|---|-------------|---------|--------|
| 1 | `PredictGenes_AllGenesHaveCdsType` | M1 | ✅ |
| 2 | `PredictGenes_AssignsSequentialGeneIds` | M2 | ✅ |
| 3 | `PredictGenes_IncludesStrandInformation` | M3 | ✅ |
| 4 | `PredictGenes_FiltersOrfsByMinLength` | M4 | ✅ |
| 5 | `PredictGenes_FindsGenesOnBothStrands` | M5 | ✅ |
| 6 | `PredictGenes_EmptySequence_ReturnsEmpty` | M8 | ✅ |
| 7 | `PredictGenes_NoValidOrfs_ReturnsEmpty` | M9 | ✅ |
| 8 | `PredictGenes_ProteinLengthAttributeIsAccurate` | M10 | ✅ |
| 9 | `PredictGenes_MultipleGenes_AllDetected` | S1 | ✅ |
| 10 | `PredictGenes_AlternativeStartCodons_Recognized` (2 cases) | S3 | ✅ |
| 11 | `FindRibosomeBindingSites_DetectsConsensusAggagg` | M6 | ✅ |
| 12 | `FindRibosomeBindingSites_TooClose_NotDetected` | M7/Edge | ✅ |
| 13 | `FindRibosomeBindingSites_OptimalSpacing_5nt_ChenEtAl1994` | M7/Boundary | ✅ |
| 14 | `FindRibosomeBindingSites_AtMinDistance_Detected` | M7/Boundary | ✅ |
| 15 | `FindRibosomeBindingSites_AtMaxDistance_Detected` | M7/Boundary | ✅ |
| 16 | `FindRibosomeBindingSites_BeyondMaxDistance_NotDetected` | M7/Boundary | ✅ |
| 17 | `FindRibosomeBindingSites_ShorterMotifs_Detected` (4 cases) | S4 | ✅ |
| 18 | `FindRibosomeBindingSites_ScoreReflectsMotifLength` | S5 | ✅ |
| 19 | `FindRibosomeBindingSites_NoOrfs_ReturnsEmpty` | Edge | ✅ |
| 20 | `PredictGenes_MixedCase_HandledCorrectly` | C2 | ✅ |
| 21 | `PredictGenes_NullSequence_HandledGracefully` | Edge | ✅ |
| 22 | `PredictGenes_DefaultPrefix_IsGene` | Edge | ✅ |
| 23 | `PredictGenes_CoordinatesAreValid` | Invariant | ✅ |
| 24 | `PredictGenes_FrameAttributeIsValid` | Invariant | ✅ |
| 25 | `PredictGenes_OverlappingGenes_BothReported` | S2 | ✅ |
| 26 | `PredictGenes_VeryLongOrf_Handled` | C1 | ✅ |
| 27 | `PredictGenes_StartCodonOnly_NoStop_ReturnsEmpty` | Edge | ✅ |
| 28 | `FindRibosomeBindingSites_MultipleUpstreamOfSameOrf` | C3 | ✅ |

### Classification Summary

- ✅ Covered: 32 tests (28 methods, 6 parametric cases)
- ❌ Missing: 0
- ⚠ Weak: 0
- 🔁 Duplicate: 0

## Deviations and Assumptions

None. All design parameters are grounded in external sources:

| Parameter | Value | Source |
|-----------|-------|--------|
| SD consensus | AGGAGG | Shine & Dalgarno (1975) |
| SD motifs | AGGAGG, GGAGG, AGGAG, GGAG, GAGG, AGGA | All contiguous ≥4-nt substrings of consensus; Prodigal sequence.c shine_dalgarno_exact / gene.c class "AGGA/GGAG/GAGG" |
| Functional range | 4-15 bp | Wikipedia: Shine-Dalgarno sequence |
| Optimal aligned spacing | 5 nt | Chen et al. (1994) |
| Start codons | ATG, GTG, TTG | Wikipedia: Gene prediction (prokaryotic) |
| Stop codons | TAA, TAG, TGA | Standard genetic code |
| Score normalization | motif.Length / 6.0 | Heuristic (declared; see Gene_Prediction.md §5.3) |
| One gene per stop | first start (longest ORF) | EMBOSS getorf -find 1 (numeric cross-check); Prodigal |
| Prokaryotic model | No introns, ORF-based | Wikipedia: Gene prediction |

## Validation Checklist

- [x] All MUST tests have evidence source
- [x] Invariants specified and tested
- [x] Edge cases documented and tested (empty, null, no-stop, no-start, too-close, beyond-max)
- [x] Boundary conditions tested (minDistance=4, maxDistance=15, beyond=16, optimal=5)
- [x] Both strands tested
- [x] All alternative start codons tested (GTG, TTG)
- [x] Sequential ID numbering verified across multiple genes
- [x] No assumptions — all design decisions backed by external sources
- [x] No duplicates — each test serves a distinct purpose
- [x] Coverage classification complete: 0 missing, 0 weak, 0 duplicate
- [x] Tests passing (32/32)

## Review 2026-09 (B11) additions

| Test | Reference |
|------|-----------|
| `FindRibosomeBindingSites_GgagSubMotif_Detected` | Prodigal shine_dalgarno_exact: GGAG cur_ctr = 9 (= GAGG), class AGGA/GGAG/GAGG |
| `FindRibosomeBindingSites_SingleSite_ReportedOnceAsMaximalMotif` | Prodigal returns one maximal sub-motif per window; AGGAGG@10 → single hit (10, AGGAGG, 1.0) |
| `FindRibosomeBindingSites_TooCloseConsensus_OnlyInRangeSubMotifReported` | AGGAGG at spacer 2 → only AGGA [10,14), spacer 4 (Prodigal limit=4 → AGGA, rdis 4, bin 11) |
| `FindRibosomeBindingSites_NestedStartsSharingSd_ReportedOnce` | One site, two nested starts → one hit |
| `PredictGenes_NestedStartsSharingStop_OneGenePerStop_MatchesGetorf` | EMBOSS getorf 6.6.0 -find 1 -table 11: [3-176], [324-190] REVERSE → +[2,179), −[186,324) |
| `AnnotationFuzzTests.PredictGenes_OverlappingGenes_AllReportedWithValidCoordinates` | updated: nested ORF sharing a stop → single gene [0, len) |

