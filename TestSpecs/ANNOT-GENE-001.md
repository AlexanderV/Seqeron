# Test Specification: ANNOT-GENE-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | ANNOT-GENE-001 |
| **Area** | Annotation |
| **Algorithm** | Gene Prediction |
| **Canonical Methods** | `GenomeAnnotator.PredictGenes`, `GenomeAnnotator.FindRibosomeBindingSites` |
| **Complexity** | O(n) |

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Gene prediction | Encyclopedia | ORF-based gene finding in prokaryotes, signal detection |
| Wikipedia: Shine-Dalgarno sequence | Encyclopedia | SD consensus AGGAGG, distance 4-15bp upstream of start |
| Wikipedia: Ribosome-binding site | Encyclopedia | RBS role in translation initiation |
| Shine & Dalgarno (1975) | Primary | Original SD sequence identification |
| Chen et al. (1994) | Primary | Optimal SD-to-start spacing: 5-9 bp |
| Laursen et al. (2005) | Review | Bacterial translation initiation mechanisms |

## Method Specifications

### 1. PredictGenes(dna, minOrfLength, prefix)

**Description:** Predicts genes using ORF-based approach on both strands.

**Invariants:**
- All returned genes have Type = "CDS"
- All genes have strand '+' or '-'
- All genes have Start < End
- Gene IDs follow pattern "{prefix}_{number:D4}"
- Protein length in attributes matches (End-Start)/3 - 1

### 2. FindRibosomeBindingSites(dna, upstreamWindow, minDistance, maxDistance)

**Description:** Finds Shine-Dalgarno sequences upstream of ORFs.

**Invariants:**
- Position is within upstream window of a valid ORF
- Sequence matches one of: AGGAGG, GGAGG, AGGAG, GAGG, AGGA
- Distance to start codon is within [minDistance, maxDistance]
- Score is normalized (motif.Length / 6.0)

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
| M7 | FindRibosomeBindingSites validates distance constraints | Optimal spacing 4-15bp | Chen et al. 1994 |
| M8 | Empty sequence returns empty result | Edge case | Implementation |
| M9 | Sequence without ORFs returns empty | No start/stop = no gene | Definition |
| M10 | PredictGenes protein length in attributes is accurate | Data integrity | Implementation |

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
| Exactly minimum | ATG + 100aa + stop | Included | Boundary |

## Test Pool Consolidation

### Current State
- Existing tests in `GenomeAnnotatorTests.cs` (lines 69-97): 4 basic tests
- Tests are minimal smoke tests, not evidence-based

### Target State
- Create `GenomeAnnotator_Gene_Tests.cs` for ANNOT-GENE-001 canonical tests
- Move/refactor gene prediction tests from GenomeAnnotatorTests.cs
- Keep only GFF3/promoter/other tests in GenomeAnnotatorTests.cs

### Files to Modify
| File | Action |
|------|--------|
| `GenomeAnnotator_Gene_Tests.cs` | CREATE: canonical evidence-based tests |
| `GenomeAnnotatorTests.cs` | MODIFY: remove PredictGenes/RBS tests (move to new file) |

## Assumptions

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Prokaryotic model is appropriate | Implementation targets bacteria |
| A2 | SD motif set is sufficient | Common variants covered |
| A3 | Default distance 4-15bp is reasonable | Matches literature |

## Open Questions

None - implementation matches documented biological behavior.
