# Test Specification: CHROM-KARYO-001

**Test Unit ID:** CHROM-KARYO-001  
**Area:** Chromosome Analysis  
**Title:** Karyotype Analysis  
**Date:** 2026-02-01  
**Status:** Complete  

---

## Scope

### Canonical Methods (Deep Testing Required)
| Method | Class | Type |
|--------|-------|------|
| `AnalyzeKaryotype(chromosomes, ploidy)` | ChromosomeAnalyzer | Canonical |
| `DetectPloidy(depths, expected)` | ChromosomeAnalyzer | Canonical |

### Wrappers/Delegates (Smoke Tests Only)
None identified.

---

## Test Classification

### Must Tests (Evidence-Backed)

| ID | Test | Rationale | Source |
|----|------|-----------|--------|
| M1 | AnalyzeKaryotype with normal diploid set returns correct counts | Standard human karyotype (46 chromosomes) | Wikipedia Karyotype |
| M2 | AnalyzeKaryotype with trisomy detects aneuploidy | Trisomy = 3 copies of chromosome | Wikipedia Aneuploidy |
| M3 | AnalyzeKaryotype with monosomy detects aneuploidy | Monosomy = 1 copy of chromosome | Wikipedia Aneuploidy |
| M4 | AnalyzeKaryotype empty input returns empty karyotype | Graceful degradation | ASSUMPTION |
| M5 | DetectPloidy with diploid depth returns ploidy 2 | Ratio ≈ 1.0 → diploid | Wikipedia Ploidy |
| M6 | DetectPloidy with tetraploid depth returns ploidy 4 | Ratio ≈ 2.0 → tetraploid | Wikipedia Ploidy |
| M7 | DetectPloidy empty input returns default | Graceful degradation | ASSUMPTION |

### Should Tests

| ID | Test | Rationale |
|----|------|-----------|
| S1 | AnalyzeKaryotype correctly separates sex chromosomes | Karyotype distinguishes autosomes from allosomes |
| S2 | AnalyzeKaryotype calculates total genome size correctly | Sum of all chromosome lengths |
| S3 | AnalyzeKaryotype calculates mean chromosome length | TotalSize / TotalChromosomes |
| S4 | DetectPloidy with haploid depth returns ploidy 1 | Ratio ≈ 0.5 → haploid |
| S5 | DetectPloidy confidence decreases with noisy data | Confidence reflects certainty |
| S6 | AnalyzeKaryotype with custom ploidy level works | Support for polyploid organisms |

### Could Tests

| ID | Test | Rationale |
|----|------|-----------|
| C1 | DetectPloidy handles extreme ploidy values | Clamp to [1, 8] |
| C2 | AnalyzeKaryotype with multiple aneuploidies | Multiple abnormalities |

---

## Invariants

### AnalyzeKaryotype Invariants
1. `TotalChromosomes == AutosomeCount + SexChromosomes.Count`
2. `TotalGenomeSize == Σ(chromosome.Length)`
3. `MeanChromosomeLength == TotalGenomeSize / TotalChromosomes` (when TotalChromosomes > 0)
4. `HasAneuploidy == Abnormalities.Count > 0`

### DetectPloidy Invariants
1. `PloidyLevel ∈ [1, 8]`
2. `Confidence ∈ [0, 1]`
3. Uniform depth data → high confidence
4. Empty input → (2, 0) default

---

## Edge Cases

| Category | Input | Expected | Covered |
|----------|-------|----------|---------|
| Empty | No chromosomes | Empty karyotype | ✓ |
| Empty | No depth values | (2, 0) | ✓ |
| Boundary | Single chromosome | Monosomy detected | ✓ |
| Boundary | Ploidy ratio at 0.5 | Haploid (1) | ✓ |
| Boundary | Ploidy ratio at 1.0 | Diploid (2) | ✓ |
| Boundary | Ploidy ratio at 2.0 | Tetraploid (4) | ✓ |
| Error | Very high ploidy | Clamped to 8 | ✓ |

---

## Test File Structure

**Canonical Test File:** `ChromosomeAnalyzer_Karyotype_Tests.cs`

### Test Organization
```
ChromosomeAnalyzer_Karyotype_Tests
├── AnalyzeKaryotype Tests
│   ├── Normal diploid karyotype
│   ├── Trisomy detection
│   ├── Monosomy detection
│   ├── Multiple autosomes
│   ├── Sex chromosome handling
│   ├── Custom ploidy level
│   ├── Empty input
│   └── Invariant checks
├── DetectPloidy Tests
│   ├── Haploid detection
│   ├── Diploid detection
│   ├── Tetraploid detection
│   ├── High ploidy clamping
│   ├── Confidence calculation
│   ├── Empty input
│   └── Noisy data handling
```

---

## Audit Notes

### Previous Test Coverage (ChromosomeAnalyzerTests.cs)
- Tests existed in mixed file with telomere, centromere, and synteny tests
- 7 karyotype-related tests identified
- Coverage: Adequate for basic scenarios, missing invariant tests

### Consolidation Actions
1. Extract karyotype tests to dedicated `ChromosomeAnalyzer_Karyotype_Tests.cs`
2. Add missing Should tests (S1-S6)
3. Add invariant verification with Assert.Multiple
4. Keep telomere/centromere tests in their respective files per prior Test Units

---

## Open Questions

None.

---

## Decisions

1. **D1:** Maintain karyotype tests separate from other chromosome analysis tests for clarity
2. **D2:** Use Assert.Multiple for invariant checking to catch all violations
3. **D3:** Default ploidy for empty input is (2, 0) as per implementation
