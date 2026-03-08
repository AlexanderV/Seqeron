# Test Specification: CHROM-KARYO-001

**Test Unit ID:** CHROM-KARYO-001  
**Area:** Chromosome Analysis  
**Title:** Karyotype Analysis  
**Date:** 2026-03-08  
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
| M4 | AnalyzeKaryotype empty input returns empty karyotype | Graceful degradation | Design Decision DD1 |
| M5 | DetectPloidy with diploid depth returns ploidy 2 | Ratio = 1.0 → ploidy 2, confidence 1.0 | Wikipedia Ploidy |
| M6 | DetectPloidy with tetraploid depth returns ploidy 4 | Ratio = 2.0 → ploidy 4, confidence 1.0 | Wikipedia Ploidy |
| M7 | DetectPloidy empty input returns default | Default diploid, zero confidence | Design Decision DD2 |
| M8 | Tetrasomy (4 copies) uses correct ISCN nomenclature | Tetrasomy = 4 copies | Wikipedia Aneuploidy |
| M9 | Pentasomy (5 copies) uses correct ISCN nomenclature | Pentasomy = 5 copies | Wikipedia Aneuploidy |
| M10 | Tetraploid context: 3 copies = Trisomy (absolute count) | Terminology is absolute, not relative | Wikipedia Aneuploidy |
| M11 | Disomy (2 copies) uses correct ISCN nomenclature in non-diploid context | 2 copies = Disomy per standard nomenclature | Wikipedia Aneuploidy |

### Should Tests

| ID | Test | Rationale |
|----|------|-----------|
| S1 | AnalyzeKaryotype correctly separates sex chromosomes | Karyotype distinguishes autosomes from allosomes |
| S2 | AnalyzeKaryotype calculates total genome size correctly | Sum of all chromosome lengths |
| S3 | AnalyzeKaryotype calculates mean chromosome length | TotalSize / TotalChromosomes |
| S4 | DetectPloidy with haploid depth returns ploidy 1 | Ratio = 0.5 → ploidy 1, confidence 1.0 |
| S5 | DetectPloidy between-ploidy ratio reduces confidence | Confidence = 1.0 − |ratio×2 − ploidy| × 2 |
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
│   ├── Normal diploid karyotype (M1)
│   ├── Trisomy detection (M2)
│   ├── Monosomy detection (M3)
│   ├── Empty input (M4)
│   ├── Sex chromosome handling (S1)
│   ├── Total genome size (S2)
│   ├── Mean chromosome length (S3)
│   ├── Custom ploidy level (S6)
│   ├── Multiple aneuploidies (C2)
│   ├── Tetrasomy nomenclature (M8)
│   ├── Pentasomy nomenclature (M9)
│   ├── Tetraploid absolute terminology (M10)
│   ├── Disomy in non-diploid context (M11)
│   ├── Invariant: Total = Autosomes + Sex
│   └── Invariant: Aneuploidy ↔ Abnormalities
├── DetectPloidy Tests
│   ├── Diploid detection (M5)
│   ├── Tetraploid detection (M6)
│   ├── Haploid detection (S4)
│   ├── Triploid detection
│   ├── Empty input (M7)
│   ├── High ploidy clamping (C1)
│   ├── Low ploidy clamping (C1)
│   ├── Between-ploidy confidence (S5)
│   ├── Single value input
│   ├── Custom expected depth
│   ├── Invariant: Ploidy in [1, 8]
│   └── Invariant: Confidence in [0, 1]
```

---

## Audit Notes

None pending.

---

## Open Questions

None.

---

## Decisions

1. **D1:** Maintain karyotype tests separate from other chromosome analysis tests for clarity
2. **D2:** Use Assert.Multiple for invariant checking to catch all violations
3. **D3:** Default ploidy for empty input is (2, 0) — diploid is most common; zero confidence signals no data
4. **D4:** Aneuploidy terminology uses absolute copy count per ISCN / Wikipedia Aneuploidy standard
