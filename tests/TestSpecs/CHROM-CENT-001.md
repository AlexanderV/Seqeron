# Test Specification: CHROM-CENT-001

## Test Unit Information
- **ID:** CHROM-CENT-001
- **Area:** Chromosome
- **Canonical Method:** `ChromosomeAnalyzer.AnalyzeCentromere(chromosomeName, sequence, windowSize, minAlphaSatelliteContent)`
- **Complexity:** O(n)

## Methods Under Test

| Method | Class | Type | Deep Test |
|--------|-------|------|-----------|
| `AnalyzeCentromere(chrName, seq, windowSize)` | ChromosomeAnalyzer | Canonical | Yes |
| `DetectAlphaSatellite(seq)` | ChromosomeAnalyzer | Canonical (alpha-satellite-specific) | Yes |
| `FindCenpBBoxes(seq)` | ChromosomeAnalyzer | Canonical (CENP-B box scan) | Yes |

## Test Classification

### Must Tests (Evidence-backed)

| ID | Test Case | Rationale | Source | Status |
|----|-----------|-----------|--------|--------|
| M1 | Empty sequence returns Unknown type with null boundaries | Standard edge case | Wikipedia | ✅ |
| M1b | Null sequence returns Unknown type with null boundaries | Standard edge case | Wikipedia | ✅ |
| M2 | Sequence shorter than window size returns Unknown | Cannot perform analysis without sufficient data | Implementation logic | ✅ |
| M3 | Chromosome name is preserved in result | Basic correctness | API contract | ✅ |
| M4 | Repetitive region is detected as centromeric (Start/End non-null, type not Unknown) | Centromeres are characterized by repetitive DNA | Wikipedia | ✅ |
| M4b | Non-repetitive sequence returns Unknown with null boundaries | Inverse of M4 | Wikipedia | ✅ |
| M5 | CentromereResult invariants: Start <= End when found | Structural invariant | Logic | ✅ |
| M6 | Length equals (End - Start) when centromere found | Structural invariant | Logic | ✅ |
| M7 | Type classification based on arm ratio per Levan (1964) | Metacentric ≤1.7, Submetacentric (1.7,3.0], Subtelocentric (3.0,7.0), Acrocentric ≥7.0 | Levan et al. (1964) | ✅ |
| M8 | IsAcrocentric flag matches type (true for Acrocentric, false otherwise) | Consistency check | API contract | ✅ |

### Should Tests (Recommended)

| ID | Test Case | Rationale | Status |
|----|-----------|-----------|--------|
| S1 | Different window sizes all detect large repetitive region | Parameter behavior | ✅ |
| S2 | Low threshold detects centromere; high threshold reduces sensitivity | Parameter behavior | ✅ |
| S3 | Case insensitivity: uppercase and lowercase produce identical results | Robustness | ✅ |

### Could Tests (Optional)

| ID | Test Case | Rationale | Status |
|----|-----------|-----------|--------|
| C1 | Performance with large sequences | Non-functional | Not implemented |

## Edge Cases

| Case | Input | Expected Output | Status |
|------|-------|-----------------|--------|
| Empty sequence | `""` | Type = "Unknown", Start = null, End = null | ✅ |
| Null sequence | `null` | Type = "Unknown", Start = null, End = null | ✅ |
| Very short sequence | `"ATCG"` (shorter than window) | Type = "Unknown" | ✅ |
| All same base | `"AAAA..."` (300kb) | Detected (maximally repetitive) | ✅ |
| No repetitive regions | Random non-repetitive | Type = "Unknown", Start/End = null | ✅ |

## Test Invariants

1. **Result chromosome name equals input name** ✅
2. **When Start/End are both null, Length = 0** ✅ (verified in M1, M4b)
3. **When Start/End are both non-null, Length = End - Start** ✅
4. **Type is one of: Metacentric, Submetacentric, Subtelocentric, Acrocentric, Telocentric, Unknown** ✅
5. **IsAcrocentric = true iff Type = "Acrocentric"** ✅
6. **AlphaSatelliteContent >= 0** ✅

## Classification Basis

Per Levan A, Fredga K, Sandberg AA (1964) "Nomenclature for centromeric position on chromosomes", Hereditas 52(2):201-220:

| Arms length ratio (q/p) | Classification |
|--------------------------|----------------|
| 1.0 – 1.7               | Metacentric    |
| (1.7) – 3.0             | Submetacentric |
| (3.0) – (7.0)           | Subtelocentric |
| ≥ 7.0                   | Acrocentric    |
| ∞ (p = 0)               | Telocentric    |

Classification uses arm ratio (q/p) where p = short arm, q = long arm, computed from the centromere midpoint position.

**Boundary values per Levan table:** 1.7 → Metacentric, 3.0 → Submetacentric, 7.0 → Acrocentric.

**Implementation note:** Telocentric (p = 0) is handled in code but unreachable through `AnalyzeCentromere` since the sliding window detection always produces a non-zero centromere midpoint.

## Constants Tests

| Test | Status |
|------|--------|
| AlphaSatelliteConsensus: non-empty, length > 50, valid DNA bases only | ✅ |
| AlphaSatelliteMonomerLength == 171 (Willard 1985 / PMC6121732) | ✅ |
| CenpBBoxConsensus == "YTTCGTTGGAARCGGGA", length 17 (Masumoto 1989) | ✅ |

## Alpha-Satellite-Specific Detection (added 2026-06-24)

Tests in `ChromosomeAnalyzer_AlphaSatellite_Tests.cs`. These verify the **alpha-satellite-specific**
methods, distinct from the generic `AnalyzeCentromere` repeat heuristic. All expected values are
hand-derived from the synthetic fixtures and the sourced parameters in the Evidence (171-bp tandem
period; AT > 0.50; CENP-B box `YTTCGTTGGAARCGGGA`, Y=C/T, R=A/G).

### Must Tests

| ID | Test Case | Expected (exact) | Source |
|----|-----------|------------------|--------|
| M-ALPHA-1 | Perfect tandem 171-bp AT-rich array detected | IsAlphaSatellite=true; PeriodicityScore=1.0; BestPeriod=171; AtContent=100/171 | 171 bp monomer (PMC6121732) |
| M-ALPHA-2 | Random 4-letter sequence not detected | IsAlphaSatellite=false; PeriodicityScore<0.50 | — |
| M-ALPHA-3 | AT-rich (~0.70) but non-repetitive not detected | AtContent>0.50; PeriodicityScore<0.50; IsAlphaSatellite=false | AT-richness alone insufficient |
| M-ALPHA-4 | GC-rich 16-bp tandem repeat not detected | AtContent=0.0; IsAlphaSatellite=false | not AT-rich, not 171-bp period |
| M-ALPHA-5 | Tandem array, 1 CENP-B box per 171-bp monomer × 10 copies | IsAlphaSatellite=true; BestPeriod=171; CenpBBoxCount=10 | Masumoto 1989 |
| C-ALPHA-1 | Canonical box `CTTCGTTGGAAACGGGA` matches at index 0 | 1 hit at position 0 | Masumoto 1989 |
| C-ALPHA-2 | All four Y/R ambiguity resolutions match | 1 hit each | IUPAC Y=C/T, R=A/G |
| C-ALPHA-3 | Mutating a fixed consensus base prevents a match | 0 hits | exact-match requirement |
| C-ALPHA-4 | A (≠ C/T) in the leading Y position prevents a match | 0 hits | Y ambiguity scope |
| C-ALPHA-5 | Box embedded after 50-bp flank reported at offset 50 | 1 hit at position 50 | 0-based offset contract |

### Edge Cases (alpha-satellite methods)

| Case | Input | Expected | Status |
|------|-------|----------|--------|
| Empty / null sequence | `""` / `null` | DetectAlphaSatellite → no-signal (false, 0,0,0,0); FindCenpBBoxes → empty | ✅ |
| Too short to measure period | 100 bp | DetectAlphaSatellite → BestPeriod=0, PeriodicityScore=0, false | ✅ |
| 16-bp input to FindCenpBBoxes | 16 bp | empty (shorter than 17-bp box) | ✅ |
| Mixed-case input | lower vs upper | identical result | ✅ |

### Residual / out of scope
Higher-order repeat (HOR) structure and suprachromosomal family classification are not modelled;
`DetectAlphaSatellite` detects the monomer-level (171-bp tandem + AT + CENP-B) signal only.
