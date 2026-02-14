# Test Specification: SEQ-GCSKEW-001

**Test Unit ID:** SEQ-GCSKEW-001
**Algorithm:** GC Skew Analysis
**Area:** Sequence Composition
**Status:** ☑ Complete
**Last Updated:** 2026-02-14
**Owner:** QA Architect

---

## 1. Scope

This Test Unit covers all GC/AT skew and related analysis methods in `GcSkewCalculator`:

| Method | Type | Description |
|--------|------|-------------|
| `CalculateGcSkew(DnaSequence)` | Canonical | Calculate overall GC skew |
| `CalculateGcSkew(string)` | Overload | String input variant |
| `CalculateWindowedGcSkew(...)` | Windowed | Sliding window analysis |
| `CalculateCumulativeGcSkew(...)` | Cumulative | Cumulative GC skew for origin detection |
| `CalculateAtSkew(DnaSequence)` | Canonical | Calculate overall AT skew |
| `CalculateAtSkew(string)` | Overload | String input variant |
| `PredictReplicationOrigin(...)` | Derived | Predict replication origin/terminus |
| `AnalyzeGcContent(...)` | Composite | Comprehensive GC analysis with windowed metrics |

---

## 2. Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia - GC Skew](https://en.wikipedia.org/wiki/GC_skew) | Authoritative | Formula: GC skew = (G − C)/(G + C); AT skew = (A − T)/(A + T); Range: [−1, +1]; Sign switch at ori/ter |
| Lobry, J.R. (1996) Mol. Biol. Evol. 13:660-665 | Primary | Original GC skew observations in bacterial genomes; deviation from [C]=[G] as (C−G)/(C+G) |
| Grigoriev, A. (1998) Nucleic Acids Res. 26:2286-2290 | Primary | Cumulative GC skew method; minimum = origin, maximum = terminus (with modern (G−C)/(G+C) definition) |
| Chargaff's Rule (1950) | Foundational | G≈C and A≈T in double-stranded DNA (parity rule 1) |
| [Biopython Bio.SeqUtils.GC_skew](https://biopython.org/docs/latest/api/Bio.SeqUtils.html) | Cross-verification | `GC_skew(seq, window)`: formula (G−C)/(G+C), returns 0.0 on ZeroDivisionError, case-insensitive |

### Key Evidence Points

1. **Formula Definition** (Wikipedia, Lobry 1996):
   - Modern definition: GC skew = (G − C) / (G + C)
   - AT skew = (A − T) / (A + T)
   - Where G, C, A, T are counts of the respective bases in a defined window

2. **Value Range** (Wikipedia):
   - Range: −1 ≤ skew ≤ +1
   - −1: C only, no G (G = 0)
   - +1: G only, no C (C = 0)
   - 0: Equal G and C, OR no G and C present

3. **Biological Significance** (Lobry 1996, Grigoriev 1998):
   - Leading strand: typically positive GC skew (G > C)
   - Lagging strand: typically negative GC skew (C > G)
   - GC skew sign changes at replication origin and terminus

4. **Cumulative GC Skew** (Grigoriev 1998, Wikipedia CGC skew section):
   - Sum of window skews from arbitrary start
   - With (G−C)/(G+C) definition: global minimum = origin of replication (oriC), global maximum = terminus (ter)

5. **Zero-division Handling** (Biopython source):
   - When G + C = 0 (e.g., all A/T or empty): returns 0.0
   - Biopython: `except ZeroDivisionError: skew = 0.0`

6. **Case Insensitivity** (Biopython source):
   - Biopython counts both cases: `s.count("G") + s.count("g")`
   - Standard bioinformatics practice for sequence analysis tools

### Biopython Cross-Verification Table

| Input | Biopython Call | Expected | Formula |
|-------|---------------|----------|---------|
| GGGGC | `calc_gc_skew("GGGGC")` | 0.6 | (4−1)/(4+1) |
| GCCC | `calc_gc_skew("GCCC")` | −0.5 | (1−3)/(1+3) |
| GGGGCCCC, w=4 | `GC_skew("GGGGCCCC", 4)` | [1.0, −1.0] | GGGG→1.0, CCCC→−1.0 |
| ATGCATGC, w=4 | `GC_skew("ATGCATGC", 4)` | [0.0, 0.0] | Each window: (1−1)/(1+1) |
| AAAAAAAA, w=4 | `GC_skew("AAAAAAAA", 4)` | [0.0, 0.0] | ZeroDivisionError → 0.0 |
| "A"×50 | `GC_skew("A"*50)[0]` | 0 | Biopython test_GC_skew: assertEqual(0) |

---

## 3. Invariants

| ID | Invariant | Source |
|----|-----------|--------|
| INV-1 | −1 ≤ GC skew ≤ +1 | Wikipedia: bounded by definition of (G−C)/(G+C) |
| INV-2 | All-G sequence → skew = +1 | Formula: (G−0)/(G+0) = 1 |
| INV-3 | All-C sequence → skew = −1 | Formula: (0−C)/(0+C) = −1 |
| INV-4 | Equal G,C → skew = 0 | Formula: (n−n)/(n+n) = 0 |
| INV-5 | No G,C (all A,T) → skew = 0 | Biopython: ZeroDivisionError → 0.0 |
| INV-6 | G↔C swap negates skew | Formula: swapping G,C in (G−C)/(G+C) yields (C−G)/(G+C) = −skew |
| INV-7 | Cumulative min = origin, max = terminus | Grigoriev 1998 (with (G−C)/(G+C) definition) |

---

## 4. Test Classification

### 4.1 MUST Tests (Evidence-Based)

| Test ID | Scenario | Expected | Source |
|---------|----------|----------|--------|
| M-01 | Formula: GGGGC (4G,1C) | (4−1)/(4+1) = 0.6 | Formula |
| M-02 | Formula: CCCCC (0G,5C) | (0−5)/(0+5) = −1.0 | Formula |
| M-03 | Formula: GGGGG (5G,0C) | (5−0)/(5+0) = +1.0 | Formula |
| M-04 | Formula: GCGC (2G,2C) | (2−2)/(2+2) = 0.0 | Formula |
| M-05 | No G/C: AAATTT | 0 (G+C=0 → zero-division → 0) | Biopython: ZeroDivisionError → 0.0 |
| M-06 | Empty sequence | 0 (G=0,C=0 → zero-division → 0) | Biopython: ZeroDivisionError → 0.0 |
| M-07 | Invariant: all results in [−1, +1] | Range check | Wikipedia |
| M-08 | Windowed: positions at window centers | Position = WindowStart + WindowSize/2 | Grigoriev 1998 |
| M-09 | Cumulative: accumulates window skews | Values sum correctly | Grigoriev 1998 |
| M-10 | Cumulative: GGGG→CCCC pattern | +1, 0, +1, 0 oscillation | Grigoriev 1998 |
| M-11 | Origin detection: minimum at expected position | Min cumulative = origin | Grigoriev 1998 |
| M-12 | Terminus detection: maximum at expected position | Max cumulative = terminus | Grigoriev 1998 |
| M-13 | Null sequence → ArgumentNullException | Exception | Guard clause |
| M-14 | Window size ≤ 0 → ArgumentOutOfRangeException | Exception | Guard clause |
| M-15 | Step size ≤ 0 → ArgumentOutOfRangeException | Exception | Guard clause |
| M-16 | Case insensitivity: lowercase handled | Same result as uppercase | Biopython: counts both cases |
| M-17 | Sequence shorter than window → empty result | No windows fit | Logic |
| M-18 | Biopython cross-verification: GGGGCCCC w=4 | [1.0, −1.0] | Biopython GC_skew |
| M-19 | Biopython cross-verification: ATGCATGC w=4 | [0.0, 0.0] | Biopython GC_skew |
| M-20 | Biopython cross-verification: AAAAAAAA w=4 | [0.0, 0.0] | Biopython: ZeroDivisionError |
| M-21 | Biopython cross-verification: GCCC single | −0.5 | Biopython calc_gc_skew |

### 4.2 SHOULD Tests

| Test ID | Scenario | Expected | Source |
|---------|----------|----------|--------|
| S-01 | AT skew formula: AAAAT → (4−1)/(4+1) = 0.6 | 0.6 | Wikipedia: AT skew formula |
| S-02 | Windowed overlapping windows produce more points | More points | Grigoriev 1998: sliding window |
| S-03 | GC analysis result contains all metrics | All fields populated | Composite method |

### 4.3 COULD Tests

| Test ID | Scenario | Expected |
|---------|----------|----------|
| C-01 | Large sequence performance | Completes in reasonable time |
| C-02 | Windowed with various step sizes | Correct point counts |

---

## 5. Edge Cases

| Category | Case | Expected Behavior | Source |
|----------|------|-------------------|--------|
| Empty | Empty string | Return 0 (G=0,C=0, G+C=0 → zero-division → 0) | Biopython: ZeroDivisionError → 0.0 |
| Empty | Windowed on empty | Return empty collection (no windows fit) | Biopython: `GC_skew("")` → `[]` |
| Null | Null DnaSequence | Throw ArgumentNullException | Guard clause |
| Boundary | All G | +1.0 | Formula: (G−0)/(G+0) |
| Boundary | All C | −1.0 | Formula: (0−C)/(0+C) |
| Boundary | No G or C | 0.0 | Biopython: ZeroDivisionError → 0.0 |
| Window | Window > sequence length | Return empty | Logic: no windows fit |
| Window | Window = sequence length | Return single point | Logic: exactly one window |

---

## 6. Coverage Classification

Applied systematic coverage classification (2026-02-14):

| Action | Tests | Details |
|--------|-------|---------|
| ✅ Covered | 30 | No changes needed |
| ⚠ Strengthened | 6 | `LowercaseInput` (exact 0.6), `OverlappingWindows` (exact 3 windows + values), `FindsMinimum` (exact pos=145), `FindsMaximum` (exact pos=145), `ReturnsValidPositions→ExactPositionsAndSkew` (deterministic), `ReturnsOverallMetrics` (exact 66.667%) |
| 🔁 Removed/Merged | 8 | `RangeIsMinusOneToOne` (FsCheck covers), `AllSequences_RespectRangeInvariant` (FsCheck covers), `ReturnsMultiplePoints` (in CorrectSkewValues), `BiopythonCrossVerification_GGGGCCCC` (= CorrectSkewValues), `IncludesWindowSkew` (in AccumulatesCorrectly), Properties: `HomopolymerC/G/EqualGandC` (in main) |
| ❌ Added | 4 | `CumulativeGcSkew_ZeroWindowSize`, `AtSkew_NullSequence`, `AtSkew_EmptySequence`, `PredictReplicationOrigin_TooShortSequence` |

---

## 7. Test Summary

| File | Tests | Coverage |
|------|-------|----------|
| `GcSkewCalculatorTests.cs` | 41 methods | Formula, windowed, cumulative, origin/terminus, AT skew, edge cases, cross-verification |
| `GcSkewProperties.cs` | 3 FsCheck properties | Range invariant (×100), AT range (×100), complement negation (×100) |
| **Total** | **44 test methods** (+ 300 FsCheck runs) | All MUST, SHOULD, COULD scenarios covered |

---

## 8. Deviations

None — implementation matches all evidence sources.

---

## 9. Assumptions

None — all behaviors sourced to external references.

---

## 10. Open Questions

None — all behaviors are well-defined by sources.

---

## 11. Validation Checklist

- [x] TestSpec created
- [x] Evidence sources documented (Wikipedia, Lobry 1996, Grigoriev 1998, Biopython)
- [x] Must/Should/Could tests defined
- [x] Edge cases enumerated and sourced
- [x] Invariants documented
- [x] No deviations
- [x] No assumptions — all behaviors sourced
- [x] Biopython cross-verification table added
- [x] Coverage classification applied
- [x] All 44 tests pass
