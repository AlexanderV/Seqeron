# Test Specification: SEQ-GCSKEW-001

**Test Unit ID:** SEQ-GCSKEW-001
**Algorithm:** GC Skew Analysis
**Area:** Sequence Composition
**Status:** ☑ Complete
**Last Updated:** 2026-01-22
**Owner:** QA Architect

---

## 1. Scope

This Test Unit covers all GC skew calculation methods in `GcSkewCalculator`:

| Method | Type | Description |
|--------|------|-------------|
| `CalculateGcSkew(DnaSequence)` | Canonical | Calculate overall GC skew |
| `CalculateGcSkew(string)` | Overload | String input variant |
| `CalculateWindowedGcSkew(...)` | Windowed | Sliding window analysis |
| `CalculateCumulativeGcSkew(...)` | Cumulative | Cumulative GC skew for origin detection |
| `FindOriginOfReplication(...)` | Derived | Predict replication origin/terminus |

---

## 2. Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia - GC Skew](https://en.wikipedia.org/wiki/GC_skew) | Authoritative | Formula: GC skew = (G - C)/(G + C); Range: [-1, +1]; Sign switch at ori/ter |
| Lobry, J.R. (1996) Mol. Biol. Evol. 13:660-665 | Primary | Original GC skew observations in bacterial genomes |
| Grigoriev, A. (1998) Nucleic Acids Res. 26:2286-2290 | Primary | Cumulative GC skew method; minimum = origin, maximum = terminus |
| Chargaff's Rule (1950) | Foundational | G≈C and A≈T in double-stranded DNA (parity rule 1) |

### Key Evidence Points

1. **Formula Definition** (Wikipedia, Lobry 1996):
   - GC skew = (G - C) / (G + C)
   - Where G and C are counts of guanine and cytosine bases

2. **Value Range** (Wikipedia):
   - Range: -1 ≤ skew ≤ +1
   - -1: C only, no G (G = 0)
   - +1: G only, no C (C = 0)
   - 0: Equal G and C, OR no G and C present

3. **Biological Significance** (Lobry 1996, Grigoriev 1998):
   - Leading strand: typically positive GC skew (G > C)
   - Lagging strand: typically negative GC skew (C > G)
   - GC skew sign changes at replication origin and terminus

4. **Cumulative GC Skew** (Grigoriev 1998):
   - Sum of window skews from arbitrary start
   - Global minimum corresponds to origin of replication (oriC)
   - Global maximum corresponds to terminus of replication (ter)

5. **Edge Cases** (Derived from formula):
   - Empty sequence: undefined mathematically; implementation returns 0
   - No G or C bases: returns 0 (division by zero protection)

---

## 3. Invariants

| ID | Invariant | Source |
|----|-----------|--------|
| INV-1 | -1 ≤ GC skew ≤ +1 | Wikipedia |
| INV-2 | All-G sequence → skew = +1 | Formula |
| INV-3 | All-C sequence → skew = -1 | Formula |
| INV-4 | Equal G,C → skew = 0 | Formula |
| INV-5 | No G,C (all A,T) → skew = 0 | Division protection |
| INV-6 | Cumulative min position < terminus position in typical bacterial genome | Grigoriev 1998 |

---

## 4. Test Classification

### 4.1 MUST Tests (Evidence-Based)

| Test ID | Scenario | Expected | Source |
|---------|----------|----------|--------|
| M-01 | Formula: GGGGC (4G,1C) | (4-1)/(4+1) = 0.6 | Formula |
| M-02 | Formula: CCCCC (0G,5C) | (0-5)/(0+5) = -1.0 | Formula |
| M-03 | Formula: GGGGG (5G,0C) | (5-0)/(5+0) = +1.0 | Formula |
| M-04 | Formula: GCGC (2G,2C) | (2-2)/(2+2) = 0.0 | Formula |
| M-05 | No G/C: AAATTT | 0 (no G+C) | Division protection |
| M-06 | Empty sequence | 0 | Implementation choice |
| M-07 | Invariant: all results in [-1, +1] | Range check | Wikipedia |
| M-08 | Windowed: correct positions at window centers | Positions correct | Grigoriev 1998 |
| M-09 | Cumulative: accumulates window skews | Values sum correctly | Grigoriev 1998 |
| M-10 | Cumulative: GGGG→CCCC pattern shows expected accumulation | 1→0 pattern | Grigoriev 1998 |
| M-11 | Origin detection: minimum at expected position | Min = origin | Grigoriev 1998 |
| M-12 | Terminus detection: maximum at expected position | Max = terminus | Grigoriev 1998 |
| M-13 | Null sequence → ArgumentNullException | Exception | Defensive |
| M-14 | Window size ≤ 0 → ArgumentOutOfRangeException | Exception | Defensive |
| M-15 | Step size ≤ 0 → ArgumentOutOfRangeException | Exception | Defensive |
| M-16 | Case insensitivity: lowercase handled | Same result | Implementation |
| M-17 | Sequence shorter than window → empty result | No windows | Logic |

### 4.2 SHOULD Tests

| Test ID | Scenario | Expected |
|---------|----------|----------|
| S-01 | AT skew formula verification | (A-T)/(A+T) correct |
| S-02 | Windowed overlapping windows produce more points | More points |
| S-03 | GC analysis result contains all metrics | All fields populated |

### 4.3 COULD Tests

| Test ID | Scenario | Expected |
|---------|----------|----------|
| C-01 | Large sequence performance | Completes in reasonable time |
| C-02 | Windowed with various step sizes | Correct point counts |

---

## 5. Edge Cases

| Category | Case | Expected Behavior | Source |
|----------|------|-------------------|--------|
| Empty | Empty string | Return 0 | ASSUMPTION |
| Empty | Windowed on empty | Return empty collection | ASSUMPTION |
| Null | Null DnaSequence | Throw ArgumentNullException | Defensive |
| Boundary | All G | +1.0 | Formula |
| Boundary | All C | -1.0 | Formula |
| Boundary | No G or C | 0.0 | Division protection |
| Window | Window > sequence length | Return empty | Logic |
| Window | Window = sequence length | Return single point | Logic |

---

## 6. Test Consolidation Plan

### 6.1 Current State

- **Canonical test file:** `GcSkewCalculatorTests.cs` (356 lines, 36 tests)
- **Coverage:** Good basic coverage of all methods
- **Issues:** None — all gaps closed
  - ~~Missing explicit invariant range test~~ ✅ Added

### 6.2 Consolidation Actions

| Action | Details |
|--------|---------|
| Enhance | Add invariant-focused assertions using Assert.Multiple |
| Add | Property-based range invariant test |
| Add | Formula verification tests with explicit calculations |
| Organize | Group tests by method and concern |
| Clean | Ensure naming follows `Method_Scenario_ExpectedResult` |

---

## 7. ASSUMPTIONS

| ID | Assumption | Rationale |
|----|------------|-----------|
| A-01 | Empty sequence returns 0 | Consistent with division-by-zero protection |
| A-02 | Case-insensitive processing | Standard bioinformatics practice |

---

## 8. Open Questions

None - all behaviors are well-defined by sources or implementation.

---

## 9. Validation Checklist

- [x] TestSpec created
- [x] Evidence sources documented
- [x] Must/Should/Could tests defined
- [x] Edge cases enumerated
- [x] Invariants documented
- [x] ASSUMPTIONS marked
- [x] Consolidation plan defined
