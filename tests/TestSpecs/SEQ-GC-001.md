# Test Specification: SEQ-GC-001

**Test Unit ID:** SEQ-GC-001
**Area:** Composition
**Algorithm:** GC Content Calculation
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-01-22

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: GC-content | https://en.wikipedia.org/wiki/GC-content | 2026-01-22 |
| Biopython Bio.SeqUtils | https://biopython.org/docs/latest/api/Bio.SeqUtils.html | 2026-01-22 |
| Madigan & Martinko (2003) | Brock Biology of Microorganisms, 10th ed. | Reference |

### 1.2 Formula (Wikipedia)

**GC Percentage:**
$$\frac{G + C}{A + T + G + C} \times 100\%$$

**GC Fraction (Biopython):**
$$\frac{G + C}{A + T + G + C}$$ (returns value 0-1)

### 1.3 Edge Cases from Evidence

| Edge Case | Expected Behavior | Source |
|-----------|-------------------|--------|
| Empty sequence | Return 0 | Biopython: "Note that this will return zero for an empty sequence" |
| All G/C | 100% (or 1.0 for fraction) | Mathematical derivation |
| All A/T | 0% (or 0.0 for fraction) | Mathematical derivation |
| Mixed case | Case-insensitive counting | Biopython: "Copes with mixed case sequences" |
| Ambiguous nucleotides (N, etc.) | Implementation-defined | Biopython offers 3 modes: remove, ignore, weighted |
| S and W ambiguity codes | Count as GC or ignore | Biopython: "S and W are ambiguous for GC content" |

### 1.4 Known Failure Modes

1. **Division by zero** - Must handle empty sequences
2. **Case sensitivity** - Must count both 'G'/'g' and 'C'/'c'
3. **Non-DNA characters** - Must define behavior (ignore or include in denominator)

### 1.5 Biological Context

- Human genome: 35-60% GC (mean ~41%) [Wikipedia, ref 20]
- Yeast (S. cerevisiae): 38% GC [Wikipedia, ref 21]
- Plasmodium falciparum: ~20% GC (extremely AT-rich) [Wikipedia, ref 23]
- Streptomyces coelicolor: 72% GC [Wikipedia, ref 29]

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateGcContent(ReadOnlySpan<char>)` | SequenceExtensions | **Canonical** | Returns percentage 0-100 |
| `CalculateGcFraction(ReadOnlySpan<char>)` | SequenceExtensions | **Canonical** | Returns fraction 0-1 |
| `CalculateGcContentFast(string)` | SequenceExtensions | Delegate | Wraps Span version |
| `CalculateGcFractionFast(string)` | SequenceExtensions | Delegate | Wraps Span version |
| `GcContent()` | DnaSequence | Delegate | Wraps CalculateGcContentFast |

---

## 3. Invariants

| ID | Invariant | Verifiable |
|----|-----------|------------|
| INV-1 | 0 ≤ CalculateGcContent ≤ 100 | Yes |
| INV-2 | 0 ≤ CalculateGcFraction ≤ 1 | Yes |
| INV-3 | CalculateGcContent = CalculateGcFraction × 100 | Yes |
| INV-4 | GcContent(lowercase) = GcContent(uppercase) | Yes |
| INV-5 | Empty sequence → 0 | Yes |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Empty sequence returns 0 | `""` | 0.0 | Biopython docs |
| M2 | All GC returns 100% | `"GCGCGC"` | 100.0 | Formula |
| M3 | All AT returns 0% | `"ATATAT"` | 0.0 | Formula |
| M4 | Equal ACGT returns 50% | `"ACGT"` | 50.0 | Formula |
| M5 | Mixed case handling | `"acgt"` vs `"ACGT"` | Same result | Biopython |
| M6 | Fraction matches percentage/100 | any | GcFraction = GcContent/100 | Formula |
| M7 | Boundary: Single G returns 100% | `"G"` | 100.0 | Formula |
| M8 | Boundary: Single A returns 0% | `"A"` | 0.0 | Formula |
| M9 | Biological: Human-like (~50%) | `"ATGCATGC"` | 50.0 | Formula |
| M10 | Biological: High GC (~75%) | `"GCGC"` | 100.0 | Formula |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | Long sequence accuracy | 1000 nt, 500 G/C | 50.0 | Formula |
| S2 | All G, no C | `"GGGG"` | 100.0 | Formula |
| S3 | All C, no G | `"CCCC"` | 100.0 | Formula |
| S4 | Delegate methods match canonical | Same input | Same result | Implementation |
| S5 | DnaSequence.GcContent() matches | Same input | Same result | Implementation |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| C1 | Very long sequence (10K) | 10000 nt | Correct % | Performance |
| C2 | Floating point precision | Edge ratios | No precision loss | Numeric |

---

## 5. Audit of Existing Tests

### 5.1 Final State After Consolidation

| File | Tests | Role |
|------|-------|------|
| `SequenceExtensions_CalculateGcContent_Tests.cs` | 54 | Canonical (all logic) |
| `DnaSequenceTests.cs` | 0 (comment only) | Delegates to canonical |
| `RnaSequenceTests.cs` | 0 (comment only) | Delegates to canonical |

### 5.2 Removed Duplicates

| Test Method | Original File | Reason |
|-------------|---------------|--------|
| `GcContent_AllGC_Returns100` | DnaSequenceTests.cs | Duplicate of canonical M2 |
| `GcContent_NoGC_Returns0` | DnaSequenceTests.cs | Duplicate of canonical M3 |
| `GcContent_HalfGC_Returns50` | DnaSequenceTests.cs | Duplicate of canonical M4 |
| `GcContent_Empty_Returns0` | DnaSequenceTests.cs | Duplicate of canonical M1 |
| `GcContent_AllGC_Returns100` | RnaSequenceTests.cs | Duplicate of canonical M2 |
| `GcContent_NoGC_Returns0` | RnaSequenceTests.cs | Duplicate of canonical M3 |
| `GcContent_HalfGC_Returns50` | RnaSequenceTests.cs | Duplicate of canonical M4 |
| `GcContent_Empty_Returns0` | RnaSequenceTests.cs | Duplicate of canonical M1 |

### 5.3 All Gaps Addressed

| Gap | Resolution |
|-----|------------|
| No direct tests for `CalculateGcContent(ReadOnlySpan<char>)` | ✓ Added in canonical file |
| No tests for `CalculateGcFraction` methods | ✓ Added in canonical file |
| No invariant relationship tests | ✓ Added in canonical file |
| No mixed case tests | ✓ Added in canonical file |
| No single-nucleotide boundary tests | ✓ Added in canonical file |

---

## 6. Test Implementation Plan

### Test File Structure (Standard Pattern)

```
┌─────────────────────────────────────────────────────────────────┐
│  Canonical Tests (deep, evidence-based)                        │
│  File: SequenceExtensions_CalculateGcContent_Tests.cs          │
│  - All MUST/SHOULD/COULD tests                                 │
│  - Invariant verification                                      │
│  - Delegate/wrapper smoke tests (DnaSequence, RnaSequence)     │
└─────────────────────────────────────────────────────────────────┘
                            │
                            ▼
┌─────────────────────────────────────────────────────────────────┐
│  Domain Class Tests (minimal, delegate verification only)       │
│  Files: DnaSequenceTests.cs, RnaSequenceTests.cs               │
│  - Comment reference to canonical tests                         │
│  - No duplicate logic tests                                     │
└─────────────────────────────────────────────────────────────────┘
```

### Removed Duplicates

| File | Removed Tests | Reason |
|------|---------------|--------|
| DnaSequenceTests.cs | GcContent_AllGC_Returns100, GcContent_NoGC_Returns0, GcContent_HalfGC_Returns50, GcContent_Empty_Returns0 | Covered in canonical file |
| RnaSequenceTests.cs | GcContent_AllGC_Returns100, GcContent_NoGC_Returns0, GcContent_HalfGC_Returns50, GcContent_Empty_Returns0 | Covered in canonical file |

### Consolidation Principle

1. **Canonical tests** test the core algorithm thoroughly
2. **Wrapper/delegate tests** verify only that delegation works correctly
3. **No duplication** - each behavior tested in exactly one place

---

## 7. Open Questions / Decisions

| Question | Decision | Rationale |
|----------|----------|-----------|
| How to handle N and other IUPAC codes? | **Not counted** (current impl) | Matches simplest Biopython mode |
| Return type double precision? | Standard double | Sufficient for biological data |

---

## 8. ASSUMPTIONS

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Non-ACGT characters are not counted toward GC | Implementation does not explicitly exclude, but only counts G/C |
| A2 | Return 0 for empty (not NaN or exception) | Matches Biopython behavior |

---

## 9. Validation Checklist

- [x] Evidence documented with sources
- [x] All MUST tests have evidence backing
- [x] Invariants identified
- [x] Existing tests audited
- [x] Gaps identified
- [x] Implementation plan created
