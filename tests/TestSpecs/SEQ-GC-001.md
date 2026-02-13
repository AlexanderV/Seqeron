# Test Specification: SEQ-GC-001

**Test Unit ID:** SEQ-GC-001
**Area:** Composition
**Algorithm:** GC Content Calculation
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: GC-content | https://en.wikipedia.org/wiki/GC-content | 2026-02-14 |
| Biopython Bio.SeqUtils.gc_fraction | https://biopython.org/docs/latest/api/Bio.SeqUtils.html | 2026-02-14 |
| Biopython source (Bio/SeqUtils/__init__.py) | https://github.com/biopython/biopython/blob/main/Bio/SeqUtils/__init__.py | 2026-02-14 |
| Madigan & Martinko (2003) | Brock Biology of Microorganisms, 10th ed. | Reference |

### 1.2 Formula (Wikipedia, ref 7)

**GC Percentage:**
$$\frac{G + C}{A + T + G + C} \times 100\%$$

**GC Fraction (Biopython `gc_fraction`):**
$$\frac{G + C}{A + T + G + C}$$ (returns value 0-1)

**Denominator**: count of valid nucleotides only (A, T, G, C, U). Non-nucleotide characters are excluded from **both** numerator and denominator. This matches Biopython `gc_fraction(seq, "remove")` (default mode).

### 1.3 Denominator Semantics (Biopython default "remove" mode)

Biopython's `gc_fraction` default behavior (source verified):

```python
gc = sum(seq.count(x) for x in "CGScgs")
length = gc + sum(seq.count(x) for x in "ATWUatwu")  # only valid nucleotides
```

Our implementation equivalently counts **A, T, G, C, U** (case-insensitive) for the denominator. Characters outside this set (e.g., N, R, Y, B, D, H, K, M, V, X) are excluded from both numerator and denominator.

Note: Biopython also counts S (Strong = G|C) as GC and W (Weak = A|T) as AT. Our implementation does not count S/W since they are not standard nucleotides per the Wikipedia formula. For standard DNA/RNA sequences this produces identical results.

### 1.4 Edge Cases — Defined Behavior

| Case | Expected | Source |
|------|----------|--------|
| Empty sequence | Return 0 | Biopython: "Note that this will return zero for an empty sequence" |
| All G/C | 100% (or 1.0) | Formula: (G+C)/(G+C) = 1 |
| All A/T | 0% (or 0.0) | Formula: 0/(A+T) = 0 |
| Mixed case | Case-insensitive | Biopython: "Copes with mixed case sequences" |
| Contains N/ambiguous chars | Excluded from calculation | Wikipedia formula + Biopython "remove" mode |
| Only non-nucleotide chars | Return 0 | No valid nucleotides → 0 (same as empty) |
| Contains U (RNA) | U counted as valid nucleotide (not GC) | Biopython: U included in valid set |

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
| INV-6 | Complement preserves GC content | Yes |
| INV-7 | Reverse complement preserves GC content | Yes |

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
| M7 | Single G returns 100% | `"G"` | 100.0 | Formula |
| M8 | Single A returns 0% | `"A"` | 0.0 | Formula |
| M9 | N excluded from denominator | `"ACTGN"` | 50.0 | Wikipedia formula, Biopython "remove" |
| M10 | Only non-nucleotides returns 0 | `"NNNNN"` | 0.0 | No valid nucleotides |
| M11 | RNA uracil is valid nucleotide | `"GCAU"` | 50.0 | Biopython: U in valid set |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | Long sequence accuracy | 1000 nt, 500 G/C | 50.0 | Formula |
| S2 | All G, no C | `"GGGG"` | 100.0 | Formula |
| S3 | All C, no G | `"CCCC"` | 100.0 | Formula |
| S4 | Delegate methods match canonical | Same input | Same result | Implementation |
| S5 | DnaSequence.GcContent() matches | Same input | Same result | Implementation |
| S6 | Multiple N excluded | `"CCTGNN"` | 75.0 | Biopython: gc_fraction("CCTGNN","remove")=0.75 |
| S7 | All valid are GC + ambiguous | `"GCNN"` | 100.0 | Formula |
| S8 | Biopython GDVV example | `"GDVV"` | 100.0 | Biopython: gc_fraction("GDVV","remove")=1.00 |
| S9 | RNA sequence | `"GGAUCUUCGGAUCU"` | 50.0 | Biopython: gc_fraction=0.50 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| C1 | Very long sequence (10K) | 10000 nt | Correct % | Performance |
| C2 | Floating point precision | Edge ratios | No precision loss | Numeric |

---

## 5. Test File Structure

| File | Tests | Role |
|------|-------|------|
| `SequenceExtensions_CalculateGcContent_Tests.cs` | 90 runs | Canonical (all logic) |
| `Properties/GcContentProperties.cs` | 7 runs | Property-based invariants |
| `DnaSequenceTests.cs` | 0 (comment only) | Delegates to canonical |
| `RnaSequenceTests.cs` | 0 (comment only) | Delegates to canonical |

**Total: 97 test runs (was 118 before coverage classification — 21 duplicates removed, 2 new tests added)**

### 5.1 Coverage Classification Summary

| Action | Count | Details |
|--------|-------|---------|
| 🔁 Removed (duplicate) | 19 | `CalculateGcContentFast_EmptyString_ReturnsZero`, `CalculateGcFractionFast_EmptyString_ReturnsZero` (covered by delegate=canonical + empty canonical), `CalculateGcContentFast_MixedCase_MatchesUppercase` (covered by delegate=canonical + case tests), `AnyValidInput_ResultInRange0To100` 6 cases (covered by property INV-1), `AnyValidInput_ResultInRange0To1` 6 cases (covered by property INV-2), `NonNucleotideExclusion_MatchesFormula` 8 parametrized cases (duplicated individual tests) |
| ⚠ Strengthened (weak) | 5 | `LowercaseInput_MatchesUppercase` → `Returns50` (added exact value assertion), `MixedCaseInput_MatchesUppercase` → `Returns100` (added exact value), 3 Biological tests → parametrized with exact `Is.EqualTo` instead of misleading names |
| ❌ Added (missing) | 2 | `SingleU_ReturnsZero` (U boundary), `GcFraction_SequenceWithMultipleN_ExcludesAllN` (fraction CCTGNN) |
| ✅ Covered (unchanged) | 64 | All remaining tests verified as properly covering their spec requirement |
| 🔀 Merged | 3 | `AllG/AllC/MixedGC/AllA/AllT/MixedAT` (6 tests) → `HomogeneousSequences_ReturnsExtrema` (1 test, 6 cases) |

### Consolidation Principle

1. **Canonical tests** test the core algorithm thoroughly
2. **Wrapper/delegate tests** verify only that delegation works correctly
3. **No duplication** — each behavior tested in exactly one place

---

## 6. Biopython Cross-Verification Table

Values verified against Biopython `gc_fraction(seq, "remove")`:

| Input | Biopython result | Our result | Match |
|-------|-----------------|------------|-------|
| `""` | 0 | 0 | ✓ |
| `"ACTG"` | 0.50 | 50% / 0.50 | ✓ |
| `"ACTGN"` | 0.50 | 50% / 0.50 | ✓ |
| `"CCTGNN"` | 0.75 | 75% / 0.75 | ✓ |
| `"GDVV"` | 1.00 | 100% / 1.00 | ✓ |
| `"GCAU"` | 0.50 | 50% / 0.50 | ✓ |
| `"GGAUCUUCGGAUCU"` | 0.50 | 50% / 0.50 | ✓ |
| `"NNNNN"` | 0 | 0 | ✓ |

---

## 7. Validation Checklist

- [x] Evidence documented with sources
- [x] All MUST tests have evidence backing
- [x] Invariants identified and property-tested
- [x] Formula matches Wikipedia exactly: (G+C)/(A+T+G+C) × 100
- [x] Non-nucleotide handling matches Biopython default "remove" mode
- [x] Cross-verified against Biopython expected values
- [x] No assumptions — all behaviors sourced
- [x] Existing tests audited, no gaps
- [x] Coverage classification complete — no duplicates, no weak assertions
