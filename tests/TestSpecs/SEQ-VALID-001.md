# Test Specification: SEQ-VALID-001

**Test Unit ID:** SEQ-VALID-001
**Area:** Composition
**Algorithm:** Sequence Validation
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-01-22

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Nucleic acid notation | https://en.wikipedia.org/wiki/Nucleic_acid_notation | 2026-01-22 |
| IUPAC-IUB Commission on Biochemical Nomenclature (1970) | doi:10.1021/bi00822a023 | Reference |
| Bioinformatics.org: IUPAC codes | https://www.bioinformatics.org/sms/iupac.html | 2026-01-22 |
| Biopython Seq class | https://biopython.org/wiki/Seq | 2026-01-22 |

### 1.2 IUPAC Standard Nucleotide Codes (Wikipedia, Bioinformatics.org)

**Standard DNA Bases:**
| Symbol | Base |
|--------|------|
| A | Adenine |
| C | Cytosine |
| G | Guanine |
| T | Thymine |

**Standard RNA Bases:**
| Symbol | Base |
|--------|------|
| A | Adenine |
| C | Cytosine |
| G | Guanine |
| U | Uracil |

**Ambiguity Codes (IUPAC extended):**
R, Y, S, W, K, M, B, D, H, V, N (any), - (gap)

### 1.3 Edge Cases from Evidence

| Edge Case | Expected Behavior | Source |
|-----------|-------------------|--------|
| Empty sequence | **ASSUMPTION**: Return `true` (empty sequence has no invalid characters) | Implementation choice |
| Standard bases only (A, C, G, T) | DNA: valid, RNA: valid only if ACGU | IUPAC 1970 |
| Case insensitivity | Both 'a' and 'A' are valid | Common practice, Biopython |
| U in DNA context | **Implementation**: Invalid for DNA | Implementation choice |
| T in RNA context | **Implementation**: Invalid for RNA | IUPAC standard |
| Ambiguity codes (N, R, Y, etc.) | **Implementation**: Invalid for strict validation | Implementation choice |
| Whitespace | **Implementation**: Invalid | Implementation choice |
| Numeric characters | Invalid | IUPAC 1970 |

### 1.4 Implementation Design Decisions

The current implementation uses **strict validation**:
- DNA: Only A, C, G, T (case-insensitive)
- RNA: Only A, C, G, U (case-insensitive)
- No ambiguity codes accepted
- Empty sequences return `true` (vacuously valid - no invalid characters exist)

This is a deliberate design choice for genomic analysis where strict validation is preferred over permissive validation.

### 1.5 Known Failure Modes

1. **Case sensitivity bugs** - Must check both upper and lower case
2. **Boundary confusion** - T valid for DNA but not RNA; U valid for RNA but not DNA
3. **Empty sequence handling** - Must define behavior explicitly

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `IsValidDna(ReadOnlySpan<char>)` | SequenceExtensions | **Canonical** | Returns true if all chars are A/C/G/T |
| `IsValidRna(ReadOnlySpan<char>)` | SequenceExtensions | **Canonical** | Returns true if all chars are A/C/G/U |
| `TryCreate(string, out DnaSequence)` | DnaSequence | Factory | Wraps validation + construction |
| `DnaSequence(string)` constructor | DnaSequence | Constructor | Throws on invalid input |

---

## 3. Invariants

| ID | Invariant | Verifiable |
|----|-----------|------------|
| INV-1 | IsValidDna(x) ⟹ all chars ∈ {A, C, G, T, a, c, g, t} | Yes |
| INV-2 | IsValidRna(x) ⟹ all chars ∈ {A, C, G, U, a, c, g, u} | Yes |
| INV-3 | IsValidDna(uppercase(x)) = IsValidDna(lowercase(x)) | Yes |
| INV-4 | IsValidRna(uppercase(x)) = IsValidRna(lowercase(x)) | Yes |
| INV-5 | TryCreate succeeds ⟺ IsValidDna returns true (for non-null input) | Yes |
| INV-6 | Empty string is valid (vacuously true) | Yes |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Empty sequence is valid DNA | `""` | true | ASSUMPTION: vacuous truth |
| M2 | Empty sequence is valid RNA | `""` | true | ASSUMPTION: vacuous truth |
| M3 | All standard DNA bases valid | `"ACGT"` | true | IUPAC 1970 |
| M4 | All standard RNA bases valid | `"ACGU"` | true | IUPAC 1970 |
| M5 | Lowercase DNA valid | `"acgt"` | true | Common practice |
| M6 | Lowercase RNA valid | `"acgu"` | true | Common practice |
| M7 | Mixed case DNA valid | `"AcGt"` | true | Common practice |
| M8 | U in DNA is invalid | `"ACGU"` | false | IUPAC: U is RNA only |
| M9 | T in RNA is invalid | `"ACGT"` | false | IUPAC: T is DNA only |
| M10 | Invalid character X | `"ACGX"` | false | IUPAC: X not standard |
| M11 | Numeric character invalid | `"ACG1"` | false | IUPAC 1970 |
| M12 | Whitespace invalid | `"AC GT"` | false | IUPAC 1970 |
| M13 | N (ambiguity) invalid for strict | `"ACGN"` | false | Implementation: strict mode |
| M14 | Single valid base A | `"A"` | true | IUPAC 1970 |
| M15 | Single invalid base X | `"X"` | false | IUPAC 1970 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | Long valid DNA sequence | 1000+ chars | true | Performance validation |
| S2 | Invalid char at start | `"XACGT"` | false | Boundary position |
| S3 | Invalid char at end | `"ACGTX"` | false | Boundary position |
| S4 | Invalid char in middle | `"ACXGT"` | false | Boundary position |
| S5 | All same valid base | `"AAAA"` | true | Degenerate case |
| S6 | Special chars (!@#) invalid | `"AC@T"` | false | Non-alpha invalid |

### 4.3 COULD Tests (Additional coverage)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| C1 | Unicode characters invalid | `"ACG日"` | false | Non-ASCII |
| C2 | Tab character invalid | `"AC\tGT"` | false | Whitespace variant |
| C3 | Newline invalid | `"AC\nGT"` | false | Whitespace variant |

---

## 5. DnaSequence Factory Tests

### 5.1 MUST Tests (TryCreate / Constructor)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| F1 | TryCreate with valid returns true | `"ACGT"` | (true, DnaSequence) | Factory pattern |
| F2 | TryCreate with invalid returns false | `"ACGX"` | (false, null) | Factory pattern |
| F3 | Constructor with valid creates instance | `"ACGT"` | DnaSequence | Normal construction |
| F4 | Constructor with invalid throws | `"ACGX"` | ArgumentException | Validation failure |
| F5 | Constructor normalizes to uppercase | `"acgt"` | Sequence = "ACGT" | Case normalization |
| F6 | Empty string creates empty sequence | `""` | DnaSequence (empty) | Edge case |

---

## 6. Audit & Consolidation

### 6.1 Existing Tests Found

| Location | Test Methods | Status | Action |
|----------|--------------|--------|--------|
| PerformanceExtensionsTests.cs | IsValidDna_ValidSequence_ReturnsTrue | Weak (single case) | Move to canonical |
| PerformanceExtensionsTests.cs | IsValidDna_InvalidSequence_ReturnsFalse | Weak (single case) | Move to canonical |
| PerformanceExtensionsTests.cs | IsValidRna_ValidSequence_ReturnsTrue | Weak (single case) | Move to canonical |
| DnaSequenceTests.cs | TryCreate_ValidSequence_ReturnsTrue | Basic coverage | Keep as smoke |
| DnaSequenceTests.cs | TryCreate_InvalidSequence_ReturnsFalse | Basic coverage | Keep as smoke |
| DnaSequenceTests.cs | Constructor_InvalidNucleotide_ThrowsArgumentException | Basic coverage | Keep as smoke |
| DnaSequenceTests.cs | Constructor_EmptySequence_CreatesEmpty | Covered | Keep as smoke |

### 6.2 Consolidation Plan

1. **Create canonical test file**: `SequenceExtensions_IsValidDna_Tests.cs`
   - Contains all deep validation tests (M1-M15, S1-S6, C1-C3)
   - Tests both `IsValidDna` and `IsValidRna` canonical methods

2. **PerformanceExtensionsTests.cs**: Remove validation tests (moved to canonical)

3. **DnaSequenceTests.cs**: Keep factory/constructor tests as smoke tests (1-2 per method)
   - These delegate to `IsValidDna` internally, so only need smoke verification

---

## 7. ASSUMPTIONS Log

| ID | Assumption | Justification |
|----|------------|---------------|
| A1 | Empty sequence returns `true` | Vacuous truth: no invalid characters exist. Common library behavior. |
| A2 | Strict validation (no IUPAC ambiguity codes) | Design decision for genomic analysis accuracy |
| A3 | Case-insensitive validation | Universal practice in bioinformatics |

---

## 8. Open Questions / Decisions

| # | Question | Resolution |
|---|----------|------------|
| 1 | Should empty return true or false? | Resolved: TRUE (vacuous truth, matches implementation) |
| 2 | Should IUPAC ambiguity codes be accepted? | Resolved: NO (strict validation mode) |
