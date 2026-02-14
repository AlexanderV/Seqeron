# Test Specification: SEQ-VALID-001

**Test Unit ID:** SEQ-VALID-001
**Area:** Composition
**Algorithm:** Sequence Validation
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Nucleic acid notation | https://en.wikipedia.org/wiki/Nucleic_acid_notation | 2026-02-14 |
| IUPAC-IUB Commission on Biochemical Nomenclature (1970) | doi:10.1021/bi00822a023 | Reference |
| NC-IUB: Nomenclature for Incompletely Specified Bases (1984) | doi:10.1093/nar/13.9.3021 | Reference |
| Bioinformatics.org: IUPAC codes | https://www.bioinformatics.org/sms/iupac.html | 2026-02-14 |
| Biopython `Bio.Data.IUPACData` | https://github.com/biopython/biopython/blob/master/Bio/Data/IUPACData.py | 2026-02-14 |
| Biopython `Bio.Seq` API | https://biopython.org/docs/latest/api/Bio.Seq.html | 2026-02-14 |

### 1.2 IUPAC Standard Nucleotide Codes

**Unambiguous DNA Bases (IUPAC 1970; Biopython `unambiguous_dna_letters = "GATC"`):**
| Symbol | Base |
|--------|------|
| A | Adenine |
| C | Cytosine |
| G | Guanine |
| T | Thymine |

**Unambiguous RNA Bases (IUPAC 1970; Biopython `unambiguous_rna_letters = "GAUC"`):**
| Symbol | Base |
|--------|------|
| A | Adenine |
| C | Cytosine |
| G | Guanine |
| U | Uracil |

**Ambiguity Codes (NC-IUB 1984, for consensus/polymorphism notation):**
R, Y, S, W, K, M, B, D, H, V, N (any), - (gap)

> Wikipedia: "Degenerate base symbols [...] are an IUPAC representation for a position
> on a DNA sequence that can have multiple possible alternatives. These should not be
> confused with non-canonical bases because each particular sequence will have in fact
> one of the regular bases." — These codes encode uncertainty in consensus sequences,
> not actual nucleotide values.

### 1.3 Edge Cases from Evidence

| Edge Case | Expected Behavior | Source |
|-----------|-------------------|--------|
| Empty sequence | Return `true` (vacuous truth) | Biopython: "Zero-length sequences are always considered to be defined" |
| Standard bases only (A, C, G, T) | DNA: valid | IUPAC 1970; Biopython `unambiguous_dna_letters` |
| Standard bases only (A, C, G, U) | RNA: valid | IUPAC 1970; Biopython `unambiguous_rna_letters` |
| Case insensitivity | Both 'a' and 'A' are valid | Wikipedia: "Lowercase versions of the IUPAC letters are used in genetic sequence files" |
| U in DNA context | Invalid for DNA | IUPAC 1970: U is RNA nucleoside |
| T in RNA context | Invalid for RNA | IUPAC 1970: T is DNA nucleoside |
| Ambiguity codes (N, R, Y, etc.) | Invalid (not actual bases) | NC-IUB 1984: represent positional variants, not nucleotides |
| Whitespace, numerics, special chars | Invalid | IUPAC 1970: not part of nucleotide notation |

### 1.4 Validation Semantics

The methods validate whether every character in the sequence is an **unambiguous nucleotide**
of the specified type:

- `IsValidDna`: char ∈ {A, C, G, T} (case-insensitive) — matches Biopython `unambiguous_dna_letters`
- `IsValidRna`: char ∈ {A, C, G, U} (case-insensitive) — matches Biopython `unambiguous_rna_letters`

This is the standard validation for actual sequence data. IUPAC ambiguity codes (R, Y, N, etc.)
represent uncertainty in consensus sequences and are correctly rejected by unambiguous validation.

### 1.5 Known Failure Modes

1. **Case sensitivity bugs** — Must normalize both upper and lower case
2. **Boundary confusion** — T valid for DNA but not RNA; U valid for RNA but not DNA

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `IsValidDna(ReadOnlySpan<char>)` | SequenceExtensions | **Canonical** | Returns true if all chars ∈ {A,C,G,T} |
| `IsValidRna(ReadOnlySpan<char>)` | SequenceExtensions | **Canonical** | Returns true if all chars ∈ {A,C,G,U} |
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
| INV-6 | Empty string is valid (vacuous truth — Biopython: "Zero-length sequences are always considered to be defined") | Yes |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Empty sequence is valid DNA | `""` | true | Biopython: zero-length sequences are defined; vacuous truth |
| M2 | Empty sequence is valid RNA | `""` | true | Biopython: zero-length sequences are defined; vacuous truth |
| M3 | All standard DNA bases valid | `"ACGT"` | true | IUPAC 1970; Biopython `unambiguous_dna_letters` |
| M4 | All standard RNA bases valid | `"ACGU"` | true | IUPAC 1970; Biopython `unambiguous_rna_letters` |
| M5 | Lowercase DNA valid | `"acgt"` | true | Wikipedia: lowercase is standard in sequence files |
| M6 | Lowercase RNA valid | `"acgu"` | true | Wikipedia: lowercase is standard in sequence files |
| M7 | Mixed case DNA valid | `"AcGt"` | true | Wikipedia: both cases represent same nucleotides |
| M8 | U in DNA is invalid | `"ACGU"` | false | IUPAC 1970: U is RNA nucleoside |
| M9 | T in RNA is invalid | `"ACGT"` | false | IUPAC 1970: T is DNA nucleoside |
| M10 | Invalid character X | `"ACGX"` | false | IUPAC: X is not a nucleotide symbol |
| M11 | Numeric character invalid | `"ACG1"` | false | IUPAC 1970: not part of nucleotide notation |
| M12 | Whitespace invalid | `"AC GT"` | false | IUPAC 1970: not part of nucleotide notation |
| M13 | N (ambiguity) invalid | `"ACGN"` | false | NC-IUB 1984: N is ambiguity code, not a nucleotide |
| M14 | Single valid base A | `"A"` | true | IUPAC 1970 — covered by EachValidBase(A,C,G,T) |
| M15 | Single invalid base X | `"X"` | false | IUPAC: not a nucleotide |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | Long valid DNA sequence | 1000+ chars | true | Boundary/performance |
| S2 | Invalid char at start | `"XACGT"` | false | Boundary position |
| S3 | Invalid char at end | `"ACGTX"` | false | Boundary position |
| S4 | Invalid char in middle | `"ACXGT"` | false | Boundary position |
| S5 | All same valid base | `"AAAA"` | true | Degenerate case |
| S6 | Special chars (!@#) invalid | `"AC@T"` | false | Not nucleotide notation |

### 4.3 COULD Tests (Additional coverage)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
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
| F6 | Empty string creates empty sequence | `""` | DnaSequence (empty) | Vacuous validity |

---

## 6. Coverage Classification

### 6.1 Summary

| Metric | Value |
|--------|-------|
| Total canonical test runs | 46 |
| FsCheck property tests | 2 |
| Total | 48 |

### 6.2 Classification (canonical file)

| Test Method | Runs | Spec IDs | Status | Notes |
|-------------|------|----------|--------|-------|
| `IsValidDna_EmptySequence_ReturnsTrue` | 1 | M1 | ✅ Covered | |
| `IsValidDna_AllStandardBases_ReturnsTrue` | 1 | M3 | ✅ Covered | |
| `IsValidDna_LowercaseBases_ReturnsTrue` | 1 | M5 | ✅ Covered | |
| `IsValidDna_MixedCase_ReturnsTrue` | 1 | M7 | ✅ Covered | |
| `IsValidDna_ContainsUracil_ReturnsFalse` | 1 | M8 | ✅ Covered | |
| `IsValidDna_InvalidCharacterX_ReturnsFalse` | 1 | M10 | ✅ Covered | |
| `IsValidDna_NumericCharacter_ReturnsFalse` | 1 | M11 | ✅ Covered | |
| `IsValidDna_Whitespace_ReturnsFalse` | 1 | M12 | ✅ Covered | |
| `IsValidDna_AmbiguityCodeN_ReturnsFalse` | 1 | M13 | ✅ Covered | |
| `IsValidDna_SingleInvalidBase_ReturnsFalse` | 1 | M15 | ✅ Covered | |
| `IsValidDna_EachValidBase_ReturnsTrue` | 4 | M14, S-bases | ✅ Covered | A, C, G, T — subsumes M14 |
| `IsValidDna_IupacAmbiguityCodes_ReturnsFalse` | 10 | M13-ext | ✅ Covered | R,Y,S,W,K,M,B,D,H,V |
| `IsValidDna_LongValidSequence_ReturnsTrue` | 1 | S1 | ✅ Covered | |
| `IsValidDna_InvalidAtStart_ReturnsFalse` | 1 | S2 | ✅ Covered | |
| `IsValidDna_InvalidAtEnd_ReturnsFalse` | 1 | S3 | ✅ Covered | |
| `IsValidDna_InvalidInMiddle_ReturnsFalse` | 1 | S4 | ✅ Covered | |
| `IsValidDna_AllSameBase_ReturnsTrue` | 1 | S5 | ✅ Covered | |
| `IsValidDna_SpecialCharacters_ReturnsFalse` | 1 | S6 | ✅ Covered | |
| `IsValidDna_UnicodeCharacter_ReturnsFalse` | 1 | C1 | ✅ Covered | |
| `IsValidDna_TabCharacter_ReturnsFalse` | 1 | C2 | ✅ Covered | |
| `IsValidDna_NewlineCharacter_ReturnsFalse` | 1 | C3 | ✅ Covered | |
| `IsValidDna_GapCharacter_ReturnsFalse` | 1 | — | ✅ Covered | |
| `IsValidRna_EmptySequence_ReturnsTrue` | 1 | M2 | ✅ Covered | |
| `IsValidRna_AllStandardBases_ReturnsTrue` | 1 | M4 | ✅ Covered | |
| `IsValidRna_LowercaseBases_ReturnsTrue` | 1 | M6 | ✅ Covered | |
| `IsValidRna_ContainsThymine_ReturnsFalse` | 1 | M9 | ✅ Covered | |
| `IsValidRna_AmbiguityCodeN_ReturnsFalse` | 1 | M13-RNA | ✅ Covered | |
| `IsValidRna_EachValidBase_ReturnsTrue` | 4 | S-bases-RNA | ✅ Covered | A, C, G, U |
| `IsValidRna_LongValidSequence_ReturnsTrue` | 1 | S1-RNA | ✅ Covered | |
| `IsValidDna_CaseInvariance_AllCasesReturnTrue` | 1 | INV-3 | ✅ Covered | Asserts true for upper/lower/mixed |
| `IsValidRna_CaseInvariance_AllCasesReturnTrue` | 1 | INV-4 | ✅ Covered | Asserts true for upper/lower/mixed |

### 6.3 Classification (FsCheck property tests)

| Test Method | Spec ID | Status | Notes |
|-------------|---------|--------|-------|
| `PureAcgt_IsValidDna` | INV-1 | ✅ Covered | Random ACGT strings → always valid DNA |
| `PureAcgu_IsValidRna` | INV-2 | ✅ Covered | Random ACGU strings → always valid RNA |

### 6.4 Changes Applied

| Action | Test | Reason |
|--------|------|--------|
| ❌→✅ Added | `IsValidDna_UnicodeCharacter_ReturnsFalse` | C1 was missing |
| ⚠→✅ Strengthened | `IsValidDna_CaseInvariance_AllCasesReturnTrue` | Was: equality check only; Now: asserts `Is.True` for each case |
| ⚠→✅ Strengthened | `IsValidRna_CaseInvariance_AllCasesReturnTrue` | Was: equality check only; Now: asserts `Is.True` for each case |
| 🔁 Removed | `IsValidDna_SingleValidBase_ReturnsTrue` (1 run) | Duplicate of `EachValidBase(TestCase "A")` |
| 🔁 Removed | `PerformanceExtensionsTests.IsValidDna_ValidSequence_SmokeTest` (1 run) | Duplicate of M3 |
| 🔁 Removed | `PerformanceExtensionsTests.IsValidRna_ValidSequence_SmokeTest` (1 run) | Duplicate of M4 |
| 🔁 Removed | `SequenceCompositionProperties.InvalidChars_NotValidDna` (3 runs) | Duplicate of M11, M12 |

---

## 7. Deviations and Assumptions

None. All behaviors are sourced from IUPAC standards and Biopython reference implementation.
