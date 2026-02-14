# Test Specification: SEQ-REVCOMP-001 - Reverse Complement

**Test Unit ID:** SEQ-REVCOMP-001
**Area:** Composition
**Algorithm:** DNA Reverse Complement
**Canonical Method:** `SequenceExtensions.TryGetReverseComplement(ReadOnlySpan<char>, Span<char>)`
**Secondary API:** `DnaSequence.GetReverseComplementString(string)`
**Test File:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceExtensions_ReverseComplement_Tests.cs`

---

## 1. Purpose

Verify that the reverse complement operation produces a sequence that is the base-by-base complement in reverse order, per Watson-Crick base-pairing rules and IUPAC NC-IUB 1984 ambiguity code complements.

---

## 2. Authoritative Sources

| Source | Used For |
|--------|----------|
| Watson-Crick base pairing (all biochemistry textbooks) | A ↔ T, C ↔ G standard complement |
| IUPAC NC-IUB 1984 — Nucleic acid notation | Ambiguity code complement table (R ↔ Y, S ↔ S, W ↔ W, K ↔ M, B ↔ V, D ↔ H, N ↔ N) |
| Wikipedia — Nucleic acid sequence | "Complementary sequence to TTAC is GTAA" |
| Wikipedia — Complementarity (molecular biology) | Reverse complement definition, biological palindromes |
| Biopython Bio.Seq | Cross-verification: `reverse_complement()`, `complement()` examples |

---

## 3. IUPAC Complement Table (NC-IUB 1984)

| Base | Meaning | Complement | Complement Meaning |
|------|---------|------------|-------------------|
| A | Adenine | T | Thymine |
| T | Thymine | A | Adenine |
| G | Guanine | C | Cytosine |
| C | Cytosine | G | Guanine |
| U | Uracil | A | Adenine |
| R | A \| G (puRine) | Y | C \| T (pYrimidine) |
| Y | C \| T (pYrimidine) | R | A \| G (puRine) |
| S | C \| G (Strong) | S | C \| G (Strong) |
| W | A \| T (Weak) | W | A \| T (Weak) |
| K | G \| T (Keto) | M | A \| C (aMino) |
| M | A \| C (aMino) | K | G \| T (Keto) |
| B | C \| G \| T (not A) | V | A \| C \| G (not T) |
| D | A \| G \| T (not C) | H | A \| C \| T (not G) |
| H | A \| C \| T (not G) | D | A \| G \| T (not C) |
| V | A \| C \| G (not T) | B | C \| G \| T (not A) |
| N | Any nucleotide | N | Any nucleotide |

---

## 4. Requirements

### MUST

| ID | Requirement | Source |
|----|-------------|--------|
| MUST-01 | `ReverseComplement(s) = Reverse(Complement(s))` — complement each base, then reverse | Watson-Crick; Wikipedia Complementarity |
| MUST-02 | `ReverseComplement(ReverseComplement(s)) = s` — involution property holds | Mathematical property of complement + reverse |
| MUST-03 | Empty sequence → empty result (returns `true`) | Biopython: `reverse_complement("") → ""` |
| MUST-04 | Single nucleotide → its complement (A→T, T→A, G→C, C→G, U→A) | Watson-Crick |
| MUST-05 | Destination smaller than source → returns `false`, no partial writes | Span-based API safety contract |
| MUST-06 | Biological palindromes (EcoRI `GAATTC`, BamHI `GGATCC`, HindIII `AAGCTT`) equal their own reverse complement | Wikipedia: restriction enzyme palindromes |
| MUST-07 | Case-insensitive input, always uppercase output | IUPAC notation convention (uppercase); `DnaSequence` normalizes to uppercase |
| MUST-08 | U (uracil) complements to A via DNA rules (`GetComplementBase` is DNA-centric) | IUPAC: U pairs with A |
| MUST-09 | IUPAC ambiguity codes complemented per NC-IUB 1984 table (R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N) | IUPAC NC-IUB 1984; Biopython cross-verified |

### SHOULD

| ID | Requirement | Source |
|----|-------------|--------|
| SHOULD-02 | Correct for sequences ≥ 100 bases | General robustness |
| SHOULD-04 | Destination larger than source — only `source.Length` chars written | Span-based API contract |
| SHOULD-05 | Gap characters (`-`) pass through unchanged and position is reversed | Biopython: gaps pass through `complement()` |

---

## 5. Biopython Cross-Verification

All examples verified against Biopython 1.84 `Bio.Seq`:

| Input | Operation | Expected Output | Verified |
|-------|-----------|-----------------|----------|
| `CCCCCGATAGNR` | `reverse_complement()` | `YNCTATCGGGGG` | ✅ |
| `ACTG-NH` | `reverse_complement()` | `DN-CAGT` | ✅ |
| `ACTG-NH` | `complement()` | `TGAC-ND` | ✅ |
| `TTAC` | `reverse_complement()` | `GTAA` | ✅ (Wikipedia) |
| `GAATTC` | `reverse_complement()` | `GAATTC` | ✅ (palindrome) |
| `""` (empty) | `reverse_complement()` | `""` | ✅ |

---

## 6. Deviations from Textbook Definitions

| ID | Description | Justification |
|----|-------------|---------------|
| D1 | Output is always uppercase regardless of input case | IUPAC notation convention is uppercase; `DnaSequence` stores uppercase internally. Sourced to IUPAC standard. |
| D2 | `GetComplementBase` returns DNA bases for RNA input (U→A, not U→A with RNA output U) | `GetComplementBase` is DNA-centric by design. RNA sequences have a dedicated `GetRnaComplementBase` method. |

---

## 7. Test-to-Requirement Map

| Requirement | Tests |
|-------------|-------|
| MUST-01 | `TryGetReverseComplement_PalindromicSequence_ReturnsSameSequence`, `TryGetReverseComplement_AsymmetricSequence_ReturnsCorrectResult`, `TryGetReverseComplement_WikipediaExample_ReturnsGTAA`, `TryGetReverseComplement_LongerSequence_ReturnsCorrectResult` |
| MUST-02 | `TryGetReverseComplement_Involution_MultipleSequences` (8 sequences incl. ACGT, AACGTTAA) |
| MUST-03 | `TryGetReverseComplement_EmptySequence_ReturnsTrueAndBufferUntouched`, `TryGetReverseComplement_EmptySourceAndDestination_ReturnsTrue` |
| MUST-04 | `TryGetReverseComplement_AllSingleBases_ReturnComplements` (A, T, G, C, U) |
| MUST-05 | `TryGetReverseComplement_DestinationTooSmall_ReturnsFalseAndBufferUntouched`, `TryGetReverseComplement_EmptyDestinationNonEmptySource_ReturnsFalse` |
| MUST-06 | `TryGetReverseComplement_BiologicalPalindrome_IsOwnReverseComplement` ×3 (EcoRI, BamHI, HindIII) |
| MUST-07 | `TryGetReverseComplement_CaseInsensitive_ReturnsUppercase` ×3 (all-lower, mixed, asymmetric) |
| MUST-08 | `TryGetReverseComplement_RnaSequence_UsesDnaComplementRules`, `TryGetReverseComplement_RnaAsymmetric_ReturnsCorrectResult` |
| MUST-09 | `TryGetReverseComplement_IupacAmbiguityCodes_ComplementedCorrectly`, `TryGetReverseComplement_IupacInvolution_HoldsForAllCodes`, `TryGetReverseComplement_LowercaseIupac_ReturnsUppercase` (all 11 codes) |
| SHOULD-02 | `TryGetReverseComplement_LongAsymmetricSequence_ReturnsCorrectResult`, `TryGetReverseComplement_LongSequence_InvolutionHolds` |
| SHOULD-04 | `TryGetReverseComplement_DestinationLarger_WritesCorrectly` |
| SHOULD-05 | `TryGetReverseComplement_WithGap_PreservesAndReverses` |
| Cross-verify | `TryGetReverseComplement_BiopythonExample_CCCCCGATAGNR`, `TryGetReverseComplement_BiopythonExample_ACTG_NH`, `TryGetComplement_BiopythonExample_ACTG_NH` |
| Static API | `GetReverseComplementString_BasicSequence_ReturnsCorrectResult`, `GetReverseComplementString_EmptyString_ReturnsEmpty`, `GetReverseComplementString_Null_ReturnsNull`, `GetReverseComplementString_Palindrome_ReturnsSame` |

---

## 8. Coverage Summary

| Category | Status |
|----------|--------|
| Watson-Crick complement (A↔T, C↔G) | ✅ Sourced — textbook |
| IUPAC ambiguity codes (R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N) | ✅ Sourced — IUPAC NC-IUB 1984 |
| Involution property | ✅ Sourced — mathematical property |
| Empty sequence | ✅ Sourced — Biopython |
| Uppercase output | ✅ Sourced — IUPAC convention + DnaSequence API |
| RNA uracil handling | ✅ Sourced — IUPAC (U pairs with A) |
| Gap pass-through | ✅ Sourced — Biopython |
| Biological palindromes | ✅ Sourced — Wikipedia |
| Biopython cross-verification | ✅ 3 examples verified |
| Span-based API safety | ✅ Sourced — API contract |

---

## 9. Coverage Classification (2026-02-14)

**Canonical file:** `SequenceExtensions_ReverseComplement_Tests.cs`

| Metric | Before | After |
|--------|--------|-------|
| Test methods | 35 | 27 |
| Test runs (incl. TestCase) | 36 | 32 |

| Action | Count | Details |
|--------|-------|---------|
| 🔁 Duplicate removed | 1 | `SingleA_ReturnsT` — subsumed by `AllSingleBases_ReturnComplements` |
| 🔁 Duplicate removed | 1 | `DestinationExactSize_Succeeds` — identical to MUST-01 palindromic (ACGT, char[4]) |
| 🔁 Duplicate merged 2→0 | 2 | `Involution_ACGT` + `Involution_AsymmetricSequence` → absorbed into `Involution_MultipleSequences` |
| 🔁 Duplicate merged 3→1 | 2 | 3 palindrome tests → 1 `[TestCase]` parametrized |
| 🔁 Duplicate merged 3→1 | 2 | 3 MUST-07 case tests → 1 `[TestCase]` parametrized |
| ⚠ Weak strengthened | 1 | `EmptySequence`: pre-fill buffer + assert untouched |
| ⚠ Weak strengthened | 1 | `DestinationTooSmall`: pre-fill buffer + assert no partial writes |
| ⚠ Weak strengthened | 1 | `LowercaseIupac`: added 3 missing codes (d, h, v) — now all 11 |
| ⚠ Weak strengthened | 2 | Long sequence tests: replaced palindromic inputs with asymmetric sequences |
| ❌ Missing | 0 | — |
