# Test Specification: REP-PALIN-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | REP-PALIN-001 |
| **Area** | Repeats |
| **Title** | Palindrome Detection |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-01-22 |

---

## Methods Under Test

| Method | Class | Type | Test Priority |
|--------|-------|------|---------------|
| `FindPalindromes(DnaSequence, minLength, maxLength)` | RepeatFinder | Canonical | Deep testing |
| `FindPalindromes(string, minLength, maxLength)` | RepeatFinder | Overload | Smoke testing |
| `FindPalindromes(DnaSequence, minLength, maxLength)` | GenomicAnalyzer | Alternate | Smoke testing |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| [Wikipedia - Palindromic sequence](https://en.wikipedia.org/wiki/Palindromic_sequence) | Definition | A sequence is palindromic if it equals its reverse complement |
| [Wikipedia - Restriction enzyme](https://en.wikipedia.org/wiki/Restriction_enzyme) | Application | Type II restriction enzymes recognize palindromic sequences; most 4-8 bp; blunt or sticky ends |
| [Rosalind - REVP (Locating Restriction Sites)](https://rosalind.info/problems/revp/) | Test data | Provides sample dataset and expected output for palindrome detection |
| REBASE Database | Reference | Database of restriction enzymes and recognition sites |

---

## Key Definitions

### DNA Palindrome (Biological Definition)

A **DNA palindrome** (also called "reverse palindrome") is a sequence that reads the same 5' to 3' on both strands. This is equivalent to saying the sequence equals its reverse complement.

**Example: EcoRI recognition site (GAATTC)**
```
5'-GAATTC-3'
3'-CTTAAG-5'
```
Reading 5'→3' on top strand: GAATTC  
Reading 5'→3' on bottom strand: GAATTC (complement reversed)  
GAATTC == ReverseComplement(GAATTC) ✓

### Invariant

For any palindrome P:
```
P == ReverseComplement(P)
```

### Why Even Length Only

DNA palindromes in the biological sense must have **even length** because:
1. Each position has a complementary position on the opposite strand
2. The center of a palindrome must be between two nucleotides (no central base)
3. For odd lengths, a center base would need to complement itself (impossible)

Example: Length 5 → positions 1,2,3,4,5 → center = position 3
- Position 3 must complement position 3 → A↔T, G↔C → impossible for single base

---

## Test Categories

### MUST Tests (Required for DoD)

All MUST tests are justified by evidence or explicitly marked.

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| M1 | EcoRI_RecognitionSite_Detected | GAATTC is a canonical 6bp palindrome | Wikipedia - Restriction enzyme |
| M2 | HindIII_RecognitionSite_Detected | AAGCTT is a canonical 6bp palindrome | Wikipedia - Restriction enzyme |
| M3 | BamHI_RecognitionSite_Detected | GGATCC is a canonical 6bp palindrome | Wikipedia - Restriction enzyme |
| M4 | FourBasePalindrome_GCGC_Detected | 4bp palindrome (GCGC = revcomp of GCGC) | Rosalind REVP |
| M5 | FourBasePalindrome_ATAT_Detected | 4bp palindrome (ATAT = revcomp of ATAT) | Standard test case |
| M6 | ReverseComplementInvariant_Holds | Sequence must equal its reverse complement | Wikipedia - Palindromic sequence |
| M7 | NoPalindromes_ReturnsEmpty | Sequence without palindromic regions | Standard edge case |
| M8 | EmptySequence_ReturnsEmpty | Boundary - empty input | Standard boundary |
| M9 | NullSequence_ThrowsArgumentNullException | Parameter validation | Implementation contract |
| M10 | OddMinLength_ThrowsException | minLength must be even | Implementation constraint |
| M11 | MinLengthLessThan4_ThrowsException | minLength must be ≥ 4 | Implementation constraint |
| M12 | MaxLengthLessThanMinLength_ThrowsException | Invalid parameter combination | Implementation contract |
| M13 | MultiplePalindromes_FindsAll | Multiple palindromes in sequence | Rosalind REVP |
| M14 | Position_IsZeroBased | Positions are 0-indexed | Implementation contract |
| M15 | Length_MatchesSequenceLength | Length property equals Sequence.Length | Invariant |
| M16 | MinLength_RespectsThreshold | Only palindromes ≥ minLength returned | Algorithm specification |
| M17 | MaxLength_RespectsThreshold | Only palindromes ≤ maxLength returned | Algorithm specification |

### SHOULD Tests (Important but not blocking)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| S1 | StringOverload_MatchesDnaSequenceOverload | API consistency | Implementation contract |
| S2 | CaseInsensitivity_HandledCorrectly | Lowercase input processed correctly | Implementation robustness |
| S3 | Rosalind_SampleDataset_CorrectOutput | Validate against published test data | Rosalind REVP |
| S4 | EightBasePalindrome_NotI_Detected | NotI site: GCGGCCGC (8bp) | Wikipedia - Restriction enzyme |
| S5 | OverlappingPalindromes_BothDetected | Overlapping palindromes at different lengths | Edge case |
| S6 | SequenceIsPalindrome_EntireSequenceReturned | Entire sequence is palindrome | Edge case |

### COULD Tests (Nice to have)

| ID | Test Name | Rationale | Evidence |
|----|-----------|-----------|----------|
| C1 | LongPalindrome_12bp_Detected | Maximum default length | Algorithm range |
| C2 | GenomicAnalyzer_SmokeTest | Alternate implementation consistency | Implementation comparison |

---

## Test Audit

### Existing Tests (Pre-Consolidation)

| File | Test Name | Classification | Action |
|------|-----------|----------------|--------|
| RepeatFinderTests.cs | FindPalindromes_EcoRISite_FindsPalindrome | Covered (M1) | Move, enhance |
| RepeatFinderTests.cs | FindPalindromes_HindIIISite_FindsPalindrome | Covered (M2) | Move, enhance |
| RepeatFinderTests.cs | FindPalindromes_ShortPalindrome_FindsFourBasePalindrome | Covered (M4) | Move, enhance |
| RepeatFinderTests.cs | FindPalindromes_MultiplePalindromes_FindsAll | Weak (M13) | Strengthen, move |
| RepeatFinderTests.cs | FindPalindromes_NoPalindromes_ReturnsEmpty | Covered (M7) | Move |
| RepeatFinderTests.cs | FindPalindromes_StringOverload_Works | Covered (S1) | Move |
| RepeatFinderTests.cs | FindPalindromes_EmptySequence_ReturnsEmpty | Covered (M8) | Move |
| RepeatFinderTests.cs | FindPalindromes_OddMinLength_ThrowsException | Covered (M10) | Move |
| GenomicAnalyzerTests.cs | FindPalindromes_EcoRI_FindsIt | Duplicate (M1) | Keep as smoke test, reduce |
| GenomicAnalyzerTests.cs | FindPalindromes_MultipleSites_FindsAll | Duplicate (M13) | Keep as smoke test, reduce |

### Consolidation Plan

1. **Canonical file:** `RepeatFinder_Palindrome_Tests.cs` (new dedicated file)
2. **Remove from RepeatFinderTests.cs:** All palindrome tests (move to dedicated file)
3. **GenomicAnalyzerTests.cs:** Keep only minimal smoke tests (1-2), update comments
4. **Add missing tests:** M3, M5, M6, M9, M11, M12, M14, M15, M16, M17, S2, S3, S4, S5, S6

---

## Test Data from Sources

### Rosalind REVP Sample Dataset

Input:
```
>Rosalind_24
TCAATGCATGCGGGTCTATATGCAT
```

Expected palindromes (position length):
- Position 4, Length 6: ATGCAT
- Position 5, Length 4: TGCA
- Position 6, Length 6: GCATGC
- Position 7, Length 4: CATG
- Position 17, Length 4: ATAT
- Position 18, Length 4: TATA
- Position 20, Length 6: ATGCAT
- Position 21, Length 4: TGCA

Note: Rosalind uses 1-based positions; our implementation uses 0-based.

### Known Restriction Enzyme Recognition Sites

| Enzyme | Sequence | Length | Type |
|--------|----------|--------|------|
| EcoRI | GAATTC | 6 | Sticky |
| BamHI | GGATCC | 6 | Sticky |
| HindIII | AAGCTT | 6 | Sticky |
| TaqI | TCGA | 4 | Sticky |
| AluI | AGCT | 4 | Blunt |
| NotI | GCGGCCGC | 8 | Sticky |
| SmaI | CCCGGG | 6 | Blunt |
| PvuII | CAGCTG | 6 | Blunt |
| EcoRV | GATATC | 6 | Blunt |

---

## Open Questions / Decisions

| Question | Decision | Justification |
|----------|----------|---------------|
| Why minLength ≥ 4? | 2bp palindromes too common | GAATTC, ATAT, etc. are biological min |
| Case sensitivity? | Case-insensitive | Implementation uses ToUpperInvariant() |
| Report at all lengths or only max? | All lengths | Implementation iterates all lengths |
| Overlapping palindromes? | Yes, all reported | Same position can have 4bp and 6bp |

---

## Assumptions

| ID | Assumption | Impact |
|----|------------|--------|
| A1 | minLength = 4 is biologically reasonable minimum for restriction sites | Test M4, M5 use 4bp palindromes |
| A2 | All palindrome lengths (stepping by 2) in range are reported | Tests expect multiple lengths at same position |

---

## Definition of Done Checklist

- [x] All MUST tests implemented
- [x] Tests pass deterministically
- [x] No duplicate tests across files
- [x] Edge cases covered (empty, boundary, invalid)
- [x] Invariants verified with Assert.Multiple
- [x] Clean Code principles applied
