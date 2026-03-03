# Test Specification: REP-PALIN-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | REP-PALIN-001 |
| **Area** | Repeats |
| **Title** | Palindrome Detection |
| **Status** | ☑ Complete |
| **Created** | 2026-01-22 |
| **Last Updated** | 2026-03-03 |

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
| M10 | OddMinLength_ThrowsException | minLength must be even | Mathematical proof: no odd-length DNA palindrome exists |
| M10b | StringOverload_OddMinLength_ThrowsException | String overload: minLength must be even | API consistency with DnaSequence overload |
| M11 | MinLengthLessThan4_ThrowsException | minLength must be ≥ 4 | Rosalind REVP: "length between 4 and 12" |
| M11b | StringOverload_MinLengthLessThan4_ThrowsException | String overload: minLength must be ≥ 4 | API consistency with DnaSequence overload |
| M12 | MaxLengthLessThanMinLength_ThrowsException | Invalid parameter combination | Implementation contract |
| M12b | StringOverload_MaxLengthLessThanMinLength_ThrowsException | String overload: maxLength ≥ minLength | API consistency with DnaSequence overload |
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
| C1 | LongPalindrome_12bp_Detected | Maximum default length | Rosalind REVP: lengths 4–12 |

---

## Test Audit

### Consolidated Test Files

| File | Scope |
|------|-------|
| `RepeatFinder_Palindrome_Tests.cs` | All MUST, SHOULD, COULD tests for `RepeatFinder.FindPalindromes` (24 tests) |
| `GenomicAnalyzerTests.cs` | Minimal smoke tests (2) for `GenomicAnalyzer.FindPalindromes` |
| `RepeatFinderProperties.cs` | Property-based invariant tests (bounds, reverse complement, even length) |
| `RepeatSnapshotTests.cs` | Snapshot/approval test for palindrome output |

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

## Resolved Decisions

| Question | Decision | Source |
|----------|----------|--------|
| Why minLength ≥ 4? | 4 is the minimum biologically relevant palindrome length | Rosalind REVP: "length between 4 and 12"; Wikipedia: Type II enzymes recognize 4–8 bp; real enzymes TaqI (TCGA), AluI (AGCT), HaeIII (GGCC) are 4bp |
| Why even length only? | DNA palindromes are mathematically impossible at odd lengths | For odd-length seq of length 2k+1, the center base at position k must equal its own complement (A↔T, G↔C). No base is its own complement. |
| Case sensitivity? | Case-insensitive | Both overloads normalize to uppercase via `ToUpperInvariant()` |
| Report at all lengths or only max? | All palindromes at all even lengths in [minLength, maxLength] | Rosalind REVP: "the position and length of **every** reverse palindrome" |
| Overlapping palindromes? | Yes, all reported | Rosalind REVP sample output: position 4 has 6bp, position 5 has 4bp (overlapping) |
| Parameter validation for string overload? | Same validation as DnaSequence overload | Both overloads enforce: minLength ≥ 4, minLength even, maxLength ≥ minLength |

---

## Definition of Done Checklist

- [x] All MUST tests implemented
- [x] Tests pass deterministically
- [x] No duplicate tests across files
- [x] Edge cases covered (empty, boundary, invalid)
- [x] Invariants verified with Assert.Multiple
- [x] Clean Code principles applied
- [x] No remaining assumptions — all decisions grounded in external sources
- [x] Both overloads (DnaSequence, string) have consistent parameter validation
