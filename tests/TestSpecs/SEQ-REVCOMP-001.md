# Test Specification: SEQ-REVCOMP-001 - Reverse Complement

**Test Unit ID:** SEQ-REVCOMP-001
**Area:** Composition
**Status:** Complete
**Created:** 2026-01-22
**Last Updated:** 2026-01-22
**Owner:** Algorithm QA Architect

---

## 1. Test Unit Definition

### Canonical Methods
| Method | Class | Type |
|--------|-------|------|
| `TryGetReverseComplement(ReadOnlySpan<char>, Span<char>)` | SequenceExtensions | Canonical (Span API) |

### Delegate/Wrapper Methods
| Method | Class | Type |
|--------|-------|------|
| `ReverseComplement()` | DnaSequence | Instance (creates new object) |
| `ReverseComplement()` | RnaSequence | Instance (creates new object) |
| `GetReverseComplementString(string)` | DnaSequence | Static helper |
| `TryWriteReverseComplement(Span<char>)` | DnaSequence | Instance (delegates to TryGetReverseComplement) |
| `TryGetReverseComplement(ReadOnlySpan, Span)` | DnaSequence | Static (delegates to SequenceExtensions) |

### Invariants
1. **Involution Property:** `ReverseComplement(ReverseComplement(x)) = x` for all valid sequences
2. **Watson-Crick Base Pairing + Reversal:** Complement each base, then reverse (or equivalently, reverse then complement each base)
3. **Length Preservation:** Output length always equals input length
4. **Palindrome Detection:** A sequence is a biological palindrome if `ReverseComplement(x) = x`
5. **Case Insensitivity:** Input can be any case; output is uppercase

### Complexity
- **Time:** O(n)
- **Space:** O(1) for Span API (output written to caller-provided buffer)

---

## 2. Evidence

### Primary Sources

#### Source 1: Wikipedia - Complementarity (molecular biology)
**URL:** https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)
**Accessed:** 2026-01-22

**Key Facts:**
- DNA strands are antiparallel (5'→3' opposite to 3'→5')
- The reverse complement is the sequence on the opposite strand read in the 5'→3' direction
- Base pairing rules: A ↔ T (or U for RNA), G ↔ C
- Complementary sequence to TTAC is GTAA (showing reversal: complement of TTAC is AATG, reversed is GTAA)
- "One sequence can be complementary to another sequence, meaning that they have the base on each position in the complementary (i.e., A to T, C to G) and in the reverse order"

#### Source 2: Wikipedia - Nucleic Acid Sequence (IUPAC Notation)
**URL:** https://en.wikipedia.org/wiki/Nucleic_acid_sequence
**Accessed:** 2026-01-22

**IUPAC Ambiguity Code Complements:**
| Symbol | Meaning | Complement |
|--------|---------|------------|
| A | Adenine | T (or U) |
| C | Cytosine | G |
| G | Guanine | C |
| T | Thymine | A |
| U | Uracil | A |
| R | Purine (A or G) | Y |
| Y | Pyrimidine (C or T) | R |
| N | Any nucleotide | N |

#### Source 3: Biopython Bio.Seq Module
**URL:** https://biopython.org/docs/1.75/api/Bio.Seq.html
**Accessed:** 2026-01-22

**Key Implementation Details:**
- `reverse_complement()` method returns a new Seq object with reversed complemented bases
- Example: `Seq("CCCCCGATAGNR").reverse_complement()` → `Seq('YNCTATCGGGGG')`
- R (purine A or G) complements to Y (pyrimidine C or T)
- Mixed case is supported: `Seq("CCCCCgatA-G").reverse_complement()` → `Seq('C-TatcGGGGG')`
- Gaps and unknown characters are preserved (but reversed in position)
- Protein sequences raise `ValueError: Proteins do not have complements!`

**Standalone function:**
```python
>>> reverse_complement("ACTG-NH")
'DN-CAGT'
```

---

## 3. Test Cases

### 3.1 Must Tests (Required for DoD)

#### MUST-01: Basic Reverse Complement
**Evidence:** Wikipedia Complementarity
**Test:** Verify ACGT → ACGT (palindrome), AACG → CGTT

#### MUST-02: Involution Property
**Evidence:** Mathematical property of reverse complement operation
**Test:** `ReverseComplement(ReverseComplement(x)) = x` for various sequences

#### MUST-03: Empty Sequence Handling
**Evidence:** Edge case
**Test:** Empty input returns true with empty output for TryGetReverseComplement

#### MUST-04: Single Nucleotide
**Evidence:** Edge case, Watson-Crick rules
**Test:** Single base produces its complement (reversed single = same position)

#### MUST-05: Destination Too Small
**Evidence:** API contract for TryGetReverseComplement
**Test:** Returns false when destination.Length < source.Length

#### MUST-06: Biological Palindrome Detection
**Evidence:** Wikipedia Complementarity (restriction enzyme sites)
**Test:** EcoRI site GAATTC is its own reverse complement

#### MUST-07: Case Insensitivity
**Evidence:** Consistent with GetComplementBase behavior
**Test:** Lowercase input produces uppercase reverse complement

#### MUST-08: RNA Uracil Support
**Evidence:** Wikipedia, Biopython
**Test:** ACGU → ACGU (RNA palindrome), RNA sequences work correctly

### 3.2 Should Tests (Recommended)

#### SHOULD-01: Asymmetric Sequences
**Test:** Various non-palindromic sequences produce correct reverse complements

#### SHOULD-02: Long Sequences
**Test:** Verify correctness for longer sequences (100+ bases)

#### SHOULD-03: Destination Exactly Equal Size
**Test:** TryGetReverseComplement succeeds when destination.Length == source.Length

#### SHOULD-04: Destination Larger Than Source
**Test:** TryGetReverseComplement writes only source.Length characters

#### SHOULD-05: Unknown Base Handling
**Test:** Unknown bases (N, X, -) are complemented per GetComplementBase rules and reversed

### 3.3 Could Tests (Optional)

#### COULD-01: IUPAC Ambiguity Codes
**Test:** R → Y, Y → R, etc. (if implemented)
**Note:** Current implementation uses GetComplementBase which returns unknown unchanged

#### COULD-02: Performance with Large Sequences
**Test:** Sequences > 100,000 bases complete in reasonable time

---

## 4. Audit of Existing Tests

### DnaSequenceTests.cs (Lines 81-102)
| Test | Coverage | Status |
|------|----------|--------|
| `ReverseComplement_ReturnsCorrectReverseComplement` | "ACGT"→"ACGT" palindrome | **Keep as smoke** |
| `ReverseComplement_AsymmetricSequence_Works` | "AACG"→"CGTT" | **Keep as smoke** |
| `ReverseComplement_EcoRI_Site` | "GAATTC" palindrome | **Keep as smoke** |

**Assessment:** Good smoke tests for DnaSequence wrapper. Keep all 3 as delegation verification.

### RnaSequenceTests.cs (Lines 88-102)
| Test | Coverage | Status |
|------|----------|--------|
| `ReverseComplement_ReturnsCorrectReverseComplement` | "ACGU"→"ACGU" palindrome | **Keep as smoke** |
| `ReverseComplement_AsymmetricSequence_Works` | "AACG"→"CGUU" | **Keep as smoke** |

**Assessment:** Good smoke tests for RnaSequence wrapper. Keep both.

### PerformanceExtensionsTests.cs (Lines 63-74)
| Test | Coverage | Status |
|------|----------|--------|
| `TryGetReverseComplement_SpanApi_SmokeTest` | "ACGT"→"ACGT" | **Keep as smoke** |

**Assessment:** Single smoke test. Keep in place with updated comment.

### PerformanceExtensionsTests.cs (Lines 158-167)
| Test | Coverage | Status |
|------|----------|--------|
| `DnaSequence_TryWriteReverseComplement_Works` | "AAAA"→"TTTT" | **Keep as smoke** |

**Assessment:** Tests DnaSequence.TryWriteReverseComplement. Keep.

### Coverage Summary
| Category | Status |
|----------|--------|
| TryGetReverseComplement canonical | ✅ Covered |
| Empty sequence | ✅ Covered |
| Single nucleotide | ✅ Covered |
| Destination too small | ✅ Covered |
| Case insensitivity | ✅ Covered |
| Involution property | ✅ Covered |
| Unknown base handling | ✅ Covered |
| Long sequences | ✅ Covered |

---

## 5. Consolidation Plan

### Canonical Test File (NEW)
**File:** `SequenceExtensions_ReverseComplement_Tests.cs`
- All MUST and SHOULD tests for TryGetReverseComplement
- Organized into regions by functionality
- Deep, evidence-based tests

### Wrapper Tests (EXISTING - no changes needed)
**DnaSequenceTests.cs:**
- Keep all 3 ReverseComplement tests as smoke delegation tests (already well-structured)

**RnaSequenceTests.cs:**
- Keep both ReverseComplement tests as smoke delegation tests

**PerformanceExtensionsTests.cs:**
- Keep `TryGetReverseComplement_SpanApi_SmokeTest` (line 63-74)
- Keep `DnaSequence_TryWriteReverseComplement_Works` (line 158-167)

### Tests to Remove
None - existing smoke tests are appropriate for wrapper verification.

---

## 6. Open Questions / Decisions

### Q1: How should unknown characters be handled in reverse complement?
**Decision:** Per implementation, unknown characters are complemented using GetComplementBase (which returns them unchanged) and then their position is reversed. This matches Biopython behavior for gaps.

### Q2: Should RNA and DNA be in the same test file?
**Decision:** The canonical TryGetReverseComplement in SequenceExtensions handles both DNA and RNA (via GetComplementBase which supports U). Tests will cover both.

---

## 7. ASSUMPTIONS

1. **ASSUMPTION:** Output is always uppercase regardless of input case (consistent with GetComplementBase behavior)
2. **ASSUMPTION:** Empty sequence is valid input (returns true with empty output)
3. **ASSUMPTION:** Unknown bases are complemented per GetComplementBase then position-reversed (N stays N, position changes)
