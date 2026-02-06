# Test Specification: SEQ-COMP-001 - DNA Complement

**Test Unit ID:** SEQ-COMP-001
**Area:** Composition
**Status:** Draft
**Created:** 2026-01-22
**Last Updated:** 2026-01-22
**Owner:** Algorithm QA Architect

---

## 1. Test Unit Definition

### Canonical Methods
| Method | Class | Type |
|--------|-------|------|
| `GetComplementBase(char)` | SequenceExtensions | Canonical (DNA) |
| `GetRnaComplementBase(char)` | SequenceExtensions | Canonical (RNA) |
| `TryGetComplement(ReadOnlySpan<char>, Span<char>)` | SequenceExtensions | Span API |

### Delegate/Wrapper Methods
| Method | Class | Type |
|--------|-------|------|
| `Complement()` | DnaSequence | Instance (delegates to GetComplementBase) |
| `Complement()` | RnaSequence | Instance (delegates to GetRnaComplementBase) |

### Invariants
1. **Involution Property:** `Complement(Complement(x)) = x` for all valid bases
2. **Watson-Crick Base Pairing:** A ↔ T, G ↔ C (DNA); A ↔ U, G ↔ C (RNA)
3. **Length Preservation:** Output length always equals input length
4. **Case Insensitivity:** Input can be any case; output is uppercase

### Complexity
- **Time:** O(n) for sequence complement
- **Space:** O(1) per base, O(n) for full sequence output

---

## 2. Evidence

### Primary Sources

#### Source 1: Wikipedia - Complementarity (molecular biology)
**URL:** https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)
**Accessed:** 2026-01-22

**Key Facts:**
- "Adenine and guanine are purines, while thymine, cytosine and uracil are pyrimidines"
- "The base complement A = T shares two hydrogen bonds, while the base pair G ≡ C has three hydrogen bonds"
- DNA: A ↔ T, G ↔ C
- RNA: A ↔ U, G ↔ C

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
| W | Weak (A or T) | W |
| S | Strong (C or G) | S |
| M | Amino (A or C) | K |
| K | Keto (G or T) | M |
| R | Purine (A or G) | Y |
| Y | Pyrimidine (C or T) | R |
| B | Not A (C, G, T) | V |
| D | Not C (A, G, T) | H |
| H | Not G (A, C, T) | D |
| V | Not T (A, C, G) | B |
| N | Any nucleotide | N |

#### Source 3: Biopython Bio.Seq Module
**URL:** https://biopython.org/docs/1.75/api/Bio.Seq.html
**Accessed:** 2026-01-22

**Key Implementation Details:**
- `complement()` method returns a new Seq object with complemented bases
- Supports mixed case sequences (preserves case pattern in result)
- Ambiguous character D (G, A, T) complements to H (C, T, A)
- Example: `Seq("CCCCCgatA-GD").complement()` → `Seq('GGGGGctaT-CH')`
- Unknown bases (like gaps `-`) are preserved unchanged

**Edge Cases from Biopython:**
- Protein sequences raise `ValueError: Proteins do not have complements!`
- Mixed case is handled: lowercase input → lowercase output (case-preserving)

---

## 3. Test Cases

### 3.1 Must Tests (Required for DoD)

#### MUST-01: Standard Watson-Crick Base Pairing (DNA)
**Evidence:** Wikipedia Complementarity
**Test:** Verify A→T, T→A, G→C, C→G for GetComplementBase

#### MUST-02: Case Insensitivity with Uppercase Output
**Evidence:** Implementation requirement, Biopython behavior
**Test:** Lowercase input (a, t, g, c) should return uppercase complement

#### MUST-03: RNA Uracil Complement (U ↔ A)
**Evidence:** Wikipedia Complementarity
**Test:** GetComplementBase('U') = 'A', GetComplementBase('u') = 'A'

#### MUST-04: Involution Property
**Evidence:** Mathematical property of complement operation
**Test:** `Complement(Complement(x)) = x` for all standard bases

#### MUST-05: Unknown Base Handling
**Evidence:** Implementation behavior (returns unchanged)
**Test:** Non-ACGTU bases (e.g., 'N', 'X', '-') should return unchanged

#### MUST-06: TryGetComplement - Destination Too Small
**Evidence:** API contract
**Test:** Returns false when destination.Length < source.Length

#### MUST-07: TryGetComplement - Correct Complement
**Evidence:** Watson-Crick rules
**Test:** Full sequence complement with sufficient destination buffer

#### MUST-08: Empty Sequence Handling
**Evidence:** Edge case
**Test:** Empty input returns empty output (TryGetComplement returns true)

#### MUST-09: RNA Complement (GetRnaComplementBase)
**Evidence:** Wikipedia Complementarity
**Test:** A→U, U→A, G→C, C→G for RNA

#### MUST-10: RNA Unknown Base → N
**Evidence:** Implementation behavior
**Test:** Unknown bases return 'N' in RNA context

### 3.2 Should Tests (Recommended)

#### SHOULD-01: Single Character Sequences
**Test:** TryGetComplement works correctly for single-character sequences

#### SHOULD-02: Mixed Case Full Sequence
**Test:** Verify entire sequence with mixed case produces correct complement

#### SHOULD-03: Destination Exactly Equal Size
**Test:** TryGetComplement succeeds when destination.Length == source.Length

#### SHOULD-04: Destination Larger Than Source
**Test:** TryGetComplement writes only source.Length characters

#### SHOULD-05: All Same Base Sequences
**Test:** Sequences like "AAAA" → "TTTT", "GGGG" → "CCCC"

### 3.3 Could Tests (Optional)

#### COULD-01: Very Long Sequences
**Test:** Performance/correctness for sequences > 10,000 bases

#### COULD-02: IUPAC Ambiguity Codes (Extended)
**Test:** W→W, S→S, M→K, K→M, R→Y, Y→R, etc.
**Note:** Current implementation returns unknown bases unchanged

---

## 4. Audit of Existing Tests

### DnaSequenceTests.cs (Lines 63-80)
| Test | Coverage | Status |
|------|----------|--------|
| `Complement_ReturnsCorrectComplement` | Basic "ACGT"→"TGCA" | Weak (wrapper only, no edge cases) |
| `Complement_LongerSequence_ReturnsCorrectComplement` | "AATTCCGG"→"TTAAGGCC" | Weak (still basic path) |

**Assessment:** Tests are smoke-level for wrapper. Keep 1-2 as delegation verification.

### PerformanceExtensionsTests.cs (Lines 46-65)
| Test | Coverage | Status |
|------|----------|--------|
| `TryGetComplement_ReturnsCorrectComplement` | "ACGT"→"TGCA" | Covered but minimal |
| `TryGetComplement_DestinationTooSmall_ReturnsFalse` | Destination too small | Covered |

**Assessment:** Comprehensive coverage in canonical file `SequenceExtensions_Complement_Tests.cs`.

### RnaSequenceTests.cs (Lines 70-105)
| Test | Coverage | Status |
|------|----------|--------|
| `Complement_ReturnsCorrectComplement` | "ACGU"→"UGCA" | Wrapper smoke |
| `Complement_LongerSequence` | Longer RNA | Wrapper smoke |
| `ReverseComplement_*` | Multiple | Different Test Unit (SEQ-REVCOMP-001) |

**Assessment:** Basic RNA wrapper tests exist.

### Coverage Summary
| Category | Status |
|----------|--------|
| GetComplementBase canonical | ✅ Covered |
| GetRnaComplementBase | ✅ Covered |
| TryGetComplement edge cases | ✅ Covered |
| Case insensitivity | ✅ Covered |
| Involution property | ✅ Covered |
| Unknown base handling | ✅ Covered |
| Empty sequence | ✅ Covered |

---

## 5. Consolidation Plan

### Canonical Test File (NEW)
**File:** `SequenceExtensions_Complement_Tests.cs`
- All MUST and SHOULD tests for GetComplementBase, GetRnaComplementBase, TryGetComplement
- Organized into regions by method

### Wrapper Tests (EXISTING - minimal changes)
**DnaSequenceTests.cs:**
- Keep `Complement_ReturnsCorrectComplement` as smoke delegation test
- Remove `Complement_LongerSequence_ReturnsCorrectComplement` (duplicate of canonical)

**RnaSequenceTests.cs:**
- Keep `Complement_ReturnsCorrectComplement` as smoke delegation test
- Remove `Complement_LongerSequence_ReturnsCorrectComplement` (duplicate)

**PerformanceExtensionsTests.cs:**
- Move TryGetComplement tests to canonical file
- Keep file focused on cancellation/async tests

---

## 6. Open Questions / Decisions

### Q1: Should IUPAC ambiguity codes be supported?
**Decision:** Current implementation returns unknown bases unchanged for DNA context and 'N' for RNA context. This is acceptable per evidence. IUPAC support could be added later (COULD test).

### Q2: What happens with gap characters ('-')?
**Decision:** Per Biopython evidence, gaps are preserved unchanged. Current implementation does this.

---

## 7. ASSUMPTIONS

1. **ASSUMPTION:** Output is always uppercase regardless of input case (implementation choice, consistent with DnaSequence normalization)
2. **ASSUMPTION:** Unknown bases in DNA context return unchanged (implementation behavior)
3. **ASSUMPTION:** Empty sequence is valid input (returns true with empty output)
