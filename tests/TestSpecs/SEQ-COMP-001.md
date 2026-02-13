# Test Specification: SEQ-COMP-001 - DNA/RNA Complement

**Test Unit ID:** SEQ-COMP-001
**Area:** Composition
**Status:** Active
**Created:** 2026-01-22
**Last Updated:** 2026-02-15
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
| `Complement()` | DnaSequence | Instance (inline Watson-Crick switch) |
| `Complement()` | RnaSequence | Instance (inline Watson-Crick switch) |

### Invariants
1. **Involution Property:** `Complement(Complement(x)) = x` for all valid bases
2. **Watson-Crick Base Pairing:** A ↔ T, G ↔ C (DNA); A ↔ U, G ↔ C (RNA)
3. **Length Preservation:** Output length always equals input length
4. **Uppercase Output:** Output is always uppercase (DnaSequence/RnaSequence normalize to uppercase on construction)

### Complexity
- **Time:** O(n) for sequence complement
- **Space:** O(1) per base, O(n) for full sequence output

---

## 2. Evidence

### Primary Sources

#### Source 1: Wikipedia — Complementarity (molecular biology)
**URL:** https://en.wikipedia.org/wiki/Complementarity_(molecular_biology)
**Accessed:** 2026-02-14

**Key Facts:**
- DNA base complement: A = T (2 H-bonds), G ≡ C (3 H-bonds)
- RNA base complement: A = U, G ≡ C

| Context | Pairs |
|---------|-------|
| DNA | A ↔ T, G ↔ C |
| RNA | A ↔ U, G ↔ C |

#### Source 2: Wikipedia — Nucleic acid notation (IUPAC)
**URL:** https://en.wikipedia.org/wiki/Nucleic_acid_notation
**Accessed:** 2026-02-14

**IUPAC Complement Table (NC-IUB 1984):**

| Symbol | Meaning | Complement |
|--------|---------|------------|
| A | Adenine | T |
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
| - | Gap | - |

#### Source 3: Biopython Bio.Seq — `complement()` / `complement_rna()`
**URL:** https://biopython.org/docs/latest/api/Bio.Seq.html
**Accessed:** 2026-02-14

**Documented Examples (from Biopython docstrings):**

| Function | Input | Output | Notes |
|----------|-------|--------|-------|
| `complement()` | `"CGA"` | `"GCT"` | Standard DNA |
| `complement()` | `"CGAUT"` | `"GCTAA"` | U treated as T in DNA context |
| `complement()` | `"ACGTUacgtuXYZxyz"` | `"TGCAAtgcaaXRZxrz"` | Case-preserving, IUPAC Y→R |
| `complement_rna()` | `"CGA"` | `"GCU"` | RNA output |
| `complement_rna()` | `"CGAUT"` | `"GCUAA"` | T→A in RNA context |
| `complement_rna()` | `"ACGTUacgtuXYZxyz"` | `"UGCAAugcaaXRZxrz"` | Case-preserving, IUPAC Y→R |

**Key Behaviors:**
- **Case preservation:** lowercase input → lowercase output (Biopython). Our API normalizes to uppercase (consistent with DnaSequence/RnaSequence contract).
- **Unknown bases:** non-IUPAC characters pass through unchanged in both DNA and RNA context.
- **T in RNA context:** T pairs with A (`complement_rna` treats T → A).
- **IUPAC ambiguity codes:** Biopython applies IUPAC complement rules (Y→R, M→K, etc.). Our API does not currently support IUPAC codes (see COULD-02).

### Behavior-to-Source Mapping

| Behavior | Source |
|----------|--------|
| A↔T, G↔C (DNA) | Wikipedia Complementarity; IUPAC table; Biopython |
| A↔U, G↔C (RNA) | Wikipedia Complementarity; IUPAC table; Biopython |
| T→A in RNA context | Biopython `complement_rna("CGAUT")` → `"GCUAA"` |
| Unknown DNA bases pass through | Biopython `complement("XYZ")` — X stays X, Z stays Z |
| Unknown RNA bases pass through | Biopython `complement_rna("XYZ")` — X stays X, Z stays Z |
| Gap characters pass through | IUPAC: gap complement is gap; Biopython preserves gaps |
| Uppercase output | DnaSequence/RnaSequence normalize to uppercase on construction |
| Empty sequence returns empty | Standard API contract (zero-length = trivially valid) |
| Involution: comp(comp(x)) = x | Mathematical property of base pairing bijection |

---

## 3. Test Cases

### 3.1 Must Tests (Required for DoD)

#### MUST-01: Standard Watson-Crick Base Pairing (DNA)
**Evidence:** Wikipedia Complementarity; IUPAC Table
**Test:** Verify A→T, T→A, G→C, C→G for GetComplementBase

#### MUST-02: Case Insensitivity with Uppercase Output
**Evidence:** DnaSequence/RnaSequence normalize to uppercase; API contract
**Test:** Lowercase input (a, t, g, c) returns uppercase complement

#### MUST-03: DNA Complement Supports Uracil (U → A)
**Evidence:** IUPAC Table: U complement = A; Biopython: U treated as T in DNA context
**Test:** GetComplementBase('U') = 'A', GetComplementBase('u') = 'A'

#### MUST-04: Involution Property
**Evidence:** Mathematical property of complement bijection
**Test:** `Complement(Complement(x)) = x` for all standard bases (DNA: ATGC; RNA: AUGC)

#### MUST-05: Unknown Base Handling (DNA)
**Evidence:** Biopython: `complement("XYZ")` — X, Z pass through unchanged
**Test:** Non-IUPAC characters (N, X, -, ., ?, *) return unchanged

#### MUST-06: TryGetComplement — Destination Too Small
**Evidence:** API contract (Try pattern)
**Test:** Returns false when destination.Length < source.Length

#### MUST-07: TryGetComplement — Correct Complement
**Evidence:** Watson-Crick rules; Biopython cross-verification
**Test:** Full sequence complement with sufficient destination buffer

#### MUST-08: Empty Sequence Handling
**Evidence:** Standard API contract
**Test:** Empty input → returns true, no output written

#### MUST-09: RNA Complement (GetRnaComplementBase)
**Evidence:** Wikipedia Complementarity; Biopython `complement_rna()`
**Test:** A→U, U→A, G→C, C→G, T→A (T pairs with A per Biopython)

#### MUST-10: RNA Unknown Base Handling
**Evidence:** Biopython: `complement_rna("XYZ")` — X, Z pass through unchanged
**Test:** Unknown bases return unchanged (not 'N')

### 3.2 Should Tests (Recommended)

#### SHOULD-01: Single Character Sequences
**Test:** TryGetComplement works correctly for single-character sequences

#### SHOULD-02: Mixed Case Full Sequence
**Test:** Verify entire sequence with mixed case produces correct uppercase complement

#### SHOULD-03: Destination Exactly Equal Size
**Test:** TryGetComplement succeeds when destination.Length == source.Length

#### SHOULD-04: Destination Larger Than Source
**Test:** TryGetComplement writes only source.Length characters

#### SHOULD-05: All Same Base Sequences
**Test:** Sequences like "AAAA" → "TTTT", "GGGG" → "CCCC"

### 3.3 Could Tests (Optional)

#### COULD-01: Very Long Sequences
**Test:** Performance/correctness for sequences > 10,000 bases

#### COULD-02: IUPAC Ambiguity Codes
**Test:** W→W, S→S, M→K, K→M, R→Y, Y→R, B→V, D→H, H→D, V→B
**Status:** Not implemented. Current implementation passes unknown characters through unchanged.

---

## 4. Biopython Cross-Verification

| Input | Biopython `complement()` | Our `GetComplementBase` per-char | Match |
|-------|--------------------------|-----------------------------------|-------|
| `"ACGT"` | `"TGCA"` | `"TGCA"` | ✅ |
| `"CGAUT"` | `"GCTAA"` | `"GCTAA"` | ✅ |

| Input | Biopython `complement_rna()` | Our `GetRnaComplementBase` per-char | Match |
|-------|-------------------------------|---------------------------------------|-------|
| `"CGAUT"` | `"GCUAA"` | `"GCUAA"` | ✅ |

**Note on case:** Biopython preserves case (`a`→`t`). Our API normalizes to uppercase (`a`→`T`). This is a deliberate design choice: `DnaSequence` and `RnaSequence` normalize input to uppercase on construction.

---

## 5. Coverage Summary

| Category | Status | Tests |
|----------|--------|-------|
| GetComplementBase Watson-Crick | ✅ Sourced | MUST-01 |
| GetComplementBase case handling | ✅ Sourced | MUST-02 |
| GetComplementBase Uracil | ✅ Sourced | MUST-03 |
| Involution property | ✅ Sourced | MUST-04 |
| DNA unknown bases (pass-through) | ✅ Sourced | MUST-05 |
| TryGetComplement buffer validation | ✅ Sourced | MUST-06 |
| TryGetComplement correctness | ✅ Sourced | MUST-07 |
| Empty sequence | ✅ Sourced | MUST-08 |
| GetRnaComplementBase (incl. T→A) | ✅ Sourced | MUST-09 |
| RNA unknown bases (pass-through) | ✅ Sourced | MUST-10 |
| Biopython cross-verification | ✅ 3 examples | Section 4 |

### Coverage Classification (2026-02-15)

**Canonical file:** `SequenceExtensions_Complement_Tests.cs`

| Metric | Before | After |
|--------|--------|-------|
| Test methods | 45 | 21 |
| Test runs | 45 | 22 |
| Complement filter total | 139 | 116 |

| Action | Count | Details |
|--------|-------|---------|
| 🔁 Duplicate removed | 22 | Individual base tests covered by `Assert.Multiple` combined tests |
| 🔁 Duplicate merged | 2→1 | `AllAdenine` + `AllGuanine` → parametrized `HomogeneousSequence` |
| 🔁 Duplicate merged | 2→1 | `Uracil` + `LowercaseUracil` → single `Assert.Multiple` |
| ⚠ Weak strengthened | 1 | `EmptySource`: added buffer pre-fill + verify destination untouched |
| ❌ Missing | 0 | — |

**Kept unchanged:** 5 FsCheck properties (`SequenceProperties.cs`), 2 wrapper smoke tests (`DnaSequenceTests`, `RnaSequenceTests`)

---

## 6. Deviations from Sources

| # | Behavior | Source Behavior | Our Behavior | Justification |
|---|----------|-----------------|--------------|---------------|
| D1 | Case handling | Biopython preserves case (`a`→`t`) | Always uppercase | DnaSequence/RnaSequence normalize to uppercase. Consistent internal API contract. |
| D2 | IUPAC ambiguity codes | Biopython complements Y→R, M→K, etc. | Passed through unchanged | Not implemented (COULD-02). Standard ACGTU bases fully handled. |
