# RNA-specific Complement

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition |
| Test Unit ID | SEQ-RNACOMP-001 |
| Related Projects | Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Returns the Watson–Crick complement of a single nucleotide **in the RNA alphabet**, i.e. emitting U rather than T. Unlike the DNA complement (`A → T`), the RNA complement maps `A → U`, and a thymine in the input is treated as a uracil (`T → A`). The mapping is IUPAC-complete: it also complements the eleven degenerate ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) [1][2]. The operation is exact and specification-driven (a fixed lookup), O(1) per base.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

RNA is single-stranded but base-pairs with a complementary strand (in transcription, RNA:DNA hybrids, and antisense design). The complementary base of a nucleotide is the one it forms a canonical pair with: A pairs with U (in RNA; A pairs with T in DNA), and C pairs with G [5]. Degenerate IUPAC symbols denote sets of bases; the complement of a set symbol is the symbol denoting the complements of its members [1].

### 2.2 Core Model

The RNA complement is the lookup defined by Biopython's `ambiguous_rna_complement` table [2], identical to the DNA table except the emitted alphabet uses U for T:

`A→U, C→G, G→C, U→A, M→K, R→Y, W→W, S→S, Y→R, K→M, V→B, H→D, D→H, B→V, N→N` [2].

Thymine, if present in an RNA context, is treated as uracil: Biopython builds the table with `ambiguous_rna_complement["T"] = ambiguous_rna_complement["U"]`, so `T → A` [3]. The documentation states this verbatim: "Any T in the sequence is treated as a U" [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Canonical pairing A↔U and C↔G are reciprocal | complement table is symmetric for these pairs [2][5] |
| INV-02 | No recognized base maps to `T`; `A → U` | RNA alphabet emits U, not T [2][3][4] |
| INV-03 | `T` maps to `A` (treated as U) | `ambiguous_rna_complement["T"]=...["U"]="A"` [3][4] |
| INV-04 | W, S, N are self-complementary | their member sets are closed under complement [1][2] |
| INV-05 | Reciprocal ambiguity pairs R↔Y, M↔K, D↔H, B↔V | symmetric entries in the table [2] |
| INV-06 | Complement is an involution on the canonical RNA bases and the eleven ambiguity codes | every recognized symbol's complement is itself recognized and maps back [2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `nucleotide` | `char` | required | A single nucleotide symbol | Any `char`; recognized values are A,C,G,U,T and IUPAC codes R,Y,S,W,K,M,B,D,H,V,N (case-insensitive) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `char` | The RNA complement. Recognized bases/codes are returned **uppercase**; `T`/`t` returns `A`; unrecognized characters are returned unchanged. |

### 3.3 Preconditions and Validation

No exceptions are thrown for any `char` input. Input is case-insensitive; recognized symbols are normalized to uppercase on output (repository convention, mirroring the DNA sibling `GetComplementBase`, SEQ-COMP-001). `T`/`t` are accepted and normalized to the RNA convention (`T → A`, treated as U) [3][4]. Characters outside the IUPAC nucleotide set (gaps `-`/`.`, digits, `Z`, whitespace) pass through unchanged, including their original case.

## 4. Algorithm

### 4.1 High-Level Steps

1. Match the input character (case-insensitive) against the IUPAC RNA complement table.
2. If recognized, return its complement (uppercase), with `T`/`t` treated as `U`.
3. Otherwise return the character unchanged.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

RNA complement reference table (origin: Biopython `ambiguous_rna_complement` [2]; cross-checked with bioinformatics.org IUPAC table [5]; standard: NC-IUB 1984 [1]):

| In | Out | In | Out | In | Out |
|----|-----|----|-----|----|-----|
| A | U | R | Y | B | V |
| U | A | Y | R | D | H |
| C | G | S | S | H | D |
| G | C | W | W | V | B |
| T | A | K | M | N | N |
| | | M | K | | |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GetRnaComplementBase(char)` | O(1) | O(1) | single `switch` expression over a fixed alphabet |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceExtensions.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs)

- `SequenceExtensions.GetRnaComplementBase(char)`: returns the IUPAC-complete RNA complement of a single nucleotide.

### 5.2 Current Behavior

Implemented as a `switch` expression with `[MethodImpl(AggressiveInlining)]`. Recognized bases and ambiguity codes are uppercased on output; the `_ => nucleotide` arm passes unrecognized characters through verbatim (preserving original case). No search/matching is involved, so the repository suffix tree is not applicable (N/A).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The full `ambiguous_rna_complement` table (A→U, C→G, G→C, U→A and all eleven ambiguity codes) [2].
- T treated as U (`T → A`) [3][4].
- Self-complementary W, S, N and reciprocal pairs R↔Y, M↔K, D↔H, B↔V [2].

**Intentionally simplified:**

- (none)

**Not implemented:**

- (none)

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Case normalization to uppercase for recognized bases | Deviation | Biopython preserves input case (`a → u`); this method returns `A → U` uppercase. Only letter casing differs; the complement identity is unchanged. | accepted | Established repo convention (SEQ-COMP-001 MUST-02; sequences normalize to uppercase). Non-IUPAC characters still pass through with original case. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `'T'` / `'t'` | `'A'` | T treated as U in RNA context [3][4] |
| lowercase recognized base (`'a'`) | uppercase complement (`'U'`) | repo uppercase convention (§5.4) |
| gap `'-'` / `'.'`, digit, `'Z'`, space | returned unchanged | not in IUPAC nucleotide set; pass-through [2][5] |
| `'A'` (RNA vs DNA) | `'U'` (vs `'T'` for `GetComplementBase`) | RNA alphabet emits U [4] |

### 6.2 Limitations

Operates on a single character only; whole-sequence RNA complement/reverse-complement is composed by the caller (or via the DNA span helpers for the DNA alphabet). It does not validate that the surrounding sequence is RNA; a stray `T` is silently mapped to `A`.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
char c1 = SequenceExtensions.GetRnaComplementBase('A'); // 'U'
char c2 = SequenceExtensions.GetRnaComplementBase('T'); // 'A' (T treated as U)
char c3 = SequenceExtensions.GetRnaComplementBase('R'); // 'Y'
```

**Numerical / biological walk-through:**

For the Biopython example string `"ACGTUacgtuXYZxyz"`, the per-base forward RNA complement is `"UGCAAugcaaXRZxrz"` in Biopython (case preserved) [4]. Under this repository's uppercase convention, recognized bases are uppercased while non-IUPAC characters (`X`, `Z`, `x`, `z`) pass through verbatim, giving `"UGCAAUGCAAXRZxRz"`.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceExtensions_GetRnaComplementBase_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceExtensions_GetRnaComplementBase_Tests.cs) — covers `INV-01`–`INV-06`
- Evidence: [SEQ-RNACOMP-001-Evidence.md](../../../docs/Evidence/SEQ-RNACOMP-001-Evidence.md)
- Related algorithms: [Sequence_Validation](./Sequence_Validation.md)

## 8. References

1. Cornish-Bowden A. 1985. Nomenclature for incompletely specified bases in nucleic acid sequences: recommendations 1984. *Nucleic Acids Research* 13(9):3021–3030. https://doi.org/10.1093/nar/13.9.3021
2. Biopython contributors. `Bio/Data/IUPACData.py` — `ambiguous_rna_complement`. https://raw.githubusercontent.com/biopython/biopython/master/Bio/Data/IUPACData.py
3. Biopython contributors. `Bio/Seq.py` — `complement_rna` / `_rna_complement_table` (`ambiguous_rna_complement["T"]=...["U"]`). https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py
4. Biopython 1.79. `Bio.Seq` module documentation — `complement_rna`, `reverse_complement_rna` worked examples. https://biopython.org/docs/1.79/api/Bio.Seq.html
5. Bioinformatics.org Sequence Manipulation Suite. IUPAC codes table. https://www.bioinformatics.org/sms/iupac.html
