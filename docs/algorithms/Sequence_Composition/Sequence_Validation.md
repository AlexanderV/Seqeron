# Sequence Validation

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition |
| Test Unit ID | SEQ-VALID-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Document Version | 1.0 |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Sequence validation determines whether a nucleic acid string contains only valid nucleotide characters according to IUPAC nomenclature standards. In this repository, the documented validation surface is a strict character-membership check for DNA and RNA alphabets, with linear time complexity in sequence length.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The authoritative basis for nucleotide validation is the IUPAC/IUB nomenclature for nucleic acids and incompletely specified bases.

- **IUPAC-IUB Commission on Biochemical Nomenclature (1970)**
  "Abbreviations and symbols for nucleic acids, polynucleotides, and their constituents"
  *Biochemistry* 9(20): 4022-4027. doi:10.1021/bi00822a023
- **NC-IUB (1984)**
  "Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences"
  *Nucleic Acids Research* 13(9): 3021-3030. doi:10.1093/nar/13.9.3021

**DNA (Deoxyribonucleic acid):**

| Symbol | Nucleobase |
|--------|------------|
| A | Adenine |
| C | Cytosine |
| G | Guanine |
| T | Thymine |

**RNA (Ribonucleic acid):**

| Symbol | Nucleobase |
|--------|------------|
| A | Adenine |
| C | Cytosine |
| G | Guanine |
| U | Uracil |

The IUPAC standard also defines ambiguity codes for incompletely specified bases:

| Symbol | Meaning | Represents |
|--------|---------|------------|
| R | Purine | A or G |
| Y | Pyrimidine | C or T/U |
| S | Strong | G or C |
| W | Weak | A or T/U |
| K | Keto | G or T/U |
| M | Amino | A or C |
| B | Not A | C, G, or T/U |
| D | Not C | A, G, or T/U |
| H | Not G | A, C, or T/U |
| V | Not T/U | A, C, or G |
| N | Any | A, C, G, or T/U |
| - | Gap | - |

### 2.2 Core Model

Sequence validation is a per-character membership test against an allowed alphabet. For unambiguous nucleic acid sequences, the canonical alphabets are `{A, C, G, T}` for DNA and `{A, C, G, U}` for RNA. Any validation rule derived from the IUPAC standard therefore reduces to checking whether each symbol belongs to the selected valid set.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | DNA validation accepts only symbols from `{A, C, G, T}` in the unambiguous alphabet | IUPAC standard nucleotide codes |
| INV-02 | RNA validation accepts only symbols from `{A, C, G, U}` in the unambiguous alphabet | IUPAC standard nucleotide codes |
| INV-03 | Validation fails as soon as a symbol is outside the allowed set | The algorithm is a character-membership scan |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[IsValidDna/IsValidRna] sequence` | `ReadOnlySpan<char>` | required | Sequence to validate | Checked character-by-character after uppercase normalization |
| `[TryCreate] sequence` | `string` | required | DNA sequence to validate and materialize as `DnaSequence` | Must contain only DNA characters accepted by `DnaSequence` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `[IsValidDna/IsValidRna] isValid` | `bool` | `true` when every character belongs to the accepted alphabet for the selected mode |
| `[TryCreate] return value` | `bool` | `true` when a `DnaSequence` instance is created successfully |
| `[TryCreate] result` | `DnaSequence?` | Created DNA sequence when validation succeeds; `null` when validation fails |

### 3.3 Preconditions and Validation

Validation is case-insensitive because characters are normalized with uppercase comparison. In the documented strict mode, DNA accepts only `A/C/G/T` and RNA accepts only `A/C/G/U`; ambiguity codes and gap characters are rejected. Empty sequences return `true` for `IsValidDna` and `IsValidRna`. `TryCreate` returns `false` and `null` when `DnaSequence` construction raises `ArgumentException` for invalid input.

## 4. Algorithm

### 4.1 High-Level Steps

1. Iterate through the input sequence one character at a time.
2. Convert each character to uppercase.
3. Return `false` immediately when a character is not in the allowed alphabet; otherwise return `true` after the scan completes.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `IsValidDna` / `IsValidRna` | `O(n)` | `O(1)` | `n` = sequence length |
| `TryCreate` | `O(n)` | `O(1)` auxiliary | Delegates validation to `DnaSequence` construction |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceExtensions.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs), [DnaSequence.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs)

- `SequenceExtensions.IsValidDna(ReadOnlySpan<char>)`: Validates DNA characters against `A/C/G/T` using uppercase comparison.
- `SequenceExtensions.IsValidRna(ReadOnlySpan<char>)`: Validates RNA characters against `A/C/G/U` using uppercase comparison.
- `DnaSequence.TryCreate(string, out DnaSequence?)`: Factory-style DNA validation and creation that returns `false` when `DnaSequence` construction throws `ArgumentException`.

### 5.2 Current Behavior

This library implements strict validation:

- DNA: Only `{A, C, G, T}` are valid.
- RNA: Only `{A, C, G, U}` are valid.
- IUPAC ambiguity codes are not accepted in strict mode.

Behavior examples documented for the current implementation:

| Input | IsValidDna | IsValidRna |
|-------|------------|------------|
| `""` (empty) | `true` | `true` |
| `"ACGT"` | `true` | `false` |
| `"ACGU"` | `false` | `true` |
| `"acgt"` | `true` | `false` |
| `"ACGN"` | `false` | `false` |
| `"AC GT"` | `false` | `false` |

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Standard DNA nucleotide validation for `A/C/G/T`.
- Standard RNA nucleotide validation for `A/C/G/U`.
- Character-wise validation as a set-membership test over the selected alphabet.

**Intentionally simplified:**

- Ambiguity codes defined by IUPAC are not accepted; **consequence:** sequences containing symbols such as `N`, `R`, or `Y` are rejected in strict mode.
- The gap character `-` is not accepted; **consequence:** aligned sequences containing gap symbols are rejected in strict mode.

**Not implemented:**

- Full IUPAC ambiguity-code validation mode; **users should rely on:** no current alternative documented in this test unit.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Strict validation vs full IUPAC alphabet | Deviation | Sequences containing ambiguity codes or gaps are rejected even though they are defined by the IUPAC standard | accepted | See deviation aspects below: ambiguity codes, gap character, and case handling |

Deviation aspects documented for the current implementation:

| Aspect | IUPAC Standard | Implementation | Reason |
|--------|----------------|----------------|--------|
| Ambiguity codes | Defined | Not accepted | Strict validation mode |
| Gap character (`-`) | Defined | Not accepted | Strict validation mode |
| Case | Not specified | Case-insensitive | Common bioinformatics practice |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | `IsValidDna` and `IsValidRna` return `true` | Vacuous truth: no invalid characters are present |
| Lowercase input such as `"acgt"` | Validation succeeds for DNA | Characters are compared after uppercase normalization |
| Sequence containing `N` | Validation fails | Ambiguity codes are not accepted in strict mode |
| Sequence containing whitespace such as `"AC GT"` | Validation fails | Whitespace is outside the accepted alphabet |

### 6.2 Limitations

The implementation is limited to strict validation. It does not accept IUPAC ambiguity codes or the gap character, even though both are defined by the standard. This is appropriate for the documented strict-mode behavior but not for workflows that must preserve or validate incompletely specified sequences.

## 8. References

1. IUPAC-IUB Commission on Biochemical Nomenclature. 1970. Abbreviations and symbols for nucleic acids, polynucleotides, and their constituents. Biochemistry 9(20). doi:10.1021/bi00822a023
2. NC-IUB. 1984. Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences. Nucleic Acids Research 13(9). doi:10.1093/nar/13.9.3021
3. Wikipedia contributors. 2026. Nucleic acid notation. Wikipedia. https://en.wikipedia.org/wiki/Nucleic_acid_notation
4. Bioinformatics.org. 2026. IUPAC codes. Bioinformatics.org. https://www.bioinformatics.org/sms/iupac.html
