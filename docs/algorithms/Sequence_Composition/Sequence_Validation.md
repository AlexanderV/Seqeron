# Sequence Validation

**Algorithm Group:** Sequence Composition
**Document Version:** 1.0
**Last Updated:** 2026-01-22

---

## 1. Overview

Sequence validation determines whether a nucleic acid string contains only valid nucleotide characters according to IUPAC nomenclature standards.

**Complexity:** O(n) where n = sequence length

---

## 2. IUPAC Standard (Authoritative)

### 2.1 Source

- **IUPAC-IUB Commission on Biochemical Nomenclature (1970)**
  "Abbreviations and symbols for nucleic acids, polynucleotides, and their constituents"
  *Biochemistry* 9(20): 4022–4027. doi:10.1021/bi00822a023

- **NC-IUB (1984)**
  "Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences"
  *Nucleic Acids Research* 13(9): 3021–3030. doi:10.1093/nar/13.9.3021

### 2.2 Standard Nucleotide Codes

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

### 2.3 IUPAC Ambiguity Codes

The IUPAC standard defines additional codes for incompletely specified bases:

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
| - | Gap | — |

---

## 3. Implementation

### 3.1 Validation Mode

This library implements **strict validation**:
- DNA: Only `{A, C, G, T}` are valid (case-insensitive)
- RNA: Only `{A, C, G, U}` are valid (case-insensitive)
- IUPAC ambiguity codes are **not** accepted in strict mode

**Rationale:** Strict validation is preferred for genomic analysis where sequence accuracy is critical. Ambiguous bases should be explicitly handled rather than silently accepted.

### 3.2 Methods

```csharp
// Canonical DNA validation
public static bool IsValidDna(this ReadOnlySpan<char> sequence)

// Canonical RNA validation
public static bool IsValidRna(this ReadOnlySpan<char> sequence)

// Factory pattern with validation
public static bool TryCreate(string sequence, out DnaSequence? result)
```

### 3.3 Behavior

| Input | IsValidDna | IsValidRna |
|-------|------------|------------|
| `""` (empty) | `true` | `true` |
| `"ACGT"` | `true` | `false` |
| `"ACGU"` | `false` | `true` |
| `"acgt"` | `true` | `false` |
| `"ACGN"` | `false` | `false` |
| `"AC GT"` | `false` | `false` |

### 3.4 Empty Sequence Handling

**Decision:** Empty sequences return `true`.

**Justification:** Vacuous truth — an empty sequence contains no invalid characters. This follows common library conventions (e.g., an empty string matches any "all characters satisfy" predicate).

---

## 4. Algorithm

```
function IsValidDna(sequence):
    for each character c in sequence:
        upper_c = ToUppercase(c)
        if upper_c not in {'A', 'C', 'G', 'T'}:
            return false
    return true

function IsValidRna(sequence):
    for each character c in sequence:
        upper_c = ToUppercase(c)
        if upper_c not in {'A', 'C', 'G', 'U'}:
            return false
    return true
```

**Time Complexity:** O(n)
**Space Complexity:** O(1)

---

## 5. Deviations from IUPAC Standard

| Aspect | IUPAC Standard | Implementation | Reason |
|--------|----------------|----------------|--------|
| Ambiguity codes | Defined | Not accepted | Strict validation mode |
| Gap character (-) | Defined | Not accepted | Strict validation mode |
| Case | Not specified | Case-insensitive | Common bioinformatics practice |

---

## 6. References

1. IUPAC-IUB Commission on Biochemical Nomenclature (1970). "Abbreviations and symbols for nucleic acids, polynucleotides, and their constituents." *Biochemistry* 9(20): 4022–4027.

2. NC-IUB (1984). "Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences." *Nucleic Acids Research* 13(9): 3021–3030.

3. Wikipedia: Nucleic acid notation. https://en.wikipedia.org/wiki/Nucleic_acid_notation

4. Bioinformatics.org: IUPAC codes. https://www.bioinformatics.org/sms/iupac.html
