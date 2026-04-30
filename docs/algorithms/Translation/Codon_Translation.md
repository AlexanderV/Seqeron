# Codon Translation

| Field | Value |
|-------|-------|
| Algorithm Group | Translation |
| Test Unit ID | TRANS-CODON-001 |
| Related Projects | N/A |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Codon translation maps a nucleotide triplet to a single-letter amino-acid code according to a selected genetic code table.[1][4] In this repository, the `GeneticCode` class provides codon translation, start/stop-codon classification, reverse lookup from amino acid to codons, and factory access to the supported NCBI tables. The implementation accepts DNA or RNA codons, is case-insensitive, and distinguishes between ambiguous IUPAC codons, which return `X`, and invalid non-IUPAC codons, which throw. Supported tables are Standard (1), Vertebrate Mitochondrial (2), Yeast Mitochondrial (3), and Bacterial/Plastid (11).[4][5]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The genetic code is the rule set that translates codons into amino acids. The original document records these core properties:[1][2][3]

| Property | Meaning |
|----------|---------|
| Triplet code | Each codon contains exactly 3 nucleotides. |
| Non-overlapping | Codons are read sequentially without overlap. |
| Degenerate | 64 codons encode 20 amino acids plus stop signals. |
| Nearly universal | Most organisms use the standard code with limited alternative tables. |

The original file also recorded the supported alternative tables:[4]

| Table | Name | Key AA Differences | Start Codons |
|-------|------|--------------------|--------------|
| 1 | Standard | Universal default | `AUG`, `UUG`, `CUG` |
| 2 | Vertebrate Mitochondrial | `AGA/AGG = *`, `AUA = M`, `UGA = W` | `AUG`, `AUA`, `AUU`, `AUC`, `GUG` |
| 3 | Yeast Mitochondrial | `CUU/CUC/CUA/CUG = T`, `AUA = M`, `UGA = W` | `AUG`, `AUA`, `GUG` |
| 11 | Bacterial/Plastid | Same amino-acid mapping as standard | `AUG`, `GUG`, `UUG`, `CUG`, `AUU`, `AUC`, `AUA` |

### 2.2 Core Model

Codon translation is a constant-time dictionary lookup after normalization of the input codon to uppercase RNA notation. Stop codons return `'*'`. The repository also exposes start-codon and stop-codon membership checks through explicit sets.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The caller selects a genetic code table appropriate for the organism or compartment being modeled. | Translation results can be correct for the chosen table but wrong biologically. |
| ASM-02 | The input to `Translate` is exactly one codon. | Longer or shorter inputs are rejected because codon translation is defined on triplets only. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every codon present in a table maps to exactly one amino acid character. | Each `GeneticCode` stores a dictionary from codon to amino acid. |
| INV-02 | Standard, Vertebrate Mitochondrial, Yeast Mitochondrial, and Bacterial/Plastid each define 64 codons. | The built-in tables are created from the standard base table with table-specific overrides. |
| INV-03 | Stop codons translate to `'*'`. | The built-in codon tables store `'*'` for stop codons. |
| INV-04 | `ATG` and `AUG` translate identically. | `Translate` normalizes `T` to `U` before lookup. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[Translate] codon` | `string` | required | Single DNA or RNA codon. | Must be exactly 3 characters or `Translate` throws `ArgumentException`. |
| `[IsStartCodon] codon` | `string` | required | Codon to classify as a start codon. | Invalid length or `null` returns `false`. |
| `[IsStopCodon] codon` | `string` | required | Codon to classify as a stop codon. | Invalid length or `null` returns `false`. |
| `[GetCodonsForAminoAcid] aminoAcid` | `char` | required | Amino-acid code for reverse lookup. | Lookup is case-insensitive. |
| `[GetByTableNumber] tableNumber` | `int` | required | NCBI table number for a supported code. | Only `1`, `2`, `3`, and `11` are supported. |

### 3.2 Output / Return Value

| Name | Type | Description |
|------|------|-------------|
| `Translate` result | `char` | Single-letter amino acid, `'*'` for stop, or `X` for ambiguous IUPAC codons. |
| `IsStartCodon` result | `bool` | Whether the normalized codon is in the table's start-codon set. |
| `IsStopCodon` result | `bool` | Whether the normalized codon is in the table's stop-codon set. |
| `GetCodonsForAminoAcid` result | `IEnumerable<string>` | All codons in the table that encode the supplied amino acid. |
| `GetByTableNumber` result | `GeneticCode` | The supported genetic code instance for the requested table number. |

### 3.3 Preconditions and Validation

`Translate` throws `ArgumentException` when the input is `null`, empty, or not exactly three characters long. After normalization, a codon found in the table is translated directly. If the codon is not present in the table but consists only of valid IUPAC nucleotide symbols, the method returns `X`. Non-IUPAC codons throw `ArgumentException`.[5]

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate that the input codon is non-empty and exactly three characters.
2. Normalize the codon to uppercase RNA notation by replacing `T` with `U`.
3. If the normalized codon exists in the active table, return its amino acid.
4. Otherwise, if all characters are valid IUPAC nucleotide symbols, return `X`.
5. Otherwise, throw `ArgumentException`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The built-in table support is fixed to the four NCBI tables listed in Section 2.1. Degeneracy observed in the original document is preserved below:[4]

| Number of codons | Amino acids |
|------------------|-------------|
| 1 | Met (`M`), Trp (`W`) |
| 2 | Phe, Tyr, His, Gln, Asn, Lys, Asp, Glu, Cys |
| 3 | Ile |
| 4 | Val, Pro, Thr, Ala, Gly |
| 6 | Leu, Ser, Arg |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `Translate` | `O(1)` | `O(1)` | Constant-time normalization and dictionary lookup per codon. |
| `IsStartCodon` / `IsStopCodon` | `O(1)` | `O(1)` | Set membership after normalization. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GeneticCode.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/GeneticCode.cs)

- `GeneticCode.Translate(string)`
- `GeneticCode.IsStartCodon(string)`
- `GeneticCode.IsStopCodon(string)`
- `GeneticCode.GetCodonsForAminoAcid(char)`
- `GeneticCode.GetByTableNumber(int)`

### 5.2 Current Behavior

The repository supports tables `1`, `2`, `3`, and `11` only. `Translate` normalizes `T` to `U`, looks up the codon in the active table, and returns `X` for ambiguous IUPAC codons such as `ANN` or `NNN`. Inputs containing non-IUPAC symbols such as `XYZ` or `12G` throw `ArgumentException`.[5] `IsStartCodon` and `IsStopCodon` return `false` rather than throwing when the input is `null`, empty, or not exactly three characters long.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The built-in tables match NCBI Tables `1`, `2`, `3`, and `11` for codon translations and documented start/stop sets.[4][5]
- Stop codons are represented as `'*'` and DNA codons are treated equivalently to RNA codons through `T -> U` normalization.[3][5]

**Intentionally simplified:**

- Ambiguous IUPAC codons are collapsed to `X` rather than expanded across all possible concrete codons; **consequence:** ambiguity is preserved as unknown amino acid instead of resolved probabilistically.
- Table support is limited to four built-in genetic codes; **consequence:** callers needing other NCBI tables must extend the implementation or select another tool.

**Not implemented:**

- Built-in support for NCBI genetic-code tables outside `1`, `2`, `3`, and `11`; **users should rely on:** no current alternative.
- Ambiguity-aware resolution that enumerates all amino acids implied by an IUPAC codon; **users should rely on:** no current alternative.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `AUG` | Returns `M`. | Standard start codon mapping. |
| `UAA`, `UAG`, `UGA` | Return `*` in the standard table. | Standard stop codons. |
| `ATG` | Returns `M`. | DNA is normalized to RNA. |
| Lowercase or mixed case input | Same result as uppercase. | Input is uppercased before lookup. |
| `ANN`, `NNN` | Return `X`. | Valid IUPAC ambiguity codons are treated as unknown amino acid. |
| `XYZ`, `12G` | Throw `ArgumentException`. | They contain invalid non-IUPAC symbols. |

### 6.2 Limitations

The repository exposes a fixed set of built-in genetic-code tables and returns only a single amino-acid character per codon. It does not represent ambiguity sets, probabilities, or broader codon-context effects.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GeneticCodeTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GeneticCodeTests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [TRANS-CODON-001.md](../../../tests/TestSpecs/TRANS-CODON-001.md)
- Related algorithms: [Protein_Translation.md](Protein_Translation.md)

## 8. References

1. Wikipedia contributors. 2026. Genetic code. Wikipedia. https://en.wikipedia.org/wiki/Genetic_code
2. Wikipedia contributors. 2026. Start codon. Wikipedia. https://en.wikipedia.org/wiki/Start_codon
3. Wikipedia contributors. 2026. Stop codon. Wikipedia. https://en.wikipedia.org/wiki/Stop_codon
4. NCBI. 2026. The Genetic Codes. https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
5. Test specification: [TRANS-CODON-001.md](../../../tests/TestSpecs/TRANS-CODON-001.md)
6. Crick FH. 1968. The origin of the genetic code. Journal of Molecular Biology. N/A
