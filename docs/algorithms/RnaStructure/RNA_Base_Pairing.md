# RNA Base Pairing

| Field | Value |
|-------|-------|
| Algorithm Group | RnaStructure |
| Test Unit ID | RNA-PAIR-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

RNA base pairing determines whether two ribonucleotides can form a hydrogen-bonded pair, classifies that pair as Watson-Crick or wobble, and computes the RNA complement of a single base. These are O(1) primitives underpinning all higher-level RNA secondary-structure algorithms (stem-loop finding, free-energy evaluation, MFE folding). The behavior is specification-driven and exact: a pair forms only for the canonical Watson-Crick pairs A•U and G•C [2] and the G•U wobble pair [1][3]; all other combinations over the standard alphabet do not pair.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

In RNA, nucleotides pair through hydrogen bonds on their Watson-Crick edges. The two canonical (Watson-Crick) pairs are A•U (2 hydrogen bonds) and G•C (3 hydrogen bonds); in RNA uracil (U) replaces thymine [2]. Beyond the canonical pairs, Crick's 1966 wobble hypothesis established that guanine can pair with uracil (G•U) in addition to cytosine, and uracil with guanine in addition to adenine [1][3]. The G•U wobble pair is thermodynamically comparable to a Watson-Crick pair but geometrically distinct and is classified separately [3].

### 2.2 Core Model

For two bases `b1`, `b2` over the RNA alphabet {A, C, G, U} (T normalized to U):

- **Watson-Crick:** {A,U}, {G,C} — pair, type WatsonCrick [2].
- **Wobble:** {G,U} — pair, type Wobble [1][3].
- **All others** (e.g., A•A, A•G, A•C, C•U, G•G, C•C) — do not pair [1] (A pairs only with U; C pairs only with G).

`CanPair(b1,b2) = true` iff the unordered pair is one of {A,U}, {G,C}, {G,U}. `GetBasePairType` returns WatsonCrick, Wobble, or null accordingly.

The RNA complement of a single base maps A→U, U→A, G→C, C→G, with T treated as U (complement A), and IUPAC degenerate codes complemented per the IUPAC-IUB table (N→N, R→Y, Y→R, S→S, W→W, K→M, M→K, B→V, D→H, H→D, V→B) [4][5].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `CanPair(x,y) == CanPair(y,x)` | Base pairing is reciprocal; the lookup table is populated symmetrically [2] |
| INV-02 | `GetBasePairType(x,y) == GetBasePairType(y,x)` | Same symmetric table [2] |
| INV-03 | `CanPair(x,y) == (GetBasePairType(x,y) != null)` | A pair exists iff it has a defined type (shared lookup) |
| INV-04 | `GetBasePairType('G','U')` and `('U','G')` are Wobble, never WatsonCrick | Wobble pairs do not follow Watson-Crick rules [1][3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| base1, base2 | char | required | RNA bases to test | case-insensitive; A/C/G/U, T treated as U |
| base_ (GetComplement) | char | required | RNA base to complement | case-insensitive; IUPAC codes supported |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| CanPair | bool | true iff the two bases form a Watson-Crick or wobble pair |
| GetBasePairType | BasePairType? | WatsonCrick, Wobble, or null (no pair) |
| GetComplement | char | RNA complement of the input base |

### 3.3 Preconditions and Validation

Inputs are upper-cased before lookup (case-insensitive). `CanPair`/`GetBasePairType` operate over the RNA alphabet {A, C, G, U}; the DNA base T is not an RNA base and does not pair (returns `false`/`null`). For `GetComplement`, T is treated as U and complements to A [5]. Characters outside the 0–127 ASCII range return `false`/`null` from `CanPair`/`GetBasePairType` (bounds-checked lookup, no exception). `GetComplement` passes unrecognized non-IUPAC characters through unchanged.

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case both inputs.
2. For `CanPair`/`GetBasePairType`: index a precomputed 128×128 lookup table by `(b1, b2)`; value 0 = no pair, 1 = Watson-Crick, 2 = Wobble.
3. For `GetComplement`: switch over the IUPAC code to its RNA complement.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The pair lookup table is seeded with the six ordered pairs A→U, U→A, G→C, C→G (Watson-Crick) and G→U, U→G (wobble) [1][2][3]. The complement table follows the IUPAC-IUB nucleotide notation, RNA variant [4][5].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CanPair / GetBasePairType | O(1) | O(1) | single table lookup; table is a fixed 16 KB byte array |
| GetComplement | O(1) | O(1) | switch expression |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [RnaSecondaryStructure.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs)

- `RnaSecondaryStructure.CanPair(char, char)`: boolean pairing test via the `PairLookup` table.
- `RnaSecondaryStructure.GetBasePairType(char, char)`: returns `BasePairType?` (WatsonCrick / Wobble / null).
- `RnaSecondaryStructure.GetComplement(char)`: delegates to [`SequenceExtensions.GetRnaComplementBase`](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs).

### 5.2 Current Behavior

`CanPair` and `GetBasePairType` share a single precomputed `byte[128*128]` lookup (`PairLookup`), giving branch-free O(1) classification and guaranteeing the symmetry invariants by construction. `GetComplement` is a thin delegate to the Core IUPAC RNA complement helper, which is cross-verified against Biopython `complement_rna` [5]. No substring search is involved, so the repository suffix tree is not applicable to this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Watson-Crick pairs A•U and G•C classified as WatsonCrick [2].
- G•U / U•G classified as Wobble, distinct from Watson-Crick [1][3].
- All non-canonical combinations over {A,C,G,U} return no pair [1].
- RNA complement table A→U, U→A, G→C, C→G, T→A, plus IUPAC degenerate codes [4][5].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Non-canonical pairs beyond G•U wobble (e.g., Hoogsteen pairs, sheared G•A, the inosine wobble pairs I•U/I•A/I•C); **users should rely on:** no current alternative in this library — only the standard four-letter RNA alphabet pairs are modeled, consistent with the nearest-neighbor folding model used by sibling methods.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| G•U / U•G | Wobble (not WatsonCrick) | Wobble pairs do not follow WC rules [1][3] |
| Lowercase input | same as uppercase | case-insensitive normalization |
| DNA T in CanPair | does not pair (T is not an RNA base) | RNA pairing defined over {A,C,G,U} [1][2] |
| DNA T in GetComplement | treated as U; complement is A | Biopython complement_rna [5] |
| Out-of-ASCII char in CanPair | false / null, no exception | bounds-checked table lookup |
| Non-IUPAC char in GetComplement | passed through unchanged | Core helper contract |

### 6.2 Limitations

Models only the standard RNA alphabet and the three pair types (Watson-Crick A•U/G•C, wobble G•U). Inosine and other modified bases, and non-canonical pair geometries (Hoogsteen, sheared, etc.), are out of scope.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
RnaSecondaryStructure.CanPair('G', 'U');          // true  (wobble)
RnaSecondaryStructure.GetBasePairType('G', 'U');  // BasePairType.Wobble
RnaSecondaryStructure.GetBasePairType('A', 'U');  // BasePairType.WatsonCrick
RnaSecondaryStructure.GetBasePairType('A', 'G');  // null  (no pair)
RnaSecondaryStructure.GetComplement('A');         // 'U'
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [RnaSecondaryStructure_CanPair_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/RnaSecondaryStructure_CanPair_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [RNA-PAIR-001-Evidence.md](../../../docs/Evidence/RNA-PAIR-001-Evidence.md)

## 8. References

1. Crick, F.H.C. (1966). Codon–anticodon pairing: the wobble hypothesis. *Journal of Molecular Biology* 19(2):548–555. https://doi.org/10.1016/S0022-2836(66)80022-0
2. Wikipedia. Base pair. https://en.wikipedia.org/wiki/Base_pair
3. Wikipedia. Wobble base pair. https://en.wikipedia.org/wiki/Wobble_base_pair
4. IUPAC-IUB Commission on Biochemical Nomenclature (1970). Abbreviations and symbols for nucleic acids, polynucleotides, and their constituents. *Biochemistry* 9(20):4022–4027. https://en.wikipedia.org/wiki/Nucleic_acid_notation
5. Biopython. Bio.Seq.complement_rna. https://biopython.org/docs/latest/api/Bio.Seq.html
