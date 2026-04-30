# IUPAC Degenerate Motif Matching

| Field | Value |
|-------|-------|
| Algorithm Group | Pattern Matching |
| Test Unit ID | PAT-IUPAC-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

IUPAC degenerate motif matching extends exact pattern search by allowing ambiguity codes at each motif position. In this repository, the main implementation scans DNA sequences against motifs containing the 15 standard IUPAC nucleotide symbols and returns exact window matches scored as `1.0`.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The IUPAC-IUB nucleotide notation standard encodes sets of accepted bases at each motif position, enabling searches for consensus sites, degenerate primers, and motifs that tolerate known allelic variation. The implementation targets DNA alphabets and the standard ambiguity-code set documented by IUPAC-IUB (1970) and NC-IUB (1984).

### 2.2 Core Model

A sequence window `S[i..i+m-1]` matches motif `P` if, for every position `j`, the sequence character belongs to the base set allowed by `P[j]`.

The repository recognizes the standard 15 DNA codes:

| Code | Represents | Mnemonic | Count |
|------|------------|----------|-------|
| A | A | Adenine | 1 |
| C | C | Cytosine | 1 |
| G | G | Guanine | 1 |
| T | T | Thymine | 1 |
| R | A, G | puRine | 2 |
| Y | C, T | pYrimidine | 2 |
| S | G, C | Strong | 2 |
| W | A, T | Weak | 2 |
| K | G, T | Keto | 2 |
| M | A, C | aMino | 2 |
| B | C, G, T | not A | 3 |
| D | A, G, T | not C | 3 |
| H | A, C, T | not G | 3 |
| V | A, C, G | not T | 3 |
| N | A, C, G, T | aNy | 4 |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported match satisfies the allowed-base constraint at every motif position | The scan exits the inner loop on the first disallowed character |
| INV-02 | Pattern validation accepts only the 15 standard IUPAC DNA codes | `ValidateIupacPattern(...)` checks membership in the internal code dictionary |
| INV-03 | Returned matches preserve the original window and normalized pattern | `MotifMatch` stores `MatchedSequence` and the uppercased motif |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `DnaSequence` or `string` | required | DNA sequence to scan | Null `DnaSequence` throws `ArgumentNullException`; null or empty raw string yields no matches |
| `motif` | `string` | required | Motif containing IUPAC codes | Empty motif yields no matches; invalid codes throw `ArgumentException` |
| `cancellationToken` | `CancellationToken` | none | Cancellation support for long scans | Checked every 1000 start positions in the core overload |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Position` | `int` | Zero-based start index of the matching window |
| `MatchedSequence` | `string` | Sequence substring that satisfied the motif |
| `Pattern` | `string` | Uppercased motif pattern |
| `Score` | `double` | Always `1.0` for this matching mode |

### 3.3 Preconditions and Validation

`MotifFinder.FindDegenerateMotif(...)` throws `ArgumentNullException` when the `DnaSequence` input is null. Both the standard and cancellation-aware implementations return no matches when the sequence or motif is empty. Patterns are uppercased and validated against the standard IUPAC code set before scanning. `IupacHelper.MatchesIupac(...)` throws `ArgumentOutOfRangeException` for unrecognized codes instead of silently falling back.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the motif to uppercase.
2. Validate each motif character against the IUPAC code table.
3. Slide the motif-length window across the sequence.
4. For each window, compare each base to the allowed set encoded by the motif character.
5. Emit a `MotifMatch` when all positions satisfy the IUPAC constraints.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The `MotifFinder` implementation uses an internal dictionary that maps each IUPAC code to the string of allowed bases. `IupacHelper.MatchesIupac(...)` exposes the same decision table as a switch expression. The cancellation-aware core checks for cancellation every 1000 starting positions.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Single-window check | `O(m)` | `O(1)` | `m` is motif length |
| Full scan | `O(n × m)` | `O(1)` auxiliary | Brute-force scan across `n - m + 1` windows |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs), [IupacHelper.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/IupacHelper.cs)

- `MotifFinder.FindDegenerateMotif(DnaSequence, string)`: Standard IUPAC motif scan.
- `MotifFinder.FindDegenerateMotif(..., CancellationToken)`: Cancellation-aware overloads for `DnaSequence` and raw strings.
- `IupacHelper.MatchesIupac(char, char)`: Public IUPAC-code matching helper.

### 5.2 Current Behavior

`FindDegenerateMotif(DnaSequence, string)` reads `sequence.Sequence` directly, uppercases the motif, validates it, and assigns `Score = 1.0` to every returned `MotifMatch`. The cancellation-aware string overload uppercases both the sequence and the motif before scanning. `IupacHelper.MatchesIupac(...)` accepts only the 15 standard codes and throws on unknown codes.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Standard IUPAC DNA ambiguity-code interpretation.
- Positionwise matching of sequence windows against degenerate motifs.
- Support for the canonical 15 nucleotide codes `A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V`.

**Intentionally simplified:**

- The search is a brute-force window scan; **consequence:** it is straightforward but not sublinear in the sequence length.
- Matches are reported with a fixed score of `1.0`; **consequence:** the result distinguishes match/non-match only and does not rank partial or probabilistic matches.

**Not implemented:**

- Nonstandard ambiguity symbols or fallback exact matching for invalid codes; **users should rely on:** valid IUPAC motifs only.

### 5.4 Deviations and Assumptions

No intentional deviation from the standard IUPAC DNA code table was confirmed in source. The notable implementation-specific assumption is that invalid motif characters are rejected before matching rather than being treated as literal characters.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null `DnaSequence` | Throws `ArgumentNullException` | Explicit guard |
| Empty motif | Returns no matches | Explicit guard |
| Pattern longer than sequence | Returns no matches | The outer loop does not run |
| Lowercase motif input | Normalized to uppercase | `ToUpperInvariant()` is applied before validation |
| Invalid IUPAC code | Throws `ArgumentException` in `MotifFinder` or `ArgumentOutOfRangeException` in `IupacHelper` | Invalid symbols are rejected |

### 6.2 Limitations

This workflow performs exact acceptance-set matching only. It does not score partial agreement, combine ambiguity with PWM-style weighting, or index motifs for faster repeated searches.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical / biological walk-through (optional):**

Common degenerate motifs preserved from the original document:

| Motif | Pattern | Description |
|-------|---------|-------------|
| E-box | `CANNTG` | Transcription factor binding site |
| TATA box | `TATAAA` | Promoter element |
| Kozak sequence | `GCCGCCRCCATG` | Translation initiation context |

## 8. References

1. IUPAC-IUB Commission on Biochemical Nomenclature (1970). "Abbreviations and symbols for nucleic acids, polynucleotides, and their constituents." *Biochemistry* 9(20):4022-4027. doi:10.1021/bi00822a023
2. NC-IUB (1984). "Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences." *Nucleic Acids Research* 13(9):3021-3030. doi:10.1093/nar/13.9.3021
3. Wikipedia contributors. "Nucleic acid notation." *Wikipedia, The Free Encyclopedia*. https://en.wikipedia.org/wiki/Nucleic_acid_notation
4. Bioinformatics.org. "IUPAC Codes." https://www.bioinformatics.org/sms/iupac.html
