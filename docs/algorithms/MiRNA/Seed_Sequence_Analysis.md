# Seed Sequence Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | MiRNA |
| Test Unit ID | MIRNA-SEED-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Seed sequence analysis extracts and compares the canonical 5' seed region of a mature miRNA, the short segment that primarily determines animal target recognition specificity [1][2][3]. The repository implements this as exact substring extraction of positions 2-8, record construction in `CreateMiRna`, and exact seed comparison in `CompareSeedRegions`. The implementation is deterministic and lightweight: it does not predict binding sites itself, but it provides the normalized seed representation consumed by the repository's target-site workflow. Functional miRNA family membership is approximated by exact equality of the stored 7-nt seed string [3][4].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

For animal miRNAs, pairing between the miRNA 5' end and the target transcript is dominated by the seed region, usually defined as nucleotides 2-8 of the mature miRNA [1][2][3]. This short sequence is the basis for canonical site classes such as 8mer, 7mer-m8, 7mer-A1, and 6mer sites in TargetScan-style prediction frameworks [1][3]. miRNAs with the same canonical seed are commonly grouped into the same family because they tend to regulate overlapping target sets [3][4].

### 2.2 Core Model

For a mature miRNA sequence $m = m_1m_2\ldots m_n$ with $n \ge 8$, the canonical seed used here is:

$$
\operatorname{seed}(m) = m_2m_3m_4m_5m_6m_7m_8
$$

using biological 1-based indexing [1][2][3]. The repository compares two seeds by exact position-wise character matching, equivalent to Hamming comparison for equal-length seeds:

$$
\operatorname{matches}(a,b) = \sum_{i=1}^{k} [a_i = b_i], \qquad
\operatorname{mismatches}(a,b) = \sum_{i=1}^{k} [a_i \ne b_i]
$$

with `IsSameFamily = true` only when the stored seed strings are identical.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `GetSeedSequence` returns either the empty string or a 7-character uppercase seed. | The method returns `""` for inputs shorter than 8 nt and otherwise calls `Substring(1, 7).ToUpperInvariant()`. |
| INV-02 | `CreateMiRna(...).SeedSequence = GetSeedSequence(CreateMiRna(...).Sequence)`. | `CreateMiRna` normalizes the stored sequence first, then computes the seed from that normalized sequence. |
| INV-03 | `CreateMiRna` stores `SeedStart = 1` and `SeedEnd = 7`. | The record is constructed with fixed indices corresponding to zero-based positions 1 through 7. |
| INV-04 | `CompareSeedRegions(...).IsSameFamily` is true if and only if the two stored seed strings are exactly equal. | The implementation sets `isSameFamily = seed1 == seed2`. |
| INV-05 | When both seeds are present and canonical, `Matches + Mismatches = 7`. | Comparison is character-by-character over the two 7-nt seed strings. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `miRnaSequence` | `string` | required | Mature miRNA sequence passed to `GetSeedSequence`. | Inputs shorter than 8 nt yield `""`; casing is normalized to uppercase only. |
| `name` | `string` | required | miRNA identifier passed to `CreateMiRna`. | Stored verbatim; no validation is performed. |
| `sequence` | `string` | required | miRNA sequence passed to `CreateMiRna`. | Expected to be non-null; `CreateMiRna` uppercases and converts `T` to `U`. |
| `mirna1`, `mirna2` | `MiRna` | required | Seed-bearing records compared by `CompareSeedRegions`. | Empty stored seeds produce a zeroed comparison result. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `seed` | `string` | Return value of `GetSeedSequence`: uppercase 7-nt canonical seed or `""` when the input is too short. |
| `Name` | `string` | miRNA identifier stored in `MiRna`. |
| `Sequence` | `string` | Normalized uppercase RNA sequence stored in `MiRna`; `T` is converted to `U` by `CreateMiRna`. |
| `SeedSequence` | `string` | Stored canonical seed used by downstream target-site prediction. |
| `SeedStart` | `int` | Zero-based inclusive start index of the stored seed in `Sequence`. |
| `SeedEnd` | `int` | Zero-based inclusive end index of the stored seed in `Sequence`. |
| `Matches` | `int` | Number of equal seed positions reported by `SeedComparison`. |
| `Mismatches` | `int` | Number of unequal positions, plus any seed-length difference. |
| `IsSameFamily` | `bool` | Exact-seed family flag reported by `SeedComparison`. |

### 3.3 Preconditions and Validation

`GetSeedSequence` is defensive: `null`, empty, and shorter-than-8-nt inputs return `""` instead of throwing. `CreateMiRna` does not guard `sequence` against `null`; it expects a valid string and immediately normalizes it with `ToUpperInvariant()` and `T→U` replacement. `CompareSeedRegions` does not throw when either stored seed is empty; instead it returns `Matches = 0`, `Mismatches = 0`, and `IsSameFamily = false`. Biological seed positions 2-8 correspond to zero-based indices 1-7 in the stored record.

## 4. Algorithm

### 4.1 High-Level Steps

1. For direct seed extraction, uppercase the input miRNA string and return characters at zero-based indices 1 through 7 when available.
2. For `CreateMiRna`, normalize the sequence to uppercase RNA by replacing `T` with `U`.
3. Extract the normalized canonical seed and store it together with fixed seed bounds `(1, 7)` in a `MiRna` record.
4. For seed comparison, read the stored seed strings from both records.
5. Count exact positional matches and mismatches, then mark the pair as the same family only if the two seed strings are identical.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `GetSeedSequence` | `O(1)` | `O(1)` | Operates on a fixed-width 7-nt window once the string is present. |
| `CreateMiRna` | `O(n)` | `O(n)` | Dominated by full-sequence normalization before seed extraction. |
| `CompareSeedRegions` | `O(k)` | `O(1)` | `k` is the shorter seed length; canonical seeds make this effectively constant-time. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MiRnaAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs)

- `MiRnaAnalyzer.GetSeedSequence(string)`: Extracts the uppercase 7-nt canonical seed from positions 2-8.
- `MiRnaAnalyzer.CreateMiRna(string, string)`: Normalizes an input sequence and materializes the `MiRna` record used elsewhere in the analyzer.
- `MiRnaAnalyzer.CompareSeedRegions(MiRna, MiRna)`: Computes exact positional matches, mismatches, and exact-seed family membership.

### 5.2 Current Behavior

`GetSeedSequence` uppercases the supplied string but does not convert DNA `T` to RNA `U`; that normalization is performed only by `CreateMiRna`. `CreateMiRna` always stores `SeedStart = 1` and `SeedEnd = 7`, so the record uses zero-based coordinates even though the biological description of the seed uses 1-based positions 2-8. `CompareSeedRegions` counts length differences as additional mismatches, even though canonical seeds generated through `CreateMiRna` are either 7 nt or empty. Empty stored seeds produce a fully zeroed comparison result rather than a partial comparison.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Canonical seed extraction from miRNA positions 2-8 [1][2][3].
- Exact use of a 7-nt seed representation for downstream canonical site classification [1][3].
- Exact character-wise seed comparison and exact-seed family equality.

**Intentionally simplified:**

- Family membership is reduced to exact equality of the stored canonical 7-mer; **consequence:** noncanonical seed classes, shifted isomiR seeds, and broader family definitions are not represented.
- Direct seed extraction does not interpret pairing or conservation context; **consequence:** the result is only the canonical seed string, not a prediction of functional targeting strength.

**Not implemented:**

- Noncanonical seed models such as offset seeds, centered pairing, or isomiR-aware family clustering; **users should rely on:** [Target_Site_Prediction.md](../../../docs/algorithms/MiRNA/Target_Site_Prediction.md) for canonical site classes and external miRNA annotation resources for broader family curation.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `GetSeedSequence(null)` | Returns `""`. | The method treats null and empty input as non-computable seed cases. |
| Input length `< 8` | Returns `""`. | Positions 2-8 do not exist unless the sequence is at least 8 nt long. |
| DNA input passed directly to `GetSeedSequence` | Returns an uppercase seed that may still contain `T`. | Direct seed extraction uppercases only. |
| DNA input passed through `CreateMiRna` | Stored sequence and seed use `U` instead of `T`. | `CreateMiRna` explicitly normalizes DNA to RNA. |
| Empty seeds in `CompareSeedRegions` | Returns `0` matches, `0` mismatches, `false` family membership. | The comparison exits early when either seed is empty. |

### 6.2 Limitations

The implementation models only the canonical 7-nt seed used by the repository's target predictor. It does not represent noncanonical seed definitions, bulged pairing, conservation, expression context, or isomiR variation. Exact-seed family membership is useful operationally inside this repository, but it is not a complete biological taxonomy of miRNA families.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
using Seqeron.Genomics.Annotation;

var mirna = MiRnaAnalyzer.CreateMiRna("hsa-let-7a-5p", "UGAGGUAGUAGGUUGUAUAGUU");
string seed = MiRnaAnalyzer.GetSeedSequence(mirna.Sequence);
var comparison = MiRnaAnalyzer.CompareSeedRegions(
    mirna,
    MiRnaAnalyzer.CreateMiRna("hsa-let-7b-5p", "UGAGGUAGUAGGUUGUGUGGUU"));

// seed == "GAGGUAG"
// comparison.IsSameFamily == true
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MiRnaAnalyzer_SeedAnalysis_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/MiRnaAnalyzer_SeedAnalysis_Tests.cs)
- Test spec: [MIRNA-SEED-001.md](../../../tests/TestSpecs/MIRNA-SEED-001.md)
- Evidence: [MIRNA-SEED-001-Evidence.md](../../../docs/Evidence/MIRNA-SEED-001-Evidence.md)
- Related property tests: [MiRnaProperties.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Properties/MiRnaProperties.cs)
- Related snapshots: [MiRnaSnapshotTests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Snapshots/MiRnaSnapshotTests.cs)

## 8. References

1. Bartel DP. 2009. MicroRNAs: target recognition and regulatory functions. Cell. 136(2):215-233.
2. Lewis BP, Burge CB, Bartel DP. 2005. Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets. Cell. 120(1):15-20.
3. TargetScan Human 8.0. 2021. Conserved and nonconserved site definitions for animal miRNA targeting. https://www.targetscan.org/
4. Kozomara A, Birgaoanu M, Griffiths-Jones S. 2019. miRBase: from microRNA sequences to function. Nucleic Acids Research. 47(D1):D155-D162.
