# Replication Origin Prediction (Cumulative GC-Skew Minimum)

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition |
| Test Unit ID | SEQ-REPLICATION-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Bacterial chromosomes are replicated bidirectionally from a single origin (*ori*) to a terminus (*ter*). The leading and lagging strands accumulate different mutational/selective biases, so the leading strand becomes relatively enriched in guanine over cytosine [1][3]. Plotting the running (cumulative) difference of G and C along the sequence produces a "cumulative skew diagram" whose global **minimum** marks the replication origin and whose global **maximum** marks the terminus [1][2]. This is an exact, deterministic O(n) computation over the sequence (it is a prediction in the biological sense, but the position it returns is the precise extremum of a well-defined function).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

GC skew quantifies strand-composition asymmetry. Lobry (1996) first reported that the two strands of bacterial genomes deviate from intrastrand A=T and C=G equifrequency, and that the skew switches sign at the origin and terminus of replication [3]. Grigoriev (1998) showed that integrating the skew into a cumulative diagram makes the origin and terminus appear as the diagram's extrema, separated by about half the chromosome length [1].

### 2.2 Core Model

For a genome `Genome` of length *n*, define the cumulative skew over the prefix `Genome[0..i)` as the running difference between the number of G and C bases [2]:

```
Skew_0 = 0
Skew_{i+1} = Skew_i + s(Genome[i]),   where s(G) = +1, s(C) = -1, s(A) = s(T) = 0
```

There are *n*+1 prefix values Skew_0 … Skew_n. The **Minimum Skew Problem** (Rosalind BA1F) asks for all positions *i* ∈ [0, n] minimizing Skew_i [2]; the minimizing position(s) predict the replication origin, and (symmetrically) the maximizing position(s) predict the terminus [1][4].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The genome has a single, bidirectionally replicated origin with a measurable leading-strand G>C bias [1][3] | The diagram is flat or multi-modal; the global extremum no longer corresponds to *ori* (e.g. linear genomes, plasmids, heavily rearranged genomes) |
| ASM-02 | The sequence is supplied 5'→3' on one strand in genome coordinates | A reverse-complemented or rotated input shifts/mirrors the predicted positions |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | PredictedOrigin is the first prefix index minimizing Skew_i | Definition of the Minimum Skew Problem [2]; ties broken by smallest index |
| INV-02 | PredictedTerminus is the first prefix index maximizing Skew_i | Maximum of the cumulative diagram = terminus [1][4] |
| INV-03 | OriginSkew ≤ 0 ≤ TerminusSkew | Skew_0 = 0 is always one of the prefix values, so the min cannot exceed 0 and the max cannot fall below 0 [2] |
| INV-04 | 0 ≤ PredictedOrigin, PredictedTerminus ≤ n | Prefix indices range over [0, n] [2] |
| INV-05 | IsSignificant ⇔ max > min (non-zero amplitude) | A flat diagram (no net G/C asymmetry) carries no origin signal [1][3] (see 5.4) |
| INV-06 | A and T bases do not change the diagram | s(A) = s(T) = 0 [2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `DnaSequence` or `string` | required | DNA sequence in genome coordinates (typically a complete chromosome) | DNA alphabet; case-insensitive; only G/C affect the result |

### 3.2 Output / Return Value

`ReplicationOriginPrediction` (record struct):

| Field | Type | Description |
|-------|------|-------------|
| PredictedOrigin | `int` | First 0-based prefix index minimizing the cumulative skew (predicted *ori*) |
| PredictedTerminus | `int` | First 0-based prefix index maximizing the cumulative skew (predicted *ter*) |
| OriginSkew | `double` | Cumulative skew value at the minimum (≤ 0) |
| TerminusSkew | `double` | Cumulative skew value at the maximum (≥ 0) |
| IsSignificant | `bool` | True when the diagram has non-zero amplitude (max > min) |

### 3.3 Preconditions and Validation

Indexing is 0-based over prefix indices [0, n] (position *i* refers to the boundary before base *i*), matching Rosalind BA1F [2]. Input is case-insensitive (lowercase is upper-cased). Only G and C contribute; A, T, and any other symbol leave the running skew unchanged. The `DnaSequence` overload throws `ArgumentNullException` for null. The `string` overload returns a zero prediction (`PredictedOrigin = PredictedTerminus = 0`, skews 0, `IsSignificant = false`) for null or empty input.

## 4. Algorithm

### 4.1 High-Level Steps

1. Initialize `cumulative = 0` (Skew_0) and track the running min/max value with their first prefix indices, both starting at 0.
2. Scan each base: G adds +1, C subtracts 1, A/T/other add 0; the value after consuming base *i* is Skew_{i+1}.
3. Update min (and its index) on a strictly smaller value, and max (and its index) on a strictly larger value — strict comparison keeps the **first** extreme index (tie-break).
4. Origin = min index, terminus = max index; `IsSignificant = max > min`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Per-base skew increment table [2]: G → +1, C → −1, A/T (and any non-G/C symbol) → 0. These are the only constants; there is no threshold or window parameter.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| PredictReplicationOrigin | O(n) | O(1) | single pass; min/max tracked in scalars, no array materialized |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GcSkewCalculator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs)

- `GcSkewCalculator.PredictReplicationOrigin(DnaSequence)`: canonical method; predicts origin/terminus from the cumulative skew diagram.
- `GcSkewCalculator.PredictReplicationOrigin(string)`: thin overload; upper-cases the input and delegates to the same core; null/empty → zero prediction.

### 5.2 Current Behavior

A single O(1)-space pass computes the diagram without materializing it. Ties for the extreme value resolve to the smallest (first) prefix index via strict `<` / `>` comparisons. This is not a substring-search/pattern-matching task (it is a running scalar fold over the sequence), so the repository suffix tree is **not** applicable and is not used.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Per-nucleotide cumulative skew with s(G)=+1, s(C)=−1, s(A)=s(T)=0 and Skew_0 = 0 [2].
- Origin = global minimum, terminus = global maximum of the cumulative diagram [1][2][4].
- 0-based prefix indexing over [0, n]; reproduces the Rosalind BA1F sample output `53 97` (first minimizer 53) [2].

**Intentionally simplified:**

- The API returns a single origin and terminus position; BA1F enumerates *all* minimizers. **Consequence:** when several positions tie for the extremum, only the first (smallest index) is reported.
- `IsSignificant` uses the threshold-free predicate `max > min` rather than a quantitative confidence measure. **Consequence:** any non-flat diagram is flagged significant; callers needing a numeric confidence should inspect `TerminusSkew − OriginSkew` directly.

**Not implemented:**

- Windowed/smoothed skew, multi-origin detection, and strand re-orientation; **users should rely on:** dedicated tools (e.g. SkewDB / oriC predictors) for noisy or non-canonical genomes.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `IsSignificant` threshold | Assumption | Determines the boolean flag for any input | accepted | No authoritative numeric cutoff exists; the previous invented `amplitude > count × 0.01` constant was removed and replaced with the threshold-free `max > min` (Evidence §Assumptions 1) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| No G/C (e.g. `AAAATTTT`) | origin = terminus = 0, skews 0, `IsSignificant = false` | Flat diagram; Skew stays at Skew_0 = 0 [2] |
| Tied extremum (e.g. `CCGGCC`) | first minimizing/maximizing index reported | INV-01/02 tie-break |
| Single base `G` | origin 0 (skew 0), terminus 1 (skew +1) | Diagram `0, +1` [2] |
| Null `DnaSequence` | `ArgumentNullException` | Input validation |
| Null/empty `string` | zero prediction, not significant | Documented overload behavior |

### 6.2 Limitations

The prediction is only meaningful for genomes that satisfy ASM-01 (single bidirectional origin, leading-strand G>C bias). Linear genomes, plasmids, eukaryotic chromosomes with many origins, and heavily rearranged genomes can produce flat or multi-modal diagrams where the global extremum does not correspond to a real origin [1]. Input must be in genome coordinates on one strand (ASM-02); a rotated or reverse-complemented sequence shifts or mirrors the result.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var genome = new DnaSequence(
    "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG");
var pred = GcSkewCalculator.PredictReplicationOrigin(genome);
// pred.PredictedOrigin == 53, pred.OriginSkew == -4   (Rosalind BA1F sample: minimizers "53 97")
```

**Numerical walk-through:** for `CCGGGG` the diagram is `0, −1, −2, −1, 0, +1, +2`; the minimum −2 first occurs at prefix index 2 (origin) and the maximum +2 at prefix index 6 (terminus).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GcSkewCalculator_PredictReplicationOrigin_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GcSkewCalculator_PredictReplicationOrigin_Tests.cs) — covers `INV-01`…`INV-06`
- Evidence: [SEQ-REPLICATION-001-Evidence.md](../../../docs/Evidence/SEQ-REPLICATION-001-Evidence.md)
- Related algorithms: [GC_Skew](./GC_Skew.md), [AT_Skew](../Extended_GC_Skew_Analysis/AT_Skew.md)

## 8. References

1. Grigoriev, A. 1998. Analyzing genomes with cumulative skew diagrams. Nucleic Acids Research 26(10):2286–2290. https://doi.org/10.1093/nar/26.10.2286
2. Rosalind. Minimum Skew Problem (BA1F). https://rosalind.info/problems/ba1f/
3. Lobry, J. R. 1996. Asymmetric substitution patterns in the two DNA strands of bacteria. Molecular Biology and Evolution 13(5):660–665. https://pubmed.ncbi.nlm.nih.gov/8676740/
4. Wikipedia. GC skew. https://en.wikipedia.org/wiki/GC_skew
