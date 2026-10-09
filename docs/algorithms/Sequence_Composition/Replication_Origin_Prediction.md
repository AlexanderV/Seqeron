# Replication Origin Prediction (Cumulative GC-Skew Minimum)

| Field | Value |
|-------|-------|
| Algorithm Group | Sequence Composition |
| Test Unit ID | SEQ-REPLICATION-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-10-09 |

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

This is Grigoriev's cumulative skew — the running sum of (G−C)/(G+C) over adjacent windows [1] — at a one-base window (each window's skew is +1, −1 or 0); the implementation computes it with the canonical `CalculateCumulativeGcSkew` kernel at window 1. There are *n*+1 prefix values Skew_0 … Skew_n. The **Minimum Skew Problem** (Rosalind BA1F) asks for all positions *i* ∈ [0, n] minimizing Skew_i [2]; the minimizing position(s) predict the replication origin, and (symmetrically) the maximizing position(s) predict the terminus [1][4].

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
| INV-07 | `FindMinimumSkewPositions` = all minimizers ascending (BA1F answer); its first element = PredictedOrigin; likewise `FindMaximumSkewPositions` / PredictedTerminus | Definition [2] |
| INV-08 | Circular mode: positions in [0, n−1] (n ≡ 0); if Skew_n = 0, rotating the input left by r maps each extremum p to (p − r) mod n | Skew'_j = Skew_{(j+r) mod n} − Skew_r when the walk closes [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `DnaSequence` or `string` | required | DNA sequence in genome coordinates (typically a complete chromosome) | DNA alphabet; case-insensitive; only G/C affect the result |
| circular | `bool` | false | (overloads `PredictReplicationOrigin(seq, bool)`, `Find{Minimum,Maximum}SkewPositions(seq, circular)`) identify prefix index n with 0 and report positions mod n | — |
| windowSize | `int` | — | (overload `PredictReplicationOrigin(seq, int)`) Grigoriev windowed cumulative skew | ≥ 1 |

### 3.2 Output / Return Value

`ReplicationOriginPrediction` (record struct):

| Field | Type | Description |
|-------|------|-------------|
| PredictedOrigin | `int` | First 0-based prefix index minimizing the cumulative skew (predicted *ori*) |
| PredictedTerminus | `int` | First 0-based prefix index maximizing the cumulative skew (predicted *ter*) |
| OriginSkew | `double` | Cumulative skew value at the minimum (≤ 0) |
| TerminusSkew | `double` | Cumulative skew value at the maximum (≥ 0) |
| IsSignificant | `bool` | True when the diagram has non-zero amplitude (max > min) |

Additional entry points (finisher 2026-10-09, additive):

- `FindMinimumSkewPositions(seq, circular = false)` / `FindMaximumSkewPositions(...)` → `IReadOnlyList<int>`: **all** minimizing / maximizing prefix indices, ascending (the full BA1F answer; never empty). Null/empty string → `[0]`.
- `PredictReplicationOrigin(seq, bool circular)`: as the base method; circular = true visits Skew_0 … Skew_{n−1} only (index n ≡ 0, the junction carries Skew_0 = 0), so positions lie in [0, n−1]. Let D = Skew_n (total #G − #C). D = 0: rotation-equivariant (each extremum moves to (p − r) mod n). D ≠ 0: the rotated walk is Skew'_j = Skew_{j+r} − Skew_r for j + r < n and Skew_{j+r−n} + D − Skew_r beyond the junction, so the extrema can change with the start (e.g. `CCGGG` → minimizers {2}; rotated by 3 → {0, 4}, i.e. {3, 2} mapped back). No detrending is applied — supply the chromosome from its coordinate 0.
- `PredictReplicationOrigin(seq, int windowSize)`: Grigoriev's windowed diagram — cumulative sum of (G−C)/(G+C) over adjacent complete windows (the points of `CalculateCumulativeGcSkew(seq, windowSize)`); origin/terminus = `Position` (window centre start + w/2) of the first min/max point, skews = those cumulative values, `IsSignificant` = max > min. There is no Skew_0 baseline point, so `OriginSkew` may be > 0; shorter-than-one-window input → zero prediction.
- `CalculateSkewIndex(seq, windowSize = 20000)` → `double?`: SkewIT Skew Index [5][6] (see 5.5).
- `ParseSkewIGenusThresholds(TextReader)` → `IReadOnlyDictionary<string, double>`, `TryGetSkewIThreshold(table, genus, out threshold)`, `IsSkewIBelowGenusThreshold(seq, genus, table, windowSize = 20000)` and `IsSkewIBelowThreshold(seq, threshold, windowSize = 20000)` → `bool?`: SkewIT's per-genus significance rule (see 5.5). The table is supplied by the caller (not bundled: SkewIT is GPL-3.0).

### 3.3 Preconditions and Validation

Indexing is 0-based over prefix indices [0, n] (position *i* refers to the boundary before base *i*), matching Rosalind BA1F [2]. Input is case-insensitive (lowercase is upper-cased). Only G and C contribute; A, T, and any other symbol leave the running skew unchanged. The `DnaSequence` overload throws `ArgumentNullException` for null. The `string` overload returns a zero prediction (`PredictedOrigin = PredictedTerminus = 0`, skews 0, `IsSignificant = false`) for null or empty input.

## 4. Algorithm

### 4.1 High-Level Steps

1. Initialize `cumulative = 0` (Skew_0) and track the running min/max value with their first prefix indices, both starting at 0.
2. Scan each base: G adds +1, C subtracts 1, A/T/other add 0; the value after consuming base *i* is Skew_{i+1}.
3. Update min (and its index) on a strictly smaller value, and max (and its index) on a strictly larger value — strict comparison keeps the **first** extreme index (tie-break).
4. Origin = min index, terminus = max index; `IsSignificant = max > min`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Per-base skew increment table [2]: G → +1, C → −1, A/T (and any non-G/C symbol) → 0. These are the only constants of the per-base walk; the windowed overload takes a window size, and the SkewI thresholds come from SkewIT's per-genus table supplied by the caller (5.5).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| PredictReplicationOrigin | O(n) | O(1) | single pass; min/max tracked in scalars, no array materialized |
| Find{Minimum,Maximum}SkewPositions | O(n) | O(#ties) | same fold, tie lists kept |
| PredictReplicationOrigin(seq, windowSize) | O(n) | O(1) | fold over windowed cumulative points |
| CalculateSkewIndex | O(n + L·r) | O(L) | L = ⌈n/k⌉ windows, r = round(0.04·L); prefix sums |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GcSkewCalculator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GcSkewCalculator.cs)

- `GcSkewCalculator.PredictReplicationOrigin(DnaSequence)`: canonical method; predicts origin/terminus from the cumulative skew diagram.
- `GcSkewCalculator.PredictReplicationOrigin(string)`: thin overload; upper-cases the input and delegates to the same core; null/empty → zero prediction.
- `GcSkewCalculator.PredictReplicationOrigin(DnaSequence|string, bool circular)`, `PredictReplicationOrigin(DnaSequence|string, int windowSize)`.
- `GcSkewCalculator.FindMinimumSkewPositions` / `FindMaximumSkewPositions(DnaSequence|string, bool circular = false)`.
- `GcSkewCalculator.CalculateSkewIndex(DnaSequence|string, int windowSize = 20000)`.
- `GcSkewCalculator.ParseSkewIGenusThresholds`, `TryGetSkewIThreshold`, `IsSkewIBelowGenusThreshold`, `IsSkewIBelowThreshold`.
- MCP `predict_replication_origin` (Analysis server): optional `circular`; optional `windowSize` → the windowed overload (linear only: `circular: true` with `windowSize` is rejected with `ArgumentException`, since the windowed Grigoriev core has no circular form; `originPositions`/`terminusPositions` are then empty); optional `skewIndexWindow` → `skewIndex` = `CalculateSkewIndex(dna, k)` (null when not requested or when skewi.py prints nothing).

### 5.2 Current Behavior

A single O(1)-space pass folds over the canonical cumulative-skew iterator (`CalculateCumulativeGcSkewCore`, window 1) without materializing the diagram. The sequence is read linearly from index 0; for a circular chromosome prefix index *n* is the same junction as 0, and when the chromosome's total #G−#C is 0 (Skew_n = 0) rotating the start only shifts the extrema by the rotation offset (verified on a synthetic genome, test R2); when Skew_n ≠ 0 the wrap-around step of size Skew_n can move the reported extremum, so supply the sequence starting at its annotated position (ASM-02). The `circular` overloads apply the same walk but report positions mod n (3.2). `Find*SkewPositions` share the same fold and additionally collect the tie list (O(#ties) space). Ties for the extreme value resolve to the smallest (first) prefix index via strict `<` / `>` comparisons. This is not a substring-search/pattern-matching task (it is a running scalar fold over the sequence), so the repository suffix tree is **not** applicable and is not used.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Per-nucleotide cumulative skew with s(G)=+1, s(C)=−1, s(A)=s(T)=0 and Skew_0 = 0 [2].
- Origin = global minimum, terminus = global maximum of the cumulative diagram [1][2][4].
- 0-based prefix indexing over [0, n]; reproduces the Rosalind BA1F sample output `53 97` (first minimizer 53) [2].
- All minimizers (BA1F answer) via `FindMinimumSkewPositions`: sample → `53 97`; BA1F extra dataset (93 523 bp) → `89969 89970 89971 90345 90346` [2].
- Grigoriev's windowed cumulative diagram prediction (`PredictReplicationOrigin(seq, windowSize)`), equal to `numpy.cumsum(Bio.SeqUtils.GC_skew(seq, w))` extrema [1].
- Circular-chromosome reporting (positions mod n) [1].
- SkewIT Skew Index, computed as the authors' `skewi.py` [5][6].
- SkewIT significance rule: SkewI below the genus threshold (mean − 2 SD) flags an atypical / possibly mis-assembled genome, using the authors' `RefSeq97_Bacteria_GenusSkewIThresholds.txt` format [5][6].

**Intentionally simplified:**

- `ReplicationOriginPrediction` holds a single origin and terminus position (first extreme index); the full tie sets are returned by `FindMinimumSkewPositions` / `FindMaximumSkewPositions`.
- `IsSignificant` uses the threshold-free predicate `max > min` rather than a quantitative confidence measure (unchanged). **Consequence:** any non-flat diagram is flagged significant; for SkewIT's sourced quantitative test use `IsSkewIBelowGenusThreshold` with SkewIT's per-genus table, or `IsSkewIBelowThreshold` with an explicit threshold (5.5) — no universal cutoff is published.

**Not implemented:**

- Multi-origin detection, detrending of an unbalanced (Skew_n ≠ 0) circular walk, and strand re-orientation; **users should rely on:** dedicated tools (e.g. SkewDB / oriC predictors) for noisy or non-canonical genomes.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | `IsSignificant` threshold | Assumption | Determines the boolean flag for any input | accepted | No authoritative numeric cutoff exists; the previous invented `amplitude > count × 0.01` constant was removed and replaced with the threshold-free `max > min` (Evidence §Assumptions 1) |

### 5.5 Skew Index (SkewIT)

`CalculateSkewIndex` reproduces `src/skewi.py` of SkewIT [6] (Lu & Salzberg 2020 [5]) step by step: windows of k bases starting at 0, k, 2k, … (last one partial) each get sign(#G − #C); with L windows, h = round(L/2) and r = round(0.04·L) (Python 3 half-to-even rounding), the sign list is doubled and maxDiff = max over i ∈ [0, L), t ∈ [i+h−r, i+h+r) of |Σ skew[i:t] − Σ skew[t:i+L]|; SkewI = min(1, maxDiff / n · k). skewi.py reports no value when maxDiff ≤ 0 (no G/C-signed window, or L ≤ 12 so r = 0) → `null`. Computed with prefix sums (O(L·r)) instead of the script's O(L²·r) slice sums — identical integers. CLI-only input filters of skewi.py are not applied (default `--min-len` 500 kb, "complete" required / "plasmid" excluded in the FASTA header; `-f` is parsed but unused by its computation). skewi.py counts only upper-case `G`/`C`; this method upper-cases (identical on RefSeq FASTA). Thresholds: SkewIT publishes per-genus thresholds = genus mean − 2 SD for genera with ≥ 10 RefSeq-97 genomes (`data/RefSeq97_Bacteria_GenusSkewIThresholds.txt`, e.g. Escherichia 0.7110, Salmonella 0.8478); a value below its genus threshold flags a possibly mis-assembled genome. There is no default cutoff, so none is built in.

**Per-genus test (finisher A2-3).** `ParseSkewIGenusThresholds(TextReader)` reads SkewIT's table as published (tab-separated, header `Genus Num_Genomes Mean STDEV Threshold`, `g__` genus prefix stripped, CRLF or LF; rows with an empty threshold — genera with < 10 genomes — are skipped; malformed threshold or duplicate genus → `FormatException`). The RefSeq-97 file has 1 147 genus rows, of which 160 carry a threshold (e.g. Escherichia 0.7110, Bordetella 0.2200, Mycobacterium 0.3959, Streptomyces 0.046; the smallest is Synechococcus −0.222). `TryGetSkewIThreshold` matches the exact genus name case-insensitively (an optional `g__` prefix is ignored). `IsSkewIBelowThreshold(seq, threshold, k)` returns `SkewI < threshold` (strict, as SkewIT's "below the threshold"), or `null` when SkewI is null. `IsSkewIBelowGenusThreshold(seq, genus, table, k)` also returns `null` when the genus has no threshold. `IsSignificant` is unchanged.

The table is **not bundled**. The SkewIT repository is licensed GPL-3.0 (its `LICENSE` file), while Seqeron is MIT, so the file cannot be embedded in the library without imposing GPL terms. Download it from https://github.com/jenniferlu717/SkewIT/blob/master/data/RefSeq97_Bacteria_GenusSkewIThresholds.txt and pass it to `ParseSkewIGenusThresholds`, or pass an explicit threshold.

**Window size.** The README gives skewi.py's default as non-overlapping 20 kb windows (`-k 20000`, minimum sequence length 500 kb). It names no other k for the RefSeq-97 data, so the published thresholds correspond to k = 20 000, the default of these methods. A SkewI computed with another k is not comparable with them.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| No G/C (e.g. `AAAATTTT`) | origin = terminus = 0, skews 0, `IsSignificant = false` | Flat diagram; Skew stays at Skew_0 = 0 [2] |
| Tied extremum (e.g. `CCGGCC`) | first minimizing/maximizing index reported | INV-01/02 tie-break |
| Single base `G` | origin 0 (skew 0), terminus 1 (skew +1) | Diagram `0, +1` [2] |
| Null `DnaSequence` | `ArgumentNullException` | Input validation |
| Null/empty `string` | zero prediction, not significant | Documented overload behavior |
| Circular `GGGCCC` | minimizers {0} (linear {0, 6}) | n ≡ 0 |
| Circular, Skew_n ≠ 0 | result depends on the start | 5.2 / 3.2 |
| SkewI with ≤ 12 windows or no signed window | `null` | skewi.py prints nothing |
| Genus without a threshold (< 10 genomes, or absent) | `IsSkewIBelowGenusThreshold` → `null` | SkewIT publishes no threshold |

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

- Tests: [GcSkewCalculator_PredictReplicationOrigin_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GcSkewCalculator_PredictReplicationOrigin_Tests.cs) — covers `INV-01`…`INV-06`; [GcSkewCalculator_ReplicationOriginExtensions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GcSkewCalculator_ReplicationOriginExtensions_Tests.cs) — `INV-07`, `INV-08`, windowed; [GcSkewCalculator_SkewIndex_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GcSkewCalculator_SkewIndex_Tests.cs) — SkewI; [GcSkewCalculator_SkewIThreshold_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/GcSkewCalculator_SkewIThreshold_Tests.cs) — per-genus thresholds
- Evidence: [SEQ-REPLICATION-001-Evidence.md](../../../docs/Evidence/SEQ-REPLICATION-001-Evidence.md)
- Related algorithms: [GC_Skew](./GC_Skew.md), [AT_Skew](../Extended_GC_Skew_Analysis/AT_Skew.md)

## 8. References

1. Grigoriev, A. 1998. Analyzing genomes with cumulative skew diagrams. Nucleic Acids Research 26(10):2286–2290. https://doi.org/10.1093/nar/26.10.2286
2. Rosalind. Minimum Skew Problem (BA1F). https://rosalind.info/problems/ba1f/
3. Lobry, J. R. 1996. Asymmetric substitution patterns in the two DNA strands of bacteria. Molecular Biology and Evolution 13(5):660–665. https://pubmed.ncbi.nlm.nih.gov/8676740/
4. Wikipedia. GC skew. https://en.wikipedia.org/wiki/GC_skew
5. Lu, J.; Salzberg, S. L. 2020. SkewIT: The Skew Index Test for large-scale GC Skew analysis of bacterial genomes. PLoS Computational Biology 16(12):e1008439. https://doi.org/10.1371/journal.pcbi.1008439
6. SkewIT source code, `src/skewi.py`, README, `data/RefSeq97_Bacteria_GenusSkewIThresholds.txt`. https://github.com/jenniferlu717/SkewIT
