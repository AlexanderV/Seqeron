# Disordered Region Detection

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinPred |
| Test Unit ID | DISORDER-REGION-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Disordered region detection in this repository converts residue-level disorder calls into contiguous intrinsically disordered regions and annotates each emitted region with a mean score, a confidence value, and a coarse composition-bias label. The region boundaries are derived from the `PredictDisorder(...)` residue predictions returned by the TOP-IDP-based disorder scorer in the same source file ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). Region classification is deterministic and rule-based. The code comments and tests show that the composition thresholds and confidence mapping are repository heuristics informed by literature, not verbatim formulas copied from a single external specification ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs); [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Intrinsically disordered regions are protein segments that do not adopt one persistent folded structure under physiological conditions. They are typically enriched in polar or charged residues and are often discussed in the context of signaling, regulation, flexible linkers, and molecular recognition (Dunker et al. 2001; van der Lee et al. 2014). The repository models IDRs as contiguous stretches of residues whose upstream TOP-IDP scores exceed the active disorder threshold ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)).

### 2.2 Core Model

Let `predictions` be the ordered residue-level output of `PredictDisorder(...)`. The region detector performs a single forward scan and emits each maximal contiguous run of residues for which `IsDisordered` is `true`, provided that the run length is at least `minRegionLength` ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). For an emitted region of length $L$ with constituent residue scores $s_1, \dots, s_L$:

$$
\mathrm{MeanScore} = \frac{1}{L} \sum_{k=1}^{L} s_k
$$

The repository then defines confidence as a clamped, normalized distance from the TOP-IDP cutoff:

$$
\mathrm{Confidence} = \min\left(1, \max\left(0, \frac{\text{MeanScore} - 0.542}{1.0 - 0.542}\right)\right)
$$

The source comment explicitly marks this confidence formula as an internal design decision, while the `0.542` cutoff itself comes from Campen et al. (2008) ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). Region labels are assigned by the ordered rule chain `Proline-rich -> Acidic -> Basic -> Ser/Thr-rich -> Long IDR -> Standard IDR`, using `> 0.25` composition thresholds for the first four classes and `length > 30` for `Long IDR` ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs); [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). The same source comment states that the `0.25` enrichment threshold is an internal heuristic, while the subtype names are inspired by van der Lee et al. (2014) and the `>30` long-segment cutoff is aligned with Ward et al. (2004).

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Emitted regions are non-overlapping and sorted by increasing start position. | The detector scans left to right and closes one run before opening the next; the tests assert sorted, non-overlapping output ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs); [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). |
| INV-02 | Every emitted region has length greater than or equal to `minRegionLength`. | Runs are yielded only after the `length >= minLength` check, and the tests verify both exact-minimum inclusion and below-minimum exclusion ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs); [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). |
| INV-03 | Region confidence is always in `[0, 1]`. | `CalculateConfidence(...)` clamps the normalized value to that interval, and the tests assert the bound across multiple region types ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs); [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). |
| INV-04 | A homopolymeric emitted region has `MeanScore` equal to the corresponding normalized TOP-IDP residue score. | Every residue in the run contributes the same score, so the region mean equals that constant value; the tests verify this for `P`, `E`, `K`, and `S` runs ([DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | Protein sequence that is first converted into residue-level disorder calls. | `null` or empty input returns no regions because `PredictDisorder(...)` returns an empty result ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |
| `windowSize` | `int` | `21` | Upstream TOP-IDP window width used before region extraction. | Region boundaries depend on the residue-level scores produced with this window size. |
| `disorderThreshold` | `double` | `0.542` | Upstream residue-level cutoff that determines which positions enter a region run. | The detector consumes `ResiduePrediction.IsDisordered`; default matches the TOP-IDP cutoff from Campen et al. (2008). |
| `minRegionLength` | `int` | `5` | Minimum contiguous run length required for emission. | Runs shorter than this value are filtered out. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `DisorderedRegions` | `IReadOnlyList<DisorderedRegion>` | Region list returned inside `DisorderPredictionResult`. |
| `DisorderedRegions[*].Start` | `int` | 0-based inclusive start index of the emitted run. |
| `DisorderedRegions[*].End` | `int` | 0-based inclusive end index of the emitted run. |
| `DisorderedRegions[*].MeanScore` | `double` | Arithmetic mean of `DisorderScore` over the region's residues. |
| `DisorderedRegions[*].Confidence` | `double` | Repository-specific confidence value derived from the region mean score and clamped to `[0, 1]`. |
| `DisorderedRegions[*].RegionType` | `string` | One of `Proline-rich`, `Acidic`, `Basic`, `Ser/Thr-rich`, `Long IDR`, or `Standard IDR`. |

### 3.3 Preconditions and Validation

Disordered region detection is not exposed as a standalone public method; it is exercised through `DisorderPredictor.PredictDisorder(...)`, which uppercases the sequence and short-circuits on `null` or empty input ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). Region coordinates are 0-based and inclusive, as shown by the tests that assert exact start and end values for full-span, leading, central, and trailing regions ([DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). Unknown residues affect region detection only through the upstream residue scorer documented in [Disorder_Prediction.md](Disorder_Prediction.md).

## 4. Algorithm

### 4.1 High-Level Steps

1. Obtain residue-level predictions from `PredictDisorder(...)`.
2. Scan the prediction list from left to right and open a region when the first `IsDisordered` residue in a run is encountered.
3. Extend the current run while `IsDisordered` remains `true`.
4. When the run ends, emit the region only if its length is at least `minRegionLength`, then compute `MeanScore`, `Confidence`, and `RegionType`.
5. After the scan, apply the same emission rule once more if the sequence ended inside a disordered run.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The detector uses `MeanScore = scoreSum / length` for each emitted run and `Confidence = (MeanScore - 0.542) / (1.0 - 0.542)` clamped to `[0, 1]` ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). Region typing follows the ordered checks `P > 0.25`, `(E + D) > 0.25`, `(K + R) > 0.25`, `(S + T) > 0.25`, `length > 30`, else `Standard IDR`. The source comment explicitly labels the `0.25` enrichment threshold as an internal heuristic and the priority chain as an internal design decision; the test suite verifies the priority cases `Proline-rich` over `Acidic` and `Acidic` over `Basic` ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs); [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Region extraction from an existing prediction list of length `n` | `O(n)` | `O(r)` | One forward scan; `r` is the number of emitted regions. |
| End-to-end region detection via `PredictDisorder(...)` | `O(nw) + O(n)` | `O(n)` | The upstream residue scorer dominates with window size `w`; region extraction is linear afterward ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)

- `DisorderPredictor.PredictDisorder(string, int, double, int)`: public API that computes residue scores and returns the final `DisorderedRegions` list.
- `DisorderPredictor.IdentifyDisorderedRegions(List<ResiduePrediction>, double, int)`: private helper that scans contiguous disordered runs and emits region records.
- `DisorderPredictor.ClassifyDisorderedRegion(List<ResiduePrediction>)`: private helper that assigns the region label from amino-acid composition.
- `DisorderPredictor.CalculateConfidence(double)`: private helper that converts `MeanScore` into a clamped confidence value.

**Supporting tests:** [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)

### 5.2 Current Behavior

The repository computes region boundaries from the `IsDisordered` flag already stored in each `ResiduePrediction`; the `threshold` argument passed into `IdentifyDisorderedRegions(...)` is forwarded from the caller but is not used inside the scan itself ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). `Start` and `End` are emitted as inclusive indices. A trailing run is handled after the loop so regions that reach the last residue are not dropped. `Confidence` depends only on `MeanScore`, not on region length, and the private helpers are tested indirectly through the public `PredictDisorder(...)` entry point ([DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Contiguous-run detection over residue-level disorder calls, with emission only for runs whose length meets `minRegionLength`.
- The `>30` cutoff used for the `Long IDR` label, which the source comments align with Ward et al. (2004) and van der Lee et al. (2014).
- Deterministic output ordering and inclusive run boundaries verified by the tests.

**Intentionally simplified:**

- The `0.25` enrichment cutoff and the confidence formula are repository heuristics rather than verbatim published formulas; **consequence:** `RegionType` and `Confidence` are stable within this codebase but should not be treated as standardized cross-tool annotations ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)).

**Not implemented:**

- Learned or probabilistic IDR subtype classifiers, or alternate region ontologies beyond the six hard-coded labels; **users should rely on:** no current alternative in this repository.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Composition enrichment cutoff `> 0.25` | Assumption | Region labels depend on a repository-defined threshold rather than a universal literature standard. | accepted | The source comment explicitly marks this value as an internal heuristic in `ClassifyDisorderedRegion(...)` ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |
| 2 | Confidence normalization from cutoff to `1.0` | Deviation | Confidence values are specific to this implementation and are not a published disorder-calibration scale. | accepted | The source comment states that the formula itself is an internal design decision in `CalculateConfidence(...)` ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| All ordered sequence such as `new string('W', 30)` | Emits no disordered regions. | Verified by the tests for the lowest TOP-IDP anchor sequence ([DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). |
| Fully disordered sequence such as `new string('P', 30)` | Emits one region spanning the full sequence. | The tests assert the exact `[0, 29]` bounds for this case ([DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). |
| Trailing disordered run | The final run is emitted after the scan ends. | Explicit end-of-loop handling in `IdentifyDisorderedRegions(...)`; tested with `W10 + P20` producing `[11, 29]` ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs); [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). |
| Run shorter than `minRegionLength` | The run is excluded from `DisorderedRegions`. | Verified by the exact-minimum and below-minimum tests ([DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). |

### 6.2 Limitations

Region detection depends entirely on the upstream TOP-IDP residue calls produced by `PredictDisorder(...)`, so any residue-level false positive or false negative changes the emitted boundaries. The region API is private and not available as a standalone entry point separate from the predictor. The six `RegionType` labels are coarse heuristic categories, not a full standardized ontology for intrinsically disordered segments.

## 7. Examples and Related Material

### 7.1 Worked Example

The region tests include two compact examples that illustrate both boundary detection and label assignment ([DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)):

- `new string('W', 10) + new string('P', 20)` yields one trailing region `[11, 29]` with the default `windowSize = 21` and `minRegionLength = 5`.
- `string.Concat(Enumerable.Repeat("EKQSP", 8))` yields one `Long IDR` region because no tracked composition exceeds `0.25` while the run length is `40 > 30`.

## 8. References

1. Campen A., et al. (2008). "TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder." Protein and Peptide Letters 15(9):956-963.
2. Das R. K., Pappu R. V. (2013). "Conformations of intrinsically disordered proteins are influenced by linear sequence distributions of oppositely charged residues." Proceedings of the National Academy of Sciences of the United States of America 110(33):13392-13397.
3. Dunker A. K., et al. (2001). "Intrinsically disordered protein." Journal of Molecular Graphics and Modelling 19(1):26-59.
4. van der Lee R., et al. (2014). "Classification of intrinsically disordered regions and proteins." Chemical Reviews 114(13):6589-6631.
5. Ward J. J., et al. (2004). "Prediction and functional analysis of native disorder in proteins from the three kingdoms of life." Journal of Molecular Biology 337(3):635-645.
6. [DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)
7. [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)
