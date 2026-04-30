# Disorder Prediction

| Field | Value |
|-------|-------|
| Algorithm Group | ProteinPred |
| Test Unit ID | DISORDER-PRED-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Disorder prediction in this repository estimates intrinsic disorder from amino-acid composition using the TOP-IDP propensity scale and a sliding window. For each residue, the score is the average of normalized TOP-IDP values in a local window, and residues at or above the cutoff are labeled disordered (Campen et al. 2008). The public API returns the uppercased sequence, per-residue predictions, contiguous disordered regions, overall disorder content, and mean disorder score through `DisorderPredictionResult` ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). The source remarks explicitly describe this implementation as a coarse, single-feature heuristic toolkit rather than a competitive multi-feature predictor ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Intrinsically disordered proteins and intrinsically disordered regions lack a single stable three-dimensional structure under physiological conditions and instead populate dynamic conformational ensembles. These regions are associated with signaling, regulation, and molecular recognition roles, especially where flexible binding or context-dependent folding is required (Dunker et al. 2001).

### 2.2 Core Model

The predictor uses the TOP-IDP amino-acid propensity scale introduced by Campen et al. (2008). For residue $i$, the disorder score is the average normalized propensity across a local window $W_i$:

$$
S_i = \frac{1}{|W_i|} \sum_{c \in W_i} \frac{p(c) - p_{\min}}{p_{\max} - p_{\min}}
$$

where $p_{\min} = -0.884$ for tryptophan and $p_{\max} = 0.987$ for proline (Campen et al. 2008). The repository default cutoff is $0.542$, matching the TOP-IDP decision threshold cited in the code comments and used throughout the tests (Campen et al. 2008; [DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs); [DisorderPredictor_DisorderPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs)). A residue is classified as disordered when $S_i \ge 0.542$.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Local amino-acid composition captured by the TOP-IDP scale is sufficient for first-pass disorder ranking. | Predictions can miss context that depends on evolutionary profiles, predicted structure, or trained machine-learning features, which the source remarks explicitly exclude from this implementation ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every per-residue disorder score is in `[0, 1]`. | Each recognized TOP-IDP value is normalized between the documented minimum and maximum before window averaging (Campen et al. 2008); the test suite checks this range across diverse inputs ([DisorderPredictor_DisorderPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs)). |
| INV-02 | The algorithm is deterministic for a fixed input sequence and parameter set. | The scoring path uses fixed lookup tables, arithmetic, and thresholding with no random state ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |
| INV-03 | Residues are marked disordered exactly when their averaged score is greater than or equal to the active threshold. | `CalculatePerResidueScores(...)` computes `isDisordered = score >= threshold` in the implementation ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | Protein sequence to score residue by residue. | `null` or empty input returns an empty result instead of throwing ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |
| `windowSize` | `int` | `21` | Width of the local scoring window used for TOP-IDP averaging. | No explicit argument validation is performed in the public method. |
| `disorderThreshold` | `double` | `0.542` | Score cutoff used to set `ResiduePrediction.IsDisordered`. | Default matches the TOP-IDP cutoff cited in the code comments and tests (Campen et al. 2008; [DisorderPredictor_DisorderPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs)). |
| `minRegionLength` | `int` | `5` | Minimum contiguous run length used when constructing `DisorderedRegions`. | Applied after residue-level scoring; shorter runs are omitted from the region list ([Disordered_Region_Detection.md](Disordered_Region_Detection.md)). |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Sequence` | `string` | Uppercased input sequence, or `""` when the input is `null` or empty. |
| `ResiduePredictions` | `IReadOnlyList<ResiduePrediction>` | One entry per input residue, each containing `Position`, `Residue`, `DisorderScore`, and `IsDisordered`. |
| `DisorderedRegions` | `IReadOnlyList<DisorderedRegion>` | Contiguous disordered runs detected from the residue predictions; region semantics are documented in [Disordered_Region_Detection.md](Disordered_Region_Detection.md). |
| `OverallDisorderContent` | `double` | Fraction of residues whose `IsDisordered` flag is `true`. |
| `MeanDisorderScore` | `double` | Arithmetic mean of all `ResiduePrediction.DisorderScore` values. |

### 3.3 Preconditions and Validation

`PredictDisorder(...)` uppercases the sequence with `ToUpperInvariant()` and returns an empty `DisorderPredictionResult` when the input is `null` or empty ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). Characters that are absent from the TOP-IDP lookup table are preserved in `ResiduePrediction.Residue`, but `CalculateDisorderScore(...)` excludes them from the propensity average; if a window contains no recognized residues, the returned score is `0.0` ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). `ResiduePrediction.Position` is 0-based, and the tests verify that mixed-case and lowercase inputs produce the same scores as uppercase input ([DisorderPredictor_DisorderPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs)).

## 4. Algorithm

### 4.1 High-Level Steps

1. Return an empty result immediately when the input sequence is `null` or empty.
2. Uppercase the sequence.
3. For each residue position, construct the local window, normalize each recognized TOP-IDP value, and average the normalized values.
4. Compare the average against the active threshold to populate `ResiduePrediction.IsDisordered`.
5. Aggregate residue predictions into contiguous disordered regions and compute summary statistics for the result object.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Per-residue scoring for sequence length `n` and window size `w` | `O(nw)` | `O(n)` | Each residue recomputes a window average in `CalculatePerResidueScores(...)`. |
| End-to-end `PredictDisorder(...)` result construction | `O(nw) + O(n)` | `O(n)` | Region extraction and summary statistics are linear after residue scoring ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)

- `DisorderPredictor.PredictDisorder(string, int, double, int)`: public entry point returning residue-level predictions, region calls, and summary statistics.
- `DisorderPredictor.CalculatePerResidueScores(string, int, double)`: private helper that creates `ResiduePrediction` records.
- `DisorderPredictor.CalculateDisorderScore(string)`: private helper that computes the normalized TOP-IDP average for one window.
- `DisorderPredictor.GetDisorderPropensity(char)`: public accessor for the raw TOP-IDP table.
- `DisorderPredictor.IsDisorderPromoting(char)`: public check for the Dunker disorder-promoting residue set.
- `DisorderPredictor.CalculateHydropathy(string)`: public Kyte-Doolittle utility for mean hydropathy.

**Supporting tests:** [DisorderPredictor_DisorderPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs), [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)

### 5.2 Current Behavior

The implementation clips edge windows to the available sequence bounds rather than padding them to a fixed width, so terminal residues can be scored from fewer than `windowSize` positions ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). `OverallDisorderContent` is computed as `disorderedCount / sequence.Length`, and `MeanDisorderScore` is the average of all residue-level scores from the returned prediction list. `DisorderedRegions` are computed internally by the private region-detection helpers documented in [Disordered_Region_Detection.md](Disordered_Region_Detection.md). The public residue-class properties expose sorted cached lists for the disorder-promoting, order-promoting, and ambiguous sets defined in the source file.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The TOP-IDP propensity table values documented in Campen et al. (2008).
- Normalization against the documented TOP-IDP minimum (`-0.884`) and maximum (`0.987`) before window averaging.
- Default residue classification at the TOP-IDP cutoff `0.542`.

**Intentionally simplified:**

- This implementation is a single-feature TOP-IDP heuristic without evolutionary profiles or trained model components; **consequence:** its residue scores are useful for coarse screening within this repository but are not competitive with modern disorder predictors named in the source remarks.

**Not implemented:**

- Predictors that incorporate evolutionary profiles, predicted secondary structure, or trained machine-learning models; **users should rely on:** external tools named in the source remarks, such as IUPred2A or MobiDB-lite, rather than this repository's single-feature heuristic ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Non-canonical residues are skipped during TOP-IDP averaging | Assumption | Windows containing unknown symbols are averaged over fewer effective residues; a window with no recognized residues scores `0.0`. | accepted | Implemented in `CalculateDisorderScore(...)`; the tests include an all-unknown input in the score-range checks ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs); [DisorderPredictor_DisorderPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs)). |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `null` or empty sequence | Returns an empty result with zero summary statistics. | Explicit short-circuit in `PredictDisorder(...)` ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |
| Lowercase or mixed-case input | Produces the same scores as uppercase input. | The method uppercases the sequence first; this is asserted in the tests ([DisorderPredictor_DisorderPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs)). |
| Window containing only unknown residues | Produces a disorder score of `0.0` for that window. | `CalculateDisorderScore(...)` returns `0` when no recognized residues are counted ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). |
| `minRegionLength` greater than every disordered run | `DisorderedRegions` is empty even when some residues score above the threshold. | Region emission is filtered after residue scoring; this behavior is exercised in the region tests ([DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)). |

### 6.2 Limitations

The source remarks state that this code is a single-feature heuristic toolkit and not a publication-grade disorder predictor ([DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)). It does not add evolutionary profiles, predicted structural context, or machine-learned features beyond the TOP-IDP composition signal. The public method also performs no explicit validation for non-default `windowSize`, `disorderThreshold`, or `minRegionLength` values.

## 7. Examples and Related Material

### 7.1 Worked Example

The prediction tests use homopolymeric anchor sequences to verify the normalization endpoints and cutoff handling ([DisorderPredictor_DisorderPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs)):

- `new string('W', 30)` produces an interior disorder score of `0.0`, the normalized minimum.
- `new string('P', 30)` produces an interior disorder score of `1.0` and `OverallDisorderContent = 1.0`.
- `new string('E', 30)` produces an interior disorder score of about `0.866`, which is above the default `0.542` cutoff.

## 8. References

1. Campen A., et al. (2008). "TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder." Protein and Peptide Letters 15(9):956-963.
2. Dunker A. K., et al. (2001). "Intrinsically disordered protein." Journal of Molecular Graphics and Modelling 19(1):26-59.
3. Kyte J., Doolittle R. F. (1982). "A simple method for displaying the hydropathic character of a protein." Journal of Molecular Biology 157(1):105-132.
4. Uversky V. N., Gillespie J. R., Fink A. L. (2000). "Why are 'natively unfolded' proteins unstructured under physiologic conditions?" Proteins 41(3):415-427.
5. [DisorderPredictor.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs)
6. [DisorderPredictor_DisorderPrediction_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderPrediction_Tests.cs)
7. [DisorderPredictor_DisorderedRegion_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/DisorderPredictor_DisorderedRegion_Tests.cs)
