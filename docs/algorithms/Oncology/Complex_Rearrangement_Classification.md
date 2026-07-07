# Somatic Complex Rearrangement Classification (Chromothripsis Inference)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Somatic Structural Variation |
| Test Unit ID | ONCO-SV-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-15 |

## 1. Overview

This algorithm classifies a chromosomal region's somatic copy-number / structural-variant profile as
**chromothripsis** (a one-off catastrophic shattering-and-repair event) or not, using the inference criteria
of Korbel & Campbell (2013) [1] together with the operational thresholds of Cortés-Ciriano et al. (2020) [2]
and the first-pass screen of Magrangeas et al. (2011) [3]. It is the oncology-specific complex-rearrangement
layer that sits above generic SV detection: generic DEL/DUP/INV/translocation typing from breakpoints is
provided elsewhere (`StructuralVariantAnalyzer`); this unit recognises the *pattern* of clustered breakpoints
with oscillating copy number that distinguishes chromothripsis from progressive amplification. The decision is
a deterministic, rule-based screen, not a probabilistic clinical caller.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Chromothripsis is localized chromosome shattering and repair in a single catastrophic event, producing tens to
hundreds of clustered rearrangements on one or a few chromosomes [1]. Korbel & Campbell (2013) proposed six
hallmark criteria for inferring it from cancer genome data: (A) clustering of breakpoints; (B) regularity of
oscillating copy-number states; (C) interspersed retention and loss of heterozygosity; (D) prevalence of
rearrangements affecting a specific haplotype; (E) randomness of DNA fragment order and fragment joins; and
(F) the ability to "walk" the derivative chromosome (invariant alternation between head and tail sequences) [1].
The signature that most cleanly separates chromothripsis from gradual amplification (e.g. breakage-fusion-bridge)
is criterion B: the copy number oscillates between (canonically) **two** states rather than climbing through many
ascending levels [1].

### 2.2 Core Model

- **Oscillating copy-number changes (criterion B).** Walking the per-segment integer copy numbers in genomic
  order, an *oscillation* is a segment whose copy-number state differs from the immediately preceding segment.
  The first-pass chromothripsis screen requires at least **10** oscillating copy-number changes [3] (the lowest
  of the "10, 20, or 50" first-pass cutoffs cited in the Korbel & Campbell framework [4]). The profile must
  oscillate between at most a small number of distinct states (canonically two, at most two-or-three) [1]; a
  monotone or many-state ascending profile is progressive amplification, not chromothripsis.
- **Structural-variant burden.** Cortés-Ciriano et al. (2020) exclude focal events "comprising fewer than six
  SVs" [2], so at least **6** clustered intrachromosomal SVs are required.
- **Confidence tier (segment count).** High-confidence calls "display oscillations between two states in at
  least seven adjacent segments"; low-confidence calls "involve between four and six segments" [2]. The tier is
  derived from the count of adjacent oscillating segments (a run of *k* transitions spans *k+1* segments).
- **Breakpoint clustering (criterion A).** Under the null of uniformly random breakpoints, inter-breakpoint
  distances are exponentially distributed [1][4]. The exponential distribution has coefficient of variation
  (CV = sd/mean) equal to 1; over-dispersion toward many short gaps with a few long gaps (a tight cluster plus
  outliers) gives CV > 1, flagging clustering.

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | An "oscillating copy-number change" is operationalised as an adjacent per-segment CN-state transition, with the profile bounded to ≤ 3 distinct states. | If a tool counts oscillations differently (e.g. only full down-up cycles), the screen count differs; calls near the threshold could change. |
| ASM-02 | Breakpoint clustering is summarised by the coefficient of variation of inter-breakpoint distances against the exponential null (CV = 1). | A different goodness-of-fit statistic (KS, χ²) could flag borderline cases differently; CV is a transparent proxy, not a clinical caller. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Oscillation count ∈ [0, n−1] for n segments. | One transition per adjacent pair; n segments have n−1 pairs. |
| INV-02 | Chromothripsis ⇒ distinct states ∈ [2, 3] AND oscillations ≥ 10 AND SV burden ≥ 6. | Conjunction of criterion B (two-state) [1], the ≥10 first-pass screen [3], and the ≥6 SV floor [2]. |
| INV-03 | Confidence = High iff oscillating segments ≥ 7; Low iff in [4, 6]; None iff < 4. | Cortés-Ciriano 2020 segment thresholds [2]. |
| INV-04 | A monotone profile with > 3 distinct states is never Chromothripsis, regardless of length. | Two-state hallmark, criterion B [1]; distinct-state cap in INV-02. |
| INV-05 | Regular-spacing breakpoints give CV ≈ 0 (not clustered); over-dispersed give CV > 1 (clustered). | CV of equal gaps is 0; exponential null CV = 1 [1][4]. |

### 2.5 Comparison with Related Methods

| Aspect | Chromothripsis (this layer) | Generic SV typing (`StructuralVariantAnalyzer`) |
|--------|-----------------------------|--------------------------------------------------|
| Output | Pattern class (chromothripsis vs not) + confidence | Individual SV type (DEL/DUP/INV/TRA) per breakpoint |
| Key signal | Oscillating CN between ≤3 states + clustered breakpoints | Read-pair / split-read breakpoint orientation |
| Scope | Oncology-specific complex-rearrangement recognition | Per-event detection/classification |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| segmentCopyNumbers | IReadOnlyList&lt;int&gt; | required | Per-segment integer copy numbers in genomic order | non-null; may be empty |
| StructuralVariantCount | int | required | Clustered intrachromosomal SV burden | ≥ 0 |
| breakpointPositions | IReadOnlyList&lt;long&gt; | required | Genomic breakpoint coordinates (any order) | non-null; may be empty |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Type | ComplexRearrangementType | NotComplex or Chromothripsis |
| Confidence | ChromothripsisConfidence | None / Low / High from oscillating-segment count |
| OscillationCount | int | Adjacent CN state transitions |
| OscillatingSegmentCount | int | Segments participating in the oscillation (k transitions → k+1) |
| DistinctStateCount | int | Number of distinct CN states |
| StructuralVariantCount | int | Echo of the input SV burden |
| BreakpointClusteringResult | record | BreakpointCount, MeanGap, CoefficientOfVariation, IsClustered |

### 3.3 Preconditions and Validation

`segmentCopyNumbers` and `breakpointPositions` must be non-null (`ArgumentNullException` otherwise). Empty or
single-element copy-number lists yield 0 oscillations and `NotComplex`. Breakpoint sets with fewer than three
positions cannot define a CV (fewer than two gaps) and return `IsClustered = false`. Copy numbers are 0-based
integer states; breakpoints are genomic coordinates and are sorted internally. No randomness is used; output is
fully deterministic.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count adjacent CN-state transitions along the segment profile (oscillations).
2. Derive oscillating-segment count (transitions + 1 when > 0) and distinct-state count.
3. Assign confidence tier from oscillating-segment count (≥7 High; 4–6 Low; else None).
4. Call Chromothripsis iff distinct states ∈ [2,3] AND oscillations ≥ 10 AND SV burden ≥ 6.
5. (Clustering) Sort breakpoints, compute inter-breakpoint gaps, CV = sd/mean; flag clustered when CV > 1.

### 4.2 Decision Rules, Scoring, Reference Tables

| Constant | Value | Source |
|----------|-------|--------|
| MinOscillatingCopyNumberChanges | 10 | Magrangeas et al. 2011 first-pass screen [3][4] |
| MaxChromothripsisCopyNumberStates | 3 | Korbel & Campbell 2013 (two-or-three state hallmark) [1] |
| MinChromothripsisSvBurden | 6 | Cortés-Ciriano et al. 2020 (focal <6 SV exclusion) [2] |
| HighConfidenceOscillatingSegments | 7 | Cortés-Ciriano et al. 2020 (≥7 adjacent segments) [2] |
| LowConfidenceOscillatingSegments | 4 | Cortés-Ciriano et al. 2020 (4–6 segments) [2] |
| Exponential-null CV | 1.0 | Korbel & Campbell 2013 (exponential inter-breakpoint distances) [1][4] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CountCopyNumberStateOscillations | O(n) | O(1) | n segments |
| ClassifyComplexRearrangement | O(n) | O(n) | distinct-state set |
| TestBreakpointClustering | O(m log m) | O(m) | m breakpoints; dominated by the sort |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CountCopyNumberStateOscillations(...)`: counts adjacent CN-state transitions.
- `OncologyAnalyzer.TestBreakpointClustering(...)`: exponential-null CV clustering summary.
- `OncologyAnalyzer.ClassifyComplexRearrangement(...)`: full chromothripsis classification with confidence tier.

### 5.2 Current Behavior

The classifier is a deterministic rule screen. It does not call SVs from raw reads (that is
`StructuralVariantAnalyzer`'s role); it consumes already-segmented per-region copy numbers and a clustered SV
count. Criteria C/D/E/F of Korbel & Campbell (LOH interspersion, haplotype prevalence, fragment-join randomness,
derivative-chromosome walk) are not evaluated; the implemented gate uses criteria A and B plus the
Cortés-Ciriano SV-burden floor. **Search reuse:** the repository suffix tree is not applicable — there is no
substring/pattern-matching step here (the work is numeric counting and a statistical summary), so it is not used.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Oscillating-CN-change count and the ≥10 first-pass screen (Magrangeas 2011 / Korbel & Campbell framework) [3][4].
- Two-/three-state oscillation hallmark distinguishing chromothripsis from progressive amplification (criterion B) [1].
- ≥6 clustered intrachromosomal SV burden floor (Cortés-Ciriano 2020) [2].
- High/Low/None confidence tiers from adjacent oscillating-segment counts (≥7 / 4–6 / <4) (Cortés-Ciriano 2020) [2].
- Breakpoint-clustering test against the exponential null via the coefficient of variation (criterion A) [1][4].

**Intentionally simplified:**

- Breakpoint clustering: summarised by CV vs the exponential null rather than a formal goodness-of-fit test;
  **consequence:** borderline clustering may be flagged differently than a KS/χ² test would.
- Oscillation = adjacent state transition (not full down-up cycle); **consequence:** the screen count is the
  transition count, which is the directly source-supported quantity.

**Not implemented:**

- Korbel & Campbell criteria C (interspersed LOH), D (haplotype prevalence), E (fragment-join randomness), and
  F (derivative-chromosome walk); **users should rely on:** haplotype-aware / breakpoint-sequence tooling
  (no current in-repo alternative for these criteria).
- Other complex patterns (chromoplexy, BFB scoring); **users should rely on:** no current alternative.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Oscillation = adjacent transition | Assumption | Affects screen count near threshold | accepted | ASM-01 |
| 2 | Clustering via CV vs exponential null | Assumption | Borderline clustering flagging | accepted | ASM-02 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty / single segment | 0 oscillations, NotComplex | INV-01 boundary |
| Monotone rising/falling CN (>3 states) | NotComplex regardless of length | Two-state hallmark, INV-04 [1] |
| 5 oscillations (below screen) | NotComplex; confidence Low if ≥4 segments | ≥10 screen [3]; tier [2] |
| SV burden < 6 | NotComplex even if oscillations ≥ 10 | Focal exclusion [2] |
| < 3 breakpoints | IsClustered = false | CV undefined with < 2 gaps |
| null inputs | ArgumentNullException | API contract |

### 6.2 Limitations

Only criteria A and B of Korbel & Campbell, plus the Cortés-Ciriano SV-burden floor, are evaluated. The screen
does not establish chromothripsis on its own (clustering is "necessary but not sufficient" [1][4]); it is a
deterministic pattern recogniser, not a validated clinical diagnostic. Multi-chromosome events, chromoplexy, and
breakage-fusion-bridge are out of scope.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var input = new OncologyAnalyzer.ComplexRearrangementInput(
    SegmentCopyNumbers: new[] { 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2 }, // 11 segments, 10 oscillations
    StructuralVariantCount: 12);
var result = OncologyAnalyzer.ClassifyComplexRearrangement(input);
// result.Type == Chromothripsis; result.Confidence == High; result.OscillationCount == 10; DistinctStateCount == 2
```

**Numerical walk-through:** profile 2,1,2,1,2,1,2,1,2,1,2 → 10 transitions; 2 distinct states {1,2}; 11
oscillating segments ≥ 7 → High; 10 ≥ 10 and 12 ≥ 6 and states ∈ [2,3] → Chromothripsis.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_ClassifyComplexRearrangement_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_ClassifyComplexRearrangement_Tests.cs) — covers INV-01..INV-05
- Evidence: [ONCO-SV-001-Evidence.md](../../../docs/Evidence/ONCO-SV-001-Evidence.md)

## 8. References

1. Korbel JO, Campbell PJ. 2013. Criteria for Inference of Chromothripsis in Cancer Genomes. Cell 152(6):1226–1236. https://doi.org/10.1016/j.cell.2013.02.023
2. Cortés-Ciriano I, Lee JJ-K, Xi R, et al. 2020. Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing. Nature Genetics 52:331–341. https://doi.org/10.1038/s41588-019-0576-7
3. Magrangeas F, Avet-Loiseau H, Munshi NC, Minvielle S. 2011. Chromothripsis identifies a rare and aggressive entity among newly diagnosed multiple myeloma patients. Blood 118(3):675–678. https://doi.org/10.1182/blood-2011-03-344069
4. Maher CA, Wilson RK. 2012. Chromothripsis and human disease (review enumerating the Korbel & Campbell criteria). https://pmc.ncbi.nlm.nih.gov/articles/PMC3861665/
