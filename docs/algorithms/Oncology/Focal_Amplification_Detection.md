# Focal Amplification Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-CNA-002 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Focal amplification detection separates highly amplified, *focal* copy-number segments from *broad*
(arm-level) amplifications and maps the focal events to recurrently amplified oncogenes. It is
specification-driven: GISTIC2.0 classifies a copy-number event as focal or arm-level purely by its
length relative to its chromosome arm, and an event is "amplified" when its log2 gain exceeds an
amplitude threshold [1][2]. Focal high-amplitude amplifications are the therapeutically actionable
targets, so the algorithm filters segments to that subset and reports the panel oncogenes resident on
their arms [3][4]. It is a deterministic O(n) filter, not a probabilistic peak caller.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Somatic copy-number alterations (SCNAs) in cancer occur at two scales: broad events that span a whole
chromosome arm (or more), and focal events confined to a sub-arm region. GISTIC2.0 observed that events
occupying exactly one chromosome arm are so frequent that arm-level and focal alterations must be modeled
separately; their lengths form a reproducible bimodal distribution that gives "a natural basis for
classifying events as 'arm-level' and 'focal' based purely on length" [1].

### 2.2 Core Model

For a segment of length L on a chromosome arm of length A with mean log2 copy ratio r:

- **Amplified** ⇔ r > t_amp, with t_amp = 0.1 (GISTIC2 `t_amp`) [2].
- **Focal** ⇔ (L / A) < c, with c = 0.98 (GISTIC2 `broad_len_cutoff`); an event occupying ≥ 98% of an
  arm is arm-level [1][2].
- **Focal amplification** ⇔ Amplified ∧ Focal.

A single-copy gain has copy ratio 3/2, i.e. log2(3/2) = 0.585, which is far above t_amp = 0.1, so the
amplitude test admits any genuine gain while rejecting low-level artifactual segments [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Every reported amplification satisfies L/A < 0.98 (strict) | focal test is `ArmFraction < BroadLengthCutoff` [1][2] |
| INV-02 | Every reported amplification satisfies r > t_amp | amplitude test is `Log2Ratio > AmplificationLog2Threshold` [2] |
| INV-03 | Output is a subset of the input in input order | the detector is a filter; it constructs no new segments |
| INV-04 | An oncogene is reported only for an arm carrying a focal amplification | the mapper consumes focal amplifications only [4] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| segments | `IEnumerable<CopyNumberArmSegment>` | required | Arm-anchored copy-number segments | not null; each ArmLength > 0, End > Start |
| thresholds | `FocalAmplificationThresholds?` | GISTIC2 defaults | Amplitude + length cutoffs | null ⇒ t_amp 0.1, broad_len_cutoff 0.98 |
| amplifications | `IEnumerable<CopyNumberArmSegment>` | required | Focal amplifications to map to oncogenes | not null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (DetectFocalAmplifications) | `IReadOnlyList<CopyNumberArmSegment>` | Input segments that are focal amplifications, in input order |
| (IdentifyAmplifiedOncogenes) | `IReadOnlyList<string>` | Distinct panel oncogene symbols on amplified arms, in panel order |

### 3.3 Preconditions and Validation

Null `segments`/`amplifications` ⇒ `ArgumentNullException`. A segment with non-positive `ArmLength` or
with `End ≤ Start` ⇒ `ArgumentException`. Arm labels are matched case-insensitively (Ordinal-ignore-case).
Coordinates are base-pair counts; segment length is `End − Start`. The boundary fraction exactly equal to
the cutoff (0.98) is arm-level (the focal test is strictly less-than).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate each segment (positive arm length, End > Start).
2. For each segment compute amplitude test `r > t_amp` and focal test `L/A < broad_len_cutoff`.
3. Emit the segment iff both tests pass, preserving input order.
4. For oncogene mapping, collect the set of arms carrying a focal amplification, then emit each panel
   oncogene whose arm is in that set, in panel order.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Parameter / Table | Value | Source |
|-------------------|-------|--------|
| t_amp (amplitude) | 0.1 (log2) | GISTIC2 `t_amp` [2] |
| broad_len_cutoff (focal/broad) | 0.98 (fraction of arm) | Mermel 2011; GISTIC2 `broad_len_cutoff` [1][2] |
| ERBB2 | 17q | NCBI Gene 2064 (17q12) [4] |
| MYC | 8q | NCBI Gene 4609 (8q24.21) [4] |
| EGFR | 7p | NCBI Gene 1956 (7p11.2) [4] |
| CCND1 | 11q | NCBI Gene 595 (11q13.3) [4] |
| MDM2 | 12q | NCBI Gene 4193 (12q15) [4] |
| CDK4 | 12q | NCBI Gene 1019 (12q14.1) [4] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| DetectFocalAmplifications | O(n) | O(k) | n segments, k focal amplifications |
| IdentifyAmplifiedOncogenes | O(n + g) | O(n) | g = fixed panel size (6) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.DetectFocalAmplifications(segments, thresholds?)`: filters segments to focal amplifications.
- `OncologyAnalyzer.IdentifyAmplifiedOncogenes(amplifications)`: maps focal amplifications to panel oncogenes.
- `OncologyAnalyzer.IsFocalAmplification(segment, thresholds)`: single-segment predicate (internal helper, public for reuse).

### 5.2 Current Behavior

Detection is a single-pass filter; no segmentation is performed here (segmentation is upstream in
`StructuralVariantAnalyzer.SegmentCopyNumber`, SV-CNV-001 / ONCO-CNA-001). No substring/pattern search is
involved, so the repository suffix tree is **not applicable** to this unit. Arm matching uses
case-insensitive ordinal comparison; arms outside the six-gene panel map to no oncogene.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Length-based focal/arm-level split at 98% of chromosome arm (Mermel 2011; GISTIC2 `broad_len_cutoff` 0.98) [1][2].
- Amplitude threshold t_amp = 0.1 for calling a gain amplified (GISTIC2 `t_amp`) [2].
- Oncogene→arm panel from NCBI Gene cytogenetic locations [4].

**Intentionally simplified:**

- Arm boundaries / arm length are supplied by the caller rather than derived from a bundled cytoband
  file; **consequence:** the caller must provide each arm's length (GISTIC2 reads this from the genome
  assembly). The 0.98 rule and amplitude test are unchanged.
- Oncogene mapping is arm-level (any focal amplification on the arm flags the gene) rather than
  coordinate-overlap of the gene locus; **consequence:** a focal amplification elsewhere on the same arm
  also flags the gene. This matches the registry panel's arm-level intent.

**Not implemented:**

- GISTIC2's probabilistic peak/q-value boundary estimation and background-rate modeling; **users should
  rely on:** the full GISTIC2 tool for genome-wide significance peaks. This unit implements only the
  deterministic focal/broad length rule and the oncogene panel mapping.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Amplitude test combined with length rule | Assumption | Defines which gains count as amplifications | accepted | t_amp from GISTIC2 docs (0.1); length rule from paper |
| 2 | Caller-supplied arm length | Assumption | Caller must provide cytoband-derived arm length | accepted | No bundled cytoband table |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Segment = 0.98 of arm | Not focal (arm-level) | Focal test is strict `< 0.98` [1][2] |
| Segment > 98% of arm | Not focal | Arm-level by GISTIC2 [1] |
| log2 ≤ t_amp | Not amplified | GISTIC2 `t_amp` [2] |
| Null input | ArgumentNullException | Guard |
| Empty input | Empty result | Guard |
| ArmLength ≤ 0 or End ≤ Start | ArgumentException | Validation |

### 6.2 Limitations

No significance testing, no background-rate modeling, and no sub-arm peak localization (these are
GISTIC2's probabilistic stages). Oncogene mapping is restricted to the six-gene registry panel and is
arm-level, not gene-locus-overlap. Deletions are out of scope (ONCO-CNA-003).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var segments = new[]
{
    new OncologyAnalyzer.CopyNumberArmSegment("17q", 100_000, 600_000, 1_000_000, 1.0), // focal amp
    new OncologyAnalyzer.CopyNumberArmSegment("8q", 0, 990_000, 1_000_000, 1.5),        // arm-level
};
var focal = OncologyAnalyzer.DetectFocalAmplifications(segments); // [17q segment]
var genes = OncologyAnalyzer.IdentifyAmplifiedOncogenes(focal);   // ["ERBB2"]
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_DetectFocalAmplifications_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectFocalAmplifications_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ONCO-CNA-002-Evidence.md](../../../docs/Evidence/ONCO-CNA-002-Evidence.md)
- Related algorithms: [Copy_Number_Alteration_Classification](./Copy_Number_Alteration_Classification.md)

## 8. References

1. Mermel CH, Schumacher SE, Hill B, Meyerson ML, Beroukhim R, Getz G. 2011. GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biology 12:R41. https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/
2. Broad Institute. GISTIC2 documentation (`broad_len_cutoff`, `t_amp`, `t_del`). https://broadinstitute.github.io/gistic2/
3. Talevich E, Shain AH, Botton T, Bastian BC. CNVkit — Calling copy number gains and losses. https://cnvkit.readthedocs.io/en/stable/calling.html
4. NCBI Gene: ERBB2 (2064), MYC (4609), EGFR (1956), CCND1 (595), MDM2 (4193), CDK4 (1019). https://www.ncbi.nlm.nih.gov/gene/
