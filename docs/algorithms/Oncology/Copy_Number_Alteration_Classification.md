# Copy-Number Alteration Classification

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-CNA-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Given a per-region log2 copy ratio `log2(tumor_depth / normal_depth)`, this algorithm converts it to an absolute copy number and assigns one of five discrete copy-number alteration (CNA) states: deep deletion, loss, neutral, gain, or amplification. It is the oncology classification layer above read-depth/log2 estimation: the integer copy-number call uses CNVkit's hard-threshold method, which is exact and deterministic for a given threshold set [2]. It is specification-driven (fixed cutoffs), not probabilistic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Somatic copy-number alterations are gains and losses of genomic segments in a tumor relative to the matched normal. Read depth is linearly related to copy number, so the depth ratio between tumor and normal, expressed in log2, is the standard quantitative signal [1]. A copy-number-neutral diploid region has tumor depth ≈ normal depth, i.e. log2 ratio ≈ 0; losses are negative, gains positive.

### 2.2 Core Model

For a pure sample with reference (germline) ploidy `p`, the absolute copy number of a region with log2 ratio `v` is

```
n = p · 2^v
```

(CNVkit `_log2_ratio_to_absolute_pure`: `ncopies = ref_copies * 2**log2_ratio`, with `ref_copies = ploidy = 2` for autosomes) [2]. So `v = 0 ⇒ n = 2`, `v = 1 ⇒ n = 4`, `v = −1 ⇒ n = 1`.

Integer copy number is called with hard thresholds. Given four ascending cutoffs `t₀ < t₁ < t₂ < t₃`, the integer copy number `CN` is the index of the first cutoff that `v ≤ tᵢ` (counting from 0); if `v` exceeds all cutoffs, `CN = ⌈p · 2^v⌉` [2]. The default tumor-sample cutoffs are `(−1.1, −0.25, 0.2, 0.7)`, stated verbatim in the CNVkit `absolute_threshold` docstring as `DEL(0) < −1.1`, `LOSS(1) < −0.25`, `GAIN(3) ≥ +0.2`, `AMP(4) ≥ +0.7`, "reasonably 'safe' for a tumor sample with purity of at least 30%" [2][3].

The integer copy number maps to a CNA state: `CN 0 → DeepDeletion`, `CN 1 → Loss`, `CN 2 → Neutral`, `CN 3 → Gain`, `CN ≥ 4 → Amplification`. GISTIC2.0 corroborates a neutral noise band (low-amplitude threshold log2 ± 0.1) with amplification/deletion at high amplitude (0.848 / −0.737) [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Reference ploidy is 2 (autosomal diploid) [2] | A non-diploid baseline (e.g. sex chromosomes, whole-genome doubling) shifts the absolute CN and amplification ceiling; the discrete cutoffs would no longer correspond to the same integer states |
| ASM-02 | Sample is treated as pure for the absolute-CN formula [2] | At low purity the same true CN yields a log2 ratio closer to 0, so impure samples under-call gains/losses (the cutoffs are calibrated for purity ≥ 30%) |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Log2RatioToCopyNumber(0, 2) = 2.0` | `2 · 2^0 = 2` [2] |
| INV-02 | Integer CN is non-decreasing in log2 ratio | cutoffs are ascending and the else-branch `⌈2·2^v⌉` is monotone in `v` [2] |
| INV-03 | Integer CN ≥ 0 for finite log2 | CN 0 is the lowest state; the ceiling of a non-negative quantity is ≥ 0 [2] |
| INV-04 | State ↔ CN: 0→DeepDeletion, 1→Loss, 2→Neutral, 3→Gain, ≥4→Amplification | CNVkit `absolute_threshold` DEL/LOSS/neutral/GAIN/AMP labeling [2] |
| INV-05 | `ClassifyCopyNumbers` output length = input length; element i ↔ input i | per-element deterministic map |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| log2Ratio | double | required | log2(tumor_depth/normal_depth) | any finite value; NaN = no-call |
| thresholds | IReadOnlyList\<double\>? | null → (−1.1, −0.25, 0.2, 0.7) | four hard cutoffs | exactly 4, strictly ascending, no NaN |
| ploidy | double | 2.0 | reference ploidy | > 0 |
| log2Ratios | IEnumerable\<double\> | required | per-region log2 ratios (batch) | not null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Log2Ratio | double | echoed input log2 ratio |
| AbsoluteCopyNumber | double | continuous `n = ploidy·2^log2` (ploidy for a NaN no-call) |
| IntegerCopyNumber | int | hard-threshold integer copy number (≥ 0) |
| State | CopyNumberState | DeepDeletion / Loss / Neutral / Gain / Amplification |

### 3.3 Preconditions and Validation

`thresholds` must be exactly four strictly ascending non-NaN values (else `ArgumentException`); `ploidy` must be positive (else `ArgumentOutOfRangeException`); the batch enumerable must not be null (else `ArgumentNullException`). A NaN log2 ratio is a no-call and returns the neutral reference copy number (rounded ploidy = 2 → Neutral) per CNVkit [2]. Threshold comparison is inclusive (`v ≤ tᵢ`), so a value exactly on a cutoff is assigned the lower state of the bin.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate thresholds (four ascending) and ploidy (> 0).
2. If log2 ratio is NaN, return neutral (rounded ploidy).
3. Otherwise scan cutoffs ascending; return the index of the first cutoff with `log2 ≤ cutoff`.
4. If no cutoff matched, return `⌈ploidy · 2^log2⌉`.
5. Map the integer copy number to a CNA state.

### 4.2 Decision Rules, Scoring, Reference Tables

Default cutoffs (CNVkit `do_call`, source-code defaults) [2]:

| log2 ratio range | Integer CN | CNA state |
|------------------|-----------|-----------|
| log2 ≤ −1.1 | 0 | DeepDeletion |
| −1.1 < log2 ≤ −0.25 | 1 | Loss |
| −0.25 < log2 ≤ 0.2 | 2 | Neutral |
| 0.2 < log2 ≤ 0.7 | 3 | Gain |
| log2 > 0.7 | ⌈2·2^log2⌉ (≥ 4) | Amplification |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| ClassifyCopyNumber | O(1) | O(1) | fixed four-cutoff scan |
| ClassifyCopyNumbers | O(m) | O(m) | m = number of regions |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.Log2RatioToCopyNumber(log2Ratio, ploidy)`: continuous `n = ploidy·2^log2`.
- `OncologyAnalyzer.CallCopyNumber(log2Ratio, thresholds, ploidy)`: hard-threshold integer copy number.
- `OncologyAnalyzer.ClassifyCopyNumber(log2Ratio, thresholds, ploidy)`: full `CopyNumberCall` (absolute, integer, state).
- `OncologyAnalyzer.ClassifyCopyNumbers(log2Ratios, thresholds, ploidy)`: per-region batch (order/length preserving).

### 5.2 Current Behavior

This unit is the oncology classification layer. SV-CNV-001 (`StructuralVariantAnalyzer.DetectCNV` / `SegmentCopyNumber`) already converts read depth → log2 → integer CN via `round(2·2^log2)` and merges segments; it does not produce the five discrete CNA states. ONCO-CNA-001 reuses the same `n = 2·2^log2` conversion formula (cited to the same CNVkit source) but uses CNVkit's hard-threshold binning (`absolute_threshold`) rather than nearest-integer rounding, and adds the state classification. No suffix tree is involved (no substring search) — see §5.3.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `n = ploidy · 2^log2` absolute copy number (CNVkit `_log2_ratio_to_absolute_pure`) [2].
- Hard-threshold integer calling: first cutoff with `log2 ≤ cutoff`, else `⌈ploidy·2^log2⌉` (CNVkit `absolute_threshold`) [2].
- Default tumor cutoffs (−1.1, −0.25, 0.2, 0.7) and the 0/1/2/3/≥4 → DeepDeletion/Loss/Neutral/Gain/Amplification mapping [2].
- NaN log2 → neutral no-call [2].

**Intentionally simplified:**

- Purity/impurity correction (`_log2_ratio_to_absolute`): only the pure-sample formula is implemented; **consequence:** at low purity the same true CN gives a log2 closer to 0 and may be under-called (cutoffs assume purity ≥ 30%) [2][3].

**Not implemented:**

- Segmentation (CBS) and allele-specific (major/minor) copy number; **users should rely on:** `StructuralVariantAnalyzer.SegmentCopyNumber` for segmentation and ONCO-LOH-001 for allele-specific LOH.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Default thresholds use CNVkit source-code tumor defaults (−0.25 / 0.2), not the docs-page germline variant (−0.4 / 0.3) | Assumption | shifts loss/gain bin edges by ≤ 0.15 log2 | accepted | callers can override via `thresholds`; ASM-01 |
| 2 | Reference ploidy fixed unless overridden | Assumption | non-diploid baselines mis-call | accepted | ASM-01; `ploidy` parameter exposed |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| log2 exactly on a cutoff | lower state of the bin | inclusive `log2 ≤ cutoff` [2] |
| log2 = NaN | Neutral, CN 2 | no-call → neutral [2] |
| very high log2 (e.g. 2.0) | Amplification, CN = ⌈2·2^2⌉ = 8 | ceiling above last cutoff [2] |
| thresholds null | default (−1.1, −0.25, 0.2, 0.7) | documented default [2] |
| empty batch | empty result | per-element map |

### 6.2 Limitations

Single-region classification only (no segmentation/joining); pure-sample absolute CN (no purity correction); diploid reference unless overridden; no allele-specific (major/minor) decomposition; cutoffs calibrated for tumor purity ≥ 30% [3].

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var call = OncologyAnalyzer.ClassifyCopyNumber(1.0);
// call.AbsoluteCopyNumber == 4.0, call.IntegerCopyNumber == 4,
// call.State == CopyNumberState.Amplification
```

**Numerical walk-through:** log2 = 1.0. Scan cutoffs (−1.1, −0.25, 0.2, 0.7): 1.0 > all four, so CN = ⌈2·2^1.0⌉ = ⌈4.0⌉ = 4 ⇒ Amplification. For log2 = 0: 0 ≤ 0.2 (third cutoff, index 2) ⇒ CN 2 ⇒ Neutral; absolute n = 2·2^0 = 2.0.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_CopyNumberClassification_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_CopyNumberClassification_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [ONCO-CNA-001-Evidence.md](../../../docs/Evidence/ONCO-CNA-001-Evidence.md)

## 8. References

1. Mermel CH, Schumacher SE, Hill B, Meyerson ML, Beroukhim R, Getz G. 2011. GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biology 12(4):R41. https://doi.org/10.1186/gb-2011-12-4-r41
2. CNVkit. `cnvlib/call.py` — `absolute_threshold`, `_log2_ratio_to_absolute_pure`, `do_call`. https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py
3. CNVkit documentation. `call` command threshold method. https://cnvkit.readthedocs.io/en/stable/pipeline.html
4. GISTIC2 documentation. `-ta` / `-td` amplification/deletion thresholds. https://broadinstitute.github.io/gistic2/
