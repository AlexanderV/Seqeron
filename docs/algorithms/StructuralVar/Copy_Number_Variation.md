# Read-Depth Copy Number Variation Detection

| Field | Value |
|-------|-------|
| Algorithm Group | StructuralVar |
| Test Unit ID | SV-CNV-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Copy number variation (CNV) detection from read depth of coverage estimates the copy number of a
genomic region from how many sequencing reads map to it: under the read-depth model the mean depth of
a region is proportional to its copy number [1]. The algorithm summarises a per-position depth track
into non-overlapping fixed-size windows, converts each window's mean depth to a log2 ratio against a
copy-number-neutral reference depth, and converts that ratio to an integer copy number for a diploid
genome [2][3]. It is a deterministic, specification-driven numeric transform; it is not a statistical
change-point/segmentation test (e.g. circular binary segmentation), and it does not model tumour
purity or GC bias beyond the median baseline.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Massively-parallel sequencing maps reads to a reference. Regions present in more copies in the sample
recruit proportionally more reads, so local read depth (RD) is a quantitative proxy for copy number.
Yoon et al. demonstrated "a linear relationship between coverage and copy number" and used "the
GC-adjusted RD within [fixed-size] windows as a quantitative measurement of genome copy number" [1].
Copy number is then expressed relative to a reference using the log2 ratio, the standard scale for CNV
callers [2][3].

### 2.2 Core Model

Let a window have mean read depth `RD_w` and let the copy-number-neutral reference depth be `RD_ref`.

- **Read depth per window** — "RD was measured by counting the number of mapped reads in [fixed-size]
  windows" [1]; here `RD_w` is the mean of the per-position depths inside the window.
- **log2 ratio** — `log2ratio = log2(RD_w / RD_ref)`. CNVkit defines this against a reference profile:
  `log2_ratio = log2(ncopies / ploidy)` ⇒ `ncopies = ploidy · 2^log2_ratio` [2].
- **Integer copy number** — `CN = round(ploidy · 2^log2ratio)`. CNVkit's pure/round method computes
  `n = r · 2^v` (`r` = reference copies) and rounds to the nearest integer; for a diploid genome
  "the absolute copy number is calculated as 2 * 2^(log2 value)" [2][3].

Reference log2 anchors (diploid), quoted from CNVkit [2][3]:

| Copy number | Ratio | log2 ratio | Interpretation |
|-------------|-------|------------|----------------|
| 1 | 1/2 | log2(1/2) = −1.0 | single-copy loss (deletion) |
| 2 | 2/2 | log2(2/2) = 0.0 | neutral |
| 3 | 3/2 | log2(3/2) = 0.585 | single-copy gain (duplication) |
| 4 | 4/2 | log2(4/2) = 1.0 | amplification |

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Reference depth = overall median of non-zero window means when no baseline is given (Yoon overall median `m` [1]). | The log2=0 anchor shifts; every copy-number call is biased up or down. |
| ASM-02 | Diploid ploidy (2) is the copy-number baseline (CNVkit diploid conversion [2][3]). | Calls are off by the ratio of true to assumed ploidy (e.g. triploid samples mis-scaled). |
| ASM-03 | Read depth is proportional to copy number after mapping (no uncorrected GC/mappability bias) [1]. | Local depth biases masquerade as copy-number change. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `RD_w = RD_ref` ⇒ log2 = 0 ⇒ CN = ploidy (2). | `round(2 · 2^0) = 2` [2]. |
| INV-02 | `CN = round(ploidy · 2^log2)`, `log2 = log2(RD_w / RD_ref)`. | CNVkit `_log2_ratio_to_absolute_pure` `n = r·2^v` [2]. |
| INV-03 | `CN ≥ 0` (non-negative integer). | CNVkit clamps `max(0, n)` [2]. |
| INV-04 | CN is monotonically non-decreasing in `RD_w`. | `2·2^log2` is strictly increasing in `RD_w` (log2 increasing) [1][2]. |
| INV-05 | Depth is partitioned into non-overlapping windows of the given size; a trailing partial window is dropped. | Fixed-size read-count windows [1]. |
| INV-06 | A zero-depth window is a no-call (excluded), not a finite call. | `log2(0)` is undefined; CNVkit treats unusable signal as a no-call [1][2]. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `depthData` | `IReadOnlyList<int>` | required | Per-position read depth along one contig, 0-based. | Non-null; values ≥ 0. |
| `windowSize` | `int` | 100 | Positions per non-overlapping window. | > 0 (else `ArgumentOutOfRangeException`). |
| `referenceDepth` | `double?` | null | Neutral reference depth (log2 anchor); null ⇒ overall median of non-zero window means. | If given, used as-is; ≤ 0 reference ⇒ no calls. |
| `chromosome` | `string` | "chr1" | Contig label recorded on each segment. | — |
| `logRatios` (`SegmentCopyNumber`) | `IEnumerable<double>` | required | Per-window log2 ratios; window i covers coordinate i. | Non-null; NaN entries are no-calls. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `CopyNumberSegment.Start` / `.End` | `int` | Inclusive 0-based window bounds (`DetectCNV`) or window indices (`SegmentCopyNumber`). |
| `CopyNumberSegment.LogRatio` | `double` | log2(RD_w / RD_ref) (per window) or mean log2 over a merged run. |
| `CopyNumberSegment.CopyNumber` | `int` | `round(2 · 2^LogRatio)`, ≥ 0. |
| `CopyNumberSegment.BAlleleFrequency` | `double` | `NaN` — BAF needs allele-specific data not provided to a depth-only caller. |
| `CopyNumberSegment.ProbeCount` | `int` | Positions in the window (`DetectCNV`) or windows merged (`SegmentCopyNumber`). |

### 3.3 Preconditions and Validation

`depthData` / `logRatios` null ⇒ `ArgumentNullException`; `windowSize ≤ 0` ⇒ `ArgumentOutOfRangeException`.
Coordinates are 0-based with inclusive end. `DetectCNV` drops a trailing partial window (INV-05) and
excludes zero-depth windows as no-calls (INV-06); `SegmentCopyNumber` drops NaN log2 ratios as no-calls.
When the computed reference depth is ≤ 0 (all windows zero), no segments are emitted.

## 4. Algorithm

### 4.1 High-Level Steps

1. Partition `depthData` into `floor(n / windowSize)` non-overlapping full windows; compute each
   window's mean read depth (INV-05).
2. Determine the reference depth: the supplied `referenceDepth`, or the overall median of the non-zero
   window means (ASM-01).
3. For each window with positive mean depth, compute `log2 = log2(RD_w / RD_ref)` and
   `CN = round(2 · 2^log2)` (INV-02); emit a `CopyNumberSegment`. Skip zero-depth windows (INV-06).
4. (`SegmentCopyNumber`) Convert each log2 ratio to CN with the same rule and merge maximal runs of
   consecutive windows of equal CN into one segment (mean log2, window count as probe count).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Copy-number call rule (diploid): `CN = max(0, round(2 · 2^log2))` [2][3]. Reference baseline: overall
median of non-zero window means [1]. No external scoring tables are required; the only numeric
constants are the diploid ploidy (2) and the literature default window size (100) [1].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DetectCNV` | O(n) | O(n / windowSize) | n = depth positions; one pass to sum windows + O(w log w) median over w = windowCount windows (w ≪ n). |
| `SegmentCopyNumber` | O(m) | O(s) | m = log2 ratios; s = longest merged run held for averaging. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [StructuralVariantAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs)

- `StructuralVariantAnalyzer.DetectCNV(depthData, windowSize, referenceDepth, chromosome)`: windowed
  read-depth → log2 ratio → integer copy number; one segment per full non-zero window (canonical).
- `StructuralVariantAnalyzer.SegmentCopyNumber(logRatios, chromosome)`: log2-ratio → CN with adjacent
  equal-CN run merging (segmentation variant).
- `LogRatioToCopyNumber` (private): `max(0, round(2 · 2^log2))`.
- `OverallMedianNonZero` (private): overall median baseline of non-zero window means.

### 5.2 Current Behavior

The reference depth, when not supplied, is the overall median of non-zero window means (ASM-01); an
explicit `referenceDepth` overrides it. Zero-depth windows are excluded (no-call, INV-06) rather than
emitting a `−∞` log2; the existing `CopyNumberSegment.BAlleleFrequency` field is set to `NaN` because a
depth-only caller has no allele-specific data. This unit performs no substring/pattern search, so the
repository suffix tree is not applicable; the work is windowed numeric aggregation over a depth array.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Windowed read-depth summary as a copy-number proxy (Yoon et al. 2009 [1]).
- log2 ratio `log2(RD_w / RD_ref)` and integer copy number `round(ploidy · 2^log2)` for a diploid
  genome (CNVkit `_log2_ratio_to_absolute_pure` / calling docs [2][3]).
- Non-negative copy number clamp and overall-median baseline (Yoon `m` [1]; CNVkit `max(0,n)` [2]).

**Intentionally simplified:**

- GC-content / mappability correction: omitted; **consequence:** windows with extreme GC are not
  bias-corrected (Yoon `r' = r·m/m_GC` [1] is not applied), so calls in such regions may be biased.
- Statistical segmentation: `DetectCNV` reports per-window calls and `SegmentCopyNumber` merges only
  exact-equal-CN runs; **consequence:** no significance test / circular binary segmentation, so noisy
  windows are not smoothed the way EWT (Yoon [1]) or CBS would.

**Not implemented:**

- Tumour-purity / heterogeneity rescaling (CNVkit impure-sample formula `n = (r·2^v − x(1−p))/p` [2]);
  **users should rely on:** a dedicated somatic CNV caller (e.g. CNVkit with `--purity`).
- B-allele-frequency / allele-specific copy number; **users should rely on:** an SNP-aware caller; the
  `BAlleleFrequency` field is reported as `NaN` here.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Reference = overall median of non-zero window means | Assumption | Sets the log2=0 anchor when no baseline is given | accepted | ASM-01; overridable via `referenceDepth` |
| 2 | Diploid ploidy = 2 | Assumption | Scales all copy-number calls | accepted | ASM-02; standard human autosomal baseline [2][3] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty `depthData` | Empty result | No window can be formed. |
| Null `depthData` / `logRatios` | `ArgumentNullException` | Input-validation contract. |
| `windowSize ≤ 0` | `ArgumentOutOfRangeException` | A window must have positive size. |
| `windowSize` > length | Empty result | No full window (INV-05). |
| Zero-depth window | Excluded (no-call) | `log2(0)` undefined (INV-06) [1][2]. |
| Window mean = reference | CN = 2 (neutral) | INV-01 [2]. |

### 6.2 Limitations

Depth-only: no GC/mappability correction, no purity/ploidy estimation, no allele-specific signal, and
no statistical change-point detection. Results are unreliable in low-mappability or extreme-GC regions
and for non-diploid or impure (tumour) samples; for those, use a dedicated caller (CNVkit/EWT). Output
is undefined-as-no-call where read depth is zero.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// Reference depth 100; windows of 4 positions.
int[] depth = { 100,100,100,100,  50,50,50,50,  150,150,150,150 };
var segments = StructuralVariantAnalyzer.DetectCNV(depth, windowSize: 4, referenceDepth: 100).ToList();
// segments[0]: LogRatio 0.0,   CopyNumber 2 (neutral)
// segments[1]: LogRatio -1.0,  CopyNumber 1 (deletion)
// segments[2]: LogRatio 0.585, CopyNumber 3 (duplication)
```

**Numerical walk-through:** window means 100, 50, 150 ⇒ log2(1.0)=0, log2(0.5)=−1, log2(1.5)=0.585 ⇒
`round(2·2^0)=2`, `round(2·2^−1)=1`, `round(2·2^0.585)=round(3.0)=3` [2][3].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [StructuralVariantAnalyzer_DetectCNV_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/StructuralVariantAnalyzer_DetectCNV_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [SV-CNV-001-Evidence.md](../../../docs/Evidence/SV-CNV-001-Evidence.md)
- Related algorithms: [SV_Detection](../StructuralVar/SV_Detection.md)

## 8. References

1. Yoon S, Xuan Z, Makarov V, Ye K, Sebat J. 2009. Sensitive and accurate detection of copy number variants using read depth of coverage. Genome Research 19(9):1586–1592. https://doi.org/10.1101/gr.092981.109 (open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC2752127/)
2. Talevich E, Shain AH, Botton T, Bastian BC. 2016. CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. PLoS Comput Biol 12(4):e1004873. https://doi.org/10.1371/journal.pcbi.1004873 (source `cnvlib/call.py`: https://raw.githubusercontent.com/etal/cnvkit/master/cnvlib/call.py)
3. CNVkit project. Calling copy number gains and losses (documentation). https://cnvkit.readthedocs.io/en/stable/calling.html (accessed 2026-06-13)
