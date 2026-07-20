---
type: concept
title: "Read-depth CNV detection (windowed depth → log2 ratio → integer copy number → del/dup)"
tags: [structural-variant, algorithm]
mcp_tools:
  - segment_copy_number
sources:
  - docs/algorithms/StructuralVar/Copy_Number_Variation.md
  - docs/Evidence/SV-CNV-001-Evidence.md
  - docs/Validation/reports/SV-CNV-001.md
source_commit: 941abc5bdd3a34832b55fcd0f9310c3f6150554a
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: sv-cnv-001-evidence
      evidence: "Test Unit ID: SV-CNV-001, Algorithm: Read-Depth Copy Number Variation Detection (windowed read depth → log2 ratio → integer copy number)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:breakpoint-detection-split-reads
      source: sv-cnv-001-evidence
      evidence: "SV-CNV-001 is the read-depth-segmentation member of the germline structural-variant family anchored by SV-BREAKPOINT-001; distinct method (read depth, not split-read junctions)."
      confidence: high
      status: current
---

# Read-depth CNV detection (windowed depth → log2 ratio → integer copy number)

The **read-depth member of the germline structural-variant (SV) family** (SV-CNV-001). Where the SV
anchor [[breakpoint-detection-split-reads]] localizes a breakpoint from **split-read junctions**, this
unit calls **copy-number variants (deletions / duplications) from read depth of coverage** — a
genuinely distinct method (aggregate depth signal, not per-read junction geometry). Validated under
test unit **SV-CNV-001** ([[sv-cnv-001-evidence]]); [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern. The two-stage validation
verdict is recorded in [[sv-cnv-001-report]] — **Stage A PASS · Stage B PASS-WITH-NOTES,
State ✅ CLEAN** (one round-half-to-even rounding defect fixed in-session; full suite 6493/0).

The evidence is the standard **read-depth (RD) CNV** paradigm — **Yoon et al. 2009** (the RD∝CN
hypothesis, windowed counting, GC correction) plus the **CNVkit** reference implementation
(`cnvlib/call.py`) for the log2-ratio → absolute-copy-number arithmetic.

## The read-depth → copy-number pipeline

1. **Read depth is proportional to copy number** (Yoon 2009 core hypothesis): the GC-adjusted mean RD
   of a region is a linear measurement of its copy number (observed across CN 1, 2, 3).
2. **Windowed counting** — depth is summarised over **non-overlapping windows of fixed size**, each
   read assigned once by its start position; the per-window mean RD is the signal.
3. **(Optional) GC correction** — per-window counts are normalised to the overall median:
   `r_i' = r_i · m / m_GC`, where `m_GC` is the median of windows at the same G+C%, and `m` is the
   overall median of all windows (Yoon 2009, verbatim).
4. **log2 ratio** — `log2(observed depth / reference depth)`, equivalently `log2(observed CN / ploidy)`
   (CNVkit). The **reference** (log2 = 0 anchor) defaults to the **overall median of the window means**
   (self-reference) when no external baseline is supplied.
5. **Integer copy number** — CNVkit "round" method: `CN = round(ploidy · 2^log2)` = `round(2 · 2^log2)`
   for a diploid (ploidy 2) reference.
6. **Call + segment** — classify each window (Deletion / Neutral / Duplication) and merge adjacent
   like-called windows into copy-number **segments** (`DetectCNV` / `SegmentCopyNumber`).

## The diploid log2 → CN anchors (verbatim from CNVkit)

| observed / reference RD | log2 ratio | CN = round(2·2^log2) | call |
|---|---|---|---|
| 0.5× | log2(1/2) = **−1.0** | 1 | Deletion (loss) |
| 1.0× | log2(2/2) = **0.0** | 2 | Neutral |
| 1.5× | log2(3/2) = **+0.585** | 3 | Duplication (gain) |
| 2.0× | log2(2.0) = **+1.0** | 4 | Duplication (amplification) |

CN = `ploidy · 2^log2` is **strictly increasing in log2** and **monotonic non-decreasing in RD**;
copy number is clamped **≥ 0** (CNVkit `max(0.0, ncopies)` — a negative CN is non-physical).

## Corner cases and failure modes

- **Zero-depth window** — RD = 0 makes `log2(0/ref) = −∞` undefined; the window is a **no-call**
  (homozygous-deletion candidate), **not** a `−∞` log2 (Yoon: RD is a raw read count).
- **NaN / no-signal window → Neutral (CN 2)** — CNVkit's `absolute_threshold` replaces a `nan` log2 with
  the neutral reference copy number ("log2=nan found; replacing with neutral copy number").
- **Negative extrapolation clamped** to CN 0.

## Assumptions and scope

- **ASSUMPTION — reference (log2 = 0) baseline = overall median of the non-zero window means** when no
  external reference is supplied, mirroring Yoon's overall-median `m` and CNVkit's "ratio against
  reference". Correctness-affecting (sets the neutral anchor) but source-supported; an explicit baseline
  may override it.
- **ASSUMPTION — diploid ploidy (2)** as the copy-number baseline (`CN = 2·2^log2`), the standard human
  autosomal baseline stated by both cited sources.

A [[research-grade-limitations|research-grade]] method, **not for clinical use.** It consumes
**windowed depth** (per-position depth summarised into windows — cf. the exact per-base
[[coverage-depth-calculation]]), not raw BAM.

## Method contract (algorithm spec)

The canonical PRIMARY spec is `docs/algorithms/StructuralVar/Copy_Number_Variation.md`
(status *Simplified*). Two entry points in `Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs`:

- **`DetectCNV(IReadOnlyList<int> depthData, int windowSize = 100, double? referenceDepth = null,
  string chromosome = "chr1")`** — the canonical path: partition per-position depth into
  `floor(n / windowSize)` **non-overlapping full windows**, take each window's **mean** depth, convert
  to `log2 = log2(RD_w / RD_ref)` and `CN = round(2·2^log2)`, and emit **one `CopyNumberSegment` per
  full non-zero window**. `windowSize` default **100** (literature default [Yoon 2009]).
- **`SegmentCopyNumber(IEnumerable<double> logRatios, string chromosome = "chr1")`** — the segmentation
  variant: convert each per-window log2 ratio to CN with the same rule and **merge maximal runs of
  consecutive equal-CN windows** into one segment (mean log2 over the run, window count as probe count).
- Private helpers: **`LogRatioToCopyNumber`** = `max(0, round(2·2^log2))` (the CN clamp+round rule);
  **`OverallMedianNonZero`** = the overall-median-of-non-zero-window-means baseline.

**Input contract** — `depthData`: `IReadOnlyList<int>`, values ≥ 0, 0-based single contig; `referenceDepth`:
`double?`, `null` ⇒ overall median of non-zero window means (an explicit value overrides, and is used
as-is); `chromosome`: label recorded on each segment. For `SegmentCopyNumber`, `logRatios` window *i*
covers coordinate *i*; **NaN entries are no-calls** (dropped) — note this is the *implemented* rule,
distinct from the CNVkit-documented `nan`→neutral-CN-2 replacement cited above as background theory.

**Output — `CopyNumberSegment`:** `Start`/`End` (inclusive 0-based window bounds for `DetectCNV`, or
window indices for `SegmentCopyNumber`); `LogRatio` (per-window log2, or mean log2 over a merged run);
`CopyNumber` = `round(2·2^LogRatio)` ≥ 0; `BAlleleFrequency` = **`NaN`** (a depth-only caller has no
allele-specific data); `ProbeCount` = positions in the window (`DetectCNV`) or windows merged
(`SegmentCopyNumber`).

**Invariants (verified against the spec):** INV-01 `RD_w = RD_ref` ⇒ log2 = 0 ⇒ CN = 2; INV-02
`CN = round(2·2^log2)`; INV-03 `CN ≥ 0`; INV-04 CN monotonically non-decreasing in `RD_w`; INV-05
depth partitioned into non-overlapping windows, **trailing partial window dropped**; INV-06 a
**zero-depth window is a no-call** (excluded), never a `−∞` log2.

**Edge cases:** empty `depthData` ⇒ empty result; `null` `depthData`/`logRatios` ⇒ `ArgumentNullException`;
`windowSize ≤ 0` ⇒ `ArgumentOutOfRangeException`; `windowSize` > length ⇒ empty result (no full window,
INV-05); computed reference depth ≤ 0 (all windows zero) ⇒ no segments emitted.

**Complexity:** `DetectCNV` O(n) time / O(n / windowSize) space (one summation pass plus an O(w log w)
median over `w` window means, `w ≪ n`); `SegmentCopyNumber` O(m) over `m` log2 ratios. No substring
search — the repository suffix tree is not applicable (windowed numeric aggregation over a depth array).

Worked oracle (`referenceDepth = 100`, `windowSize = 4`, depths `100×4, 50×4, 150×4`): window means
100/50/150 ⇒ log2 0.0/−1.0/0.585 ⇒ CN **2 / 1 / 3** (neutral / deletion / duplication).

## Relation to other copy-number units

- Distinct from the SV anchor [[breakpoint-detection-split-reads]] (split-read *junctions* → single-base
  breakpoint) — same germline-SV family, orthogonal evidence (aggregate depth vs per-read clip).
- The oncology **classification layer** [[copy-number-alteration-classification]] (ONCO-CNA-001) is built
  *on top of* this same log2 ratio: it swaps this unit's `round(2·2^log2)` for CNVkit's **hard-threshold
  binning** into five discrete CNA states (DeepDeletion / Loss / Neutral / Gain / Amplification). This
  unit stops at the integer CN and del/dup call.
- The allele-specific tumor layer [[allele-specific-copy-number-ascat]] adds BAF and a purity/ploidy fit;
  this germline unit is **total-CN only** from depth, no allelic contrast.

**Sharp edge:** [[read-depth-cnv-misses-copy-neutral-and-nan-to-cn2]] — read-depth-only: **blind to copy-neutral** events; NaN windows silently become **neutral CN2**.
