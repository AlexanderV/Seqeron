---
type: concept
title: "Copy-number alteration classification (log2 ratio → absolute CN → CNA state, CNVkit thresholds)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-CNA-001-Evidence.md
  - docs/algorithms/Oncology/Copy_Number_Alteration_Classification.md
source_commit: 1555a13272537a9a4dad42dc8da28f2e27aa2180
created: 2026-07-09
updated: 2026-07-14
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-cna-001-evidence
      evidence: "Test Unit ID: ONCO-CNA-001 ... Algorithm: Copy-Number Alteration Classification (log2 copy ratio → absolute copy number → CNA state)"
      confidence: high
      status: current
---

# Copy-number alteration classification (log2 ratio → absolute CN → CNA state)

The **oncology total-copy-number classification layer**: it takes a single **log2 copy ratio** per
segment/bin and returns (a) an **absolute integer copy number** and (b) a discrete **CNA state**
— DeepDeletion / Loss / Neutral / Gain / Amplification — using **CNVkit's hard-threshold caller**
(`absolute_threshold`), corroborated by **GISTIC2** amplitude thresholds. Validated under test unit
**ONCO-CNA-001**; the literature-traced record is [[onco-cna-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the evidence-artifact pattern.

**How it differs from its neighbours (all three are distinct layers, not duplicates):**

- [[focal-amplification-detection]] (ONCO-CNA-002) is the **GISTIC2 focal/broad + oncogene-mapping**
  sibling: it asks the orthogonal **length** question (segment < 98% of its arm → focal) and maps arms
  to oncogenes. It shares only the GISTIC2 amplitude threshold (`t_amp = 0.1`) this unit uses to place
  the Amplification bin — it does **not** bin log2 into the five discrete states.
- [[homozygous-deletion-detection]] (ONCO-CNA-003) is the **loss-side consumer** of this unit's CN-0
  call: a homozygous deletion is *exactly* the `DeepDeletion` (integer CN 0) state defined here, which it
  filters for and then maps arm→tumour-suppressor. It invents no new threshold — it reuses this binning.


- [[allele-specific-copy-number-ascat]] (ONCO-ASCAT-001) is the **allele-specific** layer — it derives
  nA/nB and jointly fits purity ρ / ploidy ψ from logR **and** BAF. This unit is **total-CN only**:
  no allelic contrast, no purity fit, just a single log2 ratio → state.
- [[aneuploidy-detection]] (CHROM-ANEU-001) applies the same `n = 2·2^log2` conversion but at
  **whole-chromosome** granularity with a ≥80%-of-bins vote; this unit classifies **per segment** into
  the five discrete oncology CNA states.
- [[chromothripsis-inference]] (ONCO-SV-001) is a **downstream consumer**: it reads a run of the
  per-segment integer CN states this unit produces and classifies the *pattern* (oscillating two-state
  → Chromothripsis vs monotone rising → not) — a region-level pattern classifier, not a per-segment
  caller.
- The in-repo `StructuralVariantAnalyzer.DetectCNV` / `SegmentCopyNumber` — the germline read-depth
  [[read-depth-cnv-segmentation]] (SV-CNV-001) — already converts read-depth → log2 → **integer** CN via
  `round(2·2^log2)` and merges segments, but does **not** classify into discrete CNA states. This unit is
  the classification layer on top of that log2 ratio, swapping `round` for CNVkit's hard-threshold binning.

## 1. Absolute copy number from log2 ratio

Reference `_log2_ratio_to_absolute_pure` (CNVkit `cnvlib/call.py`), for a pure sample at ploidy 2:

```
n = ref_copies · 2^log2     # ref_copies = ploidy = 2 (diploid autosomal reference)
                            # → n = 2 · 2^log2
```

Continuous values: log2 0 → 2.0, log2 1 → 4.0, log2 −1 → 1.0, log2(3/2) → 3.0, log2 0.8 → 3.482…

## 2. Hard-threshold CNA-state calling (`absolute_threshold`)

Default thresholds `(−1.1, −0.25, 0.2, 0.7)` map to states `[0, 1, 2, 3, 4+]`. The called integer CN is
the **index of the FIRST threshold the log2 value is `<=`** (boundary **inclusive** → the *lower* bin);
if log2 exceeds every threshold, CN = `ceil(2·2^log2)`:

```
cnum = 0
for cnum, thresh in enumerate(thresholds):
    if log2 <= thresh: break            # boundary-inclusive: exactly-on-threshold → lower state
else:
    cnum = ceil(2 · 2^log2)             # above the last threshold → progressively larger integer
```

| log2 | bin | integer CN | CNA state |
|------|-----|-----------|-----------|
| −2.0 | ≤ −1.1 | 0 | **DeepDeletion** |
| −1.1 | ≤ −1.1 (boundary) | 0 | DeepDeletion |
| −1.0 | (−1.1, −0.25] | 1 | **Loss** |
| −0.25 | boundary | 1 | Loss |
| 0.0 | (−0.25, 0.2] | 2 | **Neutral** |
| 0.2 | boundary | 2 | Neutral |
| 0.585 = log2(3/2) | (0.2, 0.7] | 3 | **Gain** |
| 0.7 | boundary | 3 | Gain |
| 1.0 | > 0.7 → ceil(4.0) | 4 | **Amplification** |
| 2.0 | > 0.7 → ceil(8.0) | 8 | Amplification |

Above the last threshold the CN is **not** fixed — high amplifications get progressively larger integers
(all in the Amplification class).

## 3. Threshold derivation and corroboration

- The cutoffs derive from log-transforming integer CN over ploidy 2 (`np.log2((copy_nums + .5) / 2)`,
  the +0.5 for round-half). The docstring heuristic (50% tumor clonality, ±0.1 noise band):
  `DEL(0) < −1.1`, `LOSS(1) < −0.25`, `GAIN(3) ≥ +0.2`, `AMP(4) ≥ +0.7`.
- **GISTIC2 (Mermel 2011)** corroborates independently: a low-amplitude **noise band of ±0.1 log2**
  (subsumed by CNVkit's neutral bin (−0.25, 0.2]) and high-amplitude cutoffs **+0.848 / −0.737** for
  amplification/deletion — well outside the noise band, consistent with the AMP ≥ +0.7 / DEL ≤ −1.1
  extremes. GISTIC's `-ta` / `-td` default to 0.1.

## Corner cases and failure modes

- **NaN log2 → no-call → Neutral (CN 2):** "log2=nan found; replacing with neutral copy number."
- **Boundary inclusivity** (`log2 <= thresh`): a value exactly on a threshold takes the **lower** state.
- **Amplification is open-ended:** `ceil(2·2^log2)`, not a fixed value.
- **Purity floor:** the default thresholds are "reasonably safe" only for tumor **purity ≥ 30%**.

## Assumptions and scope

- **ASSUMPTION — diploid reference ploidy = 2** for both the classification and the absolute-CN formula
  (`ref_copies = ploidy = 2`, CNVkit autosomal default). Sex-chromosome / haploid references are out of
  scope; the caller cannot currently override the baseline, so a non-diploid genome would change every
  output.
- **Threshold source note:** the CNVkit docs page shows the germline-tuned variant (`−0.4 / 0.3`); the
  `call.py` source-code default the implementation follows is the **tumor-sample** heuristic
  (`−0.25 / 0.2`). Custom thresholds are exposable (CNVkit `-t/--thresholds`).

## Implementation surface (ONCO-CNA-001)

Four `OncologyAnalyzer` (`OncologyAnalyzer.cs`) entry points build the layer bottom-up:

- `Log2RatioToCopyNumber(log2, ploidy=2)` → continuous `n = ploidy·2^log2`.
- `CallCopyNumber(log2, thresholds, ploidy)` → the hard-threshold **integer** CN.
- `ClassifyCopyNumber(log2, thresholds, ploidy)` → a full `CopyNumberCall`
  (`Log2Ratio`, `AbsoluteCopyNumber`, `IntegerCopyNumber`, `State`).
- `ClassifyCopyNumbers(log2Ratios, …)` → per-region batch, **order- and length-preserving**
  (output element *i* ↔ input *i*).

**Validation contract:** `thresholds` must be **exactly four strictly-ascending non-NaN** values
(else `ArgumentException`); `ploidy > 0` (else `ArgumentOutOfRangeException`); a null batch enumerable
throws `ArgumentNullException`. A `null` thresholds list falls back to the default `(−1.1, −0.25, 0.2, 0.7)`.
Per-region classification is O(1) (fixed four-cutoff scan); the batch is O(m). This is total-CN
**single-region** classification — no segmentation/joining (use `StructuralVariantAnalyzer.SegmentCopyNumber`)
and no allele-specific major/minor decomposition (see ONCO-LOH-001).

A [[scientific-rigor|research-grade]] correctness reference — **not for clinical or diagnostic use.**
No source contradictions: CNVkit `call.py` + docs (thresholds, formula, binning) and GISTIC2 (Mermel
2011 + tool docs) each corroborate the ±0.1 noise band and the high-amplitude cutoffs.
