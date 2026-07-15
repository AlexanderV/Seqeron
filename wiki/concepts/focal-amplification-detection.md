---
type: concept
title: "Focal amplification detection (GISTIC2 length-based focal/broad split + oncogene mapping)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-CNA-002-Evidence.md
  - docs/algorithms/Oncology/Focal_Amplification_Detection.md
source_commit: e8b2df0e8025da9158f1fd12db29d170e96ceeb3
created: 2026-07-09
updated: 2026-07-14
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-cna-002-evidence
      evidence: "Test Unit ID: ONCO-CNA-002 ... Algorithm: Focal Amplification Detection (GISTIC2 length-based focal/broad split + oncogene mapping)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:copy-number-alteration-classification
      source: onco-cna-002-evidence
      evidence: "Both are GISTIC2-informed CNA units; CNA-002's amplitude test reuses the GISTIC2 t_amp=0.1 threshold that CNA-001's Amplification state also cites, but CNA-002 classifies focal-vs-broad by LENGTH and maps to oncogenes rather than binning log2 into five states."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:homozygous-deletion-detection
      source: onco-cna-002-evidence
      evidence: "ONCO-CNA-003 is the deletion counterpart: IdentifyDeletedTumorSuppressors mirrors this unit's IdentifyAmplifiedOncogenes (arm→gene panel), filtering to CN-0 homozygous deletions instead of focal amplifications."
      confidence: high
      status: current
---

# Focal amplification detection (length-based focal/broad split + oncogene mapping)

The **oncology focal-amplification layer**: given per-segment copy-number data it (a) keeps only
segments that are **amplified** (gain amplitude above a threshold) **and** **focal** (short relative
to their chromosome arm), then (b) **maps** each focal amplification's arm label to a panel of known
**oncogenes**. Two operations — `DetectFocalAmplifications` and `IdentifyAmplifiedOncogenes`.
Validated under test unit **ONCO-CNA-002**; the literature-traced record is [[onco-cna-002-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern.

**How it differs from its neighbours (distinct layer, not a duplicate):**

- [[copy-number-alteration-classification]] (ONCO-CNA-001) bins a single **log2 ratio** into five
  discrete **CNA states** (DeepDeletion…Amplification) with no notion of segment length or genomic
  extent. This unit is orthogonal: it asks **"focal or broad?"** — a **length** question GISTIC2 answers
  by comparing the segment to its chromosome **arm**, not by amplitude alone. The two share only the
  GISTIC2 amplitude threshold (`t_amp = 0.1`) that CNA-001 uses to place the Amplification bin.
- [[allele-specific-copy-number-ascat]] (ONCO-ASCAT-001) and [[aneuploidy-detection]] (CHROM-ANEU-001)
  produce copy-number values; this unit is a **downstream filter + annotator** over segments, not a
  caller.

## 1. The two-part focal-amplification predicate

A segment is a **focal amplification** iff BOTH hold:

```
amplified:  log2 > t_amp                       # t_amp = 0.1 (GISTIC2 amplification threshold)
focal:      SegLen / ArmLength < broad_len_cutoff   # broad_len_cutoff = 0.98 (fraction of arm)
```

- **Amplitude gate (`t_amp = 0.1`, GISTIC2 docs):** "Regions with a copy number gain above this
  positive value are considered amplified." Because any single-copy gain has log2 ≥ log2(3/2) ≈ 0.585
  (CNVkit) — well above 0.1 — the 0.1 gate admits all true gains and rejects only near-neutral noise.
- **Length gate (`broad_len_cutoff = 0.98`, Mermel 2011 / GISTIC2 docs):** GISTIC2.0's key move is to
  separate **focal** from **arm-level/broad** SCNAs **purely by length**: an event occupying
  **> 98% of a chromosome arm** is arm-level (broad); **< 98%** is focal. Length "provides a natural
  basis for classifying events" — GISTIC2.0 shifted away from amplitude-based filtering to length.

**Boundary rule:** the comparison is **strict** — fraction `< 0.98` is focal; a segment at **exactly
0.98** (or above) is arm-level and NOT reported. (Paper wording: "more than 98% ⇒ arm-level.")

## 2. Oncogene mapping (`IdentifyAmplifiedOncogenes`)

Each focal amplification carries an **arm label** (chromosome + p/q arm, e.g. `17q`). The mapper matches
that arm prefix against a small registry of recurrently amplified **oncogenes** and their cytogenetic
locations (NCBI Gene):

| Oncogene | Arm | Cytoband |
|----------|-----|----------|
| ERBB2 | 17q | 17q12 |
| MYC | 8q | 8q24.21 |
| EGFR | 7p | 7p11.2 |
| CCND1 | 11q | 11q13.3 |
| MDM2 | 12q | 12q15 |
| CDK4 | 12q | 12q14.1 |

A single arm can carry **multiple** oncogenes (12q → both MDM2 and CDK4). Only **focal amplifications**
feed the mapper — a neutral/loss or a broad arm-level event never yields an oncogene amplification. The
arm→oncogene panel is the algorithm's built-in registry (the arm+p/q prefix is what a segment is matched
against), not a caller-supplied knowledgebase.

## Implementation surface (ONCO-CNA-002)

The spec (`docs/algorithms/Oncology/Focal_Amplification_Detection.md`) fixes the API on
`OncologyAnalyzer` (`Seqeron.Genomics.Oncology`):

- `DetectFocalAmplifications(segments, thresholds?)` → `IReadOnlyList<CopyNumberArmSegment>`: single-pass
  **filter** returning the input segments that are focal amplifications, **in input order** (constructs no
  new segments; INV-03). `thresholds` is a `FocalAmplificationThresholds?`; `null` ⇒ GISTIC2 defaults
  (`t_amp = 0.1`, `broad_len_cutoff = 0.98`).
- `IdentifyAmplifiedOncogenes(amplifications)` → `IReadOnlyList<string>`: distinct panel oncogene symbols
  on amplified arms, **in panel order**.
- `IsFocalAmplification(segment, thresholds)`: single-segment predicate (public helper for reuse).

`CopyNumberArmSegment` carries the arm label, `Start`/`End`, `ArmLength`, and mean `Log2Ratio`; segment
length is `End − Start` (base-pair counts). **Validation:** null `segments`/`amplifications` ⇒
`ArgumentNullException`; a segment with `ArmLength ≤ 0` or `End ≤ Start` ⇒ `ArgumentException`. Arm labels
match **case-insensitively (Ordinal-ignore-case)**; arms outside the six-gene panel map to no oncogene.

**Complexity:** `DetectFocalAmplifications` is `O(n)` time / `O(k)` space (n segments, k focal);
`IdentifyAmplifiedOncogenes` is `O(n + g)` with fixed panel size `g = 6`. No substring/pattern search is
involved, so the repository suffix tree is **not applicable**. **Segmentation is upstream**, not performed
here — `StructuralVariantAnalyzer.SegmentCopyNumber` (SV-CNV-001 / ONCO-CNA-001) produces the segments this
filter consumes.

**Invariants** (INV-01…04): every reported amplification satisfies `L/A < 0.98` (strict) and `r > t_amp`;
the output is an in-order subset of the input; an oncogene is reported only for an arm carrying a focal
amplification.

## Worked dataset (arm length 1,000,000 bp; t_amp 0.1; cutoff 0.98)

| Segment | Arm | SegLen | SegLen/Arm | log2 | Amplified? | Focal? | Focal amp? |
|---------|-----|--------|-----------|------|-----------|--------|-----------|
| A (ERBB2) | 17q | 500,000 | 0.50 | 1.0 | yes | yes | **yes** |
| B (whole arm) | 8q | 990,000 | 0.99 | 1.5 | yes | no | no (arm-level) |
| C (low amp) | 7p | 300,000 | 0.30 | 0.05 | no | yes | no (not amplified) |
| D (boundary) | 11q | 980,000 | 0.98 | 1.0 | yes | no (= cutoff) | no (not < 0.98) |

Segment A on 17q maps to **ERBB2**; B/C/D produce no oncogene amplification.

## Corner cases and failure modes

- **Whole-arm event excluded as focal:** a segment at ≥ 98% of its arm is arm-level/broad — never
  reported as focal even when highly amplified (segment B).
- **Boundary at exactly 0.98:** treated as arm-level (strict `< 0.98` → focal); segment D fails.
- **Amplitude below `t_amp`:** a segment whose log2 gain does not exceed 0.1 is not amplified and is
  excluded regardless of length (segment C).
- **Non-amplified segments never map:** only focal amplifications feed `IdentifyAmplifiedOncogenes`.
- **Null / empty input:** documented handling (no focal amplifications).

## Assumptions and scope

- **ASSUMPTION — amplitude gate combined with the length rule.** GISTIC2.0 itself classifies
  focal-vs-broad **purely by length**; the "amplified" gate (`t_amp = 0.1`, gain above +0.1) is taken
  from the GISTIC2 `t_amp` parameter. Combining the length rule (paper) with the amplitude rule (docs)
  into a single `DetectFocalAmplifications` predicate is this unit's documented **integration choice**,
  not an invented one — both halves are source-backed.
- **ASSUMPTION — arm fraction supplied as input.** GISTIC2 derives arm boundaries from a cytoband file;
  this unit bundles no cytoband table. The **caller supplies each segment's arm label and the arm's
  length**, and the algorithm computes `SegLen / ArmLength`. The 0.98 cutoff and the amplitude rule are
  unchanged.
- **Scope note:** deletions are out of scope here (GISTIC2 `t_del`); the deletion mirror is
  [[homozygous-deletion-detection]] (ONCO-CNA-003) — CN-0 homozygous-deletion filter + tumour-suppressor
  mapping, the loss-side counterpart of this unit's `IdentifyAmplifiedOncogenes`.

A [[scientific-rigor|research-grade]] correctness reference — **not for clinical or diagnostic use.**
No source contradictions: Mermel 2011 (98%-of-arm length rule), GISTIC2 docs (`broad_len_cutoff` 0.98,
`t_amp` 0.1), CNVkit (single-copy gain log2 0.585 > 0.1; focal = extent-limited), and NCBI Gene (oncogene
cytobands) corroborate one another across the whole predicate.
