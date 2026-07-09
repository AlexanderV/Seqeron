---
type: concept
title: "Chromothripsis inference (Korbel & Campbell hallmark criteria — oscillating CN states + breakpoint clustering)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-SV-001-Evidence.md
source_commit: 1d2674a92b8ed1d1afe16f362ae4457c56435ff8
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-sv-001-evidence
      evidence: "Test Unit ID: ONCO-SV-001 ... Algorithm: Somatic Complex Rearrangement Classification (Chromothripsis Inference)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:copy-number-alteration-classification
      source: onco-sv-001-evidence
      evidence: "The oscillation-count / two-state-hallmark test operates on the per-segment integer copy-number states that ONCO-CNA-001 produces; chromothripsis inference is a region-level pattern classifier over that CN profile, not a per-segment CN caller."
      confidence: high
      status: current
---

# Chromothripsis inference (complex-rearrangement classification)

The **oncology complex-rearrangement layer**: given a **per-segment integer copy-number profile** (plus
breakpoint positions / a count of clustered intrachromosomal SVs), it classifies an affected region as
**Chromothripsis** vs **NotComplex** using the **Korbel & Campbell 2013 hallmark criteria** — a single
catastrophic shattering-and-reassembly event, distinguished from progressive amplification by its
**oscillation between (canonically) two copy-number states** and **clustered breakpoints**. Validated
under test unit **ONCO-SV-001**; the literature-traced record is [[onco-sv-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern.

**How it differs from its neighbours (distinct region-level pattern classifier, not a duplicate):**

- [[copy-number-alteration-classification]] (ONCO-CNA-001) is the **per-segment** CN-state caller
  (log2 → integer CN → DeepDeletion…Amplification). This unit is a **region-level pattern classifier
  that consumes that profile**: it counts how the CN state *oscillates* across adjacent segments and
  whether it alternates between a bounded number of states — it invents no per-segment CN threshold.
- [[focal-amplification-detection]] (ONCO-CNA-002) asks a **length** question (focal vs arm-level) about
  a *single amplified* segment; this unit asks a **pattern** question (oscillating vs monotone) about a
  *run* of segments. A monotone rising CN profile — the very thing GISTIC-style amplification detection
  flags — is explicitly **NotComplex** here (progressive amplification / BFB, not chromothripsis).
- [[gene-fusion-detection-read-evidence]] (ONCO-FUSION-001) detects a **single** rearrangement junction
  from read evidence; chromothripsis is the **many-clustered-breakpoints** regime (≥ 6 SVs) whose
  fragment order and joins look **random**.
- [[aneuploidy-detection]] (CHROM-ANEU-001) is whole-chromosome CN gain/loss; this unit is a
  focal/intrachromosomal shattering signature within a region.
- Distinct from the gene-order [[genome-rearrangement-breakpoint-distance]] — that "breakpoint" is an
  adjacency broken between two gene orders; here a **breakpoint** is a genomic SV junction, and
  **clustering** of those junctions (vs an exponential null) is criterion A.

## 1. The six hallmark criteria (Korbel & Campbell 2013)

| # | Criterion | Operationalised signal |
|---|-----------|------------------------|
| A | Clustering of breakpoints | inter-breakpoint distances deviate from the **exponential null** toward many short gaps |
| B | Regularity of **oscillating CN states** | CN alternates between (canonically) **two** states, not progressively rising |
| C | Interspersed retained / lost heterozygosity | LOH toggles along the region (SNP-array / genotype) |
| D | Haplotype-specific prevalence | rearrangements concentrate on one haplotype |
| E | Randomness of fragment order / joins | shuffled fragment order, random join orientations |
| F | Walk the derivative chromosome | invariant head↔tail alternation |

Criteria A and B are the two the unit computes directly from a CN profile + breakpoint list; C–F are
recorded as the full hallmark set. **Clustering (A) alone is necessary but NOT sufficient** — the
oscillating-CN and randomness criteria must also hold.

## 2. Oscillation counting and the two-state hallmark (criterion B)

An **oscillation** is operationalised as a **segment whose CN state differs from its predecessor** (a
state transition); the profile must additionally alternate between a **bounded number of states (≤ 3,
canonically 2)**. Thresholds (source-traced):

- **First-pass screen (Magrangeas 2011, via PMC3861665):** the lowest commonly cited cutoff is
  **≥ 10 oscillating CN changes**.
- **Confidence tiers (Cortés-Ciriano 2020):** **≥ 7** adjacent oscillating segments → **high-confidence**;
  **4–6** → **low-confidence**; **< 4** → not called.
- **Canonical fraction (Cortés-Ciriano 2020):** **> 60%** of the region's CN segments oscillate between
  two states for a canonical (two-state) event.
- **SV-burden floor (Cortés-Ciriano 2020):** events with **< 6** clustered intrachromosomal SVs are
  excluded (not eligible for a chromothripsis call).

## 3. Breakpoint clustering test (criterion A)

The null hypothesis of **random** breakpoints predicts **exponentially distributed inter-breakpoint
distances** (CV = 1). The unit exposes the distances and their mean and **flags clustering when the
coefficient of variation exceeds 1** — over-dispersion toward many short gaps with a few large gaps
gives CV > 1; regular spacing gives CV ≈ 0 (not flagged). A transparent summary, not a fixed
goodness-of-fit test.

## Worked dataset (per-segment integer CN → oscillation count)

| Profile (integer CN per segment) | Oscillations | States | Verdict |
|----------------------------------|-------------|--------|---------|
| `2,1,2,1,2,1,2,1,2,1,2` (11 seg) | 10 | {1,2} (2) | **Chromothripsis** (10 ≥ 10, two-state) |
| `2,1,2,1,2,1` (6 seg) | 5 | {1,2} | **NotComplex** (5 < 10 screen) |
| `2,3,4,5,6,7` (monotone rising) | 0 (no down-up reversals) | 6 | **NotComplex** (progressive amplification / BFB, not two-state) |

## Corner cases and failure modes

- **Clustering necessary but not sufficient:** clustered breakpoints with **progressively rising CN**
  (> 2 ascending states, gradual amplification) is **NOT** chromothripsis — it suggests
  breakage-fusion-bridge / progressive amplification.
- **Too few segments / transitions:** below the ≥ 10 first-pass screen (or < 4 adjacent oscillating
  segments) → NotComplex (insufficient evidence).
- **Confidence tiering:** 4–6 oscillating segments → low-confidence; ≥ 7 → high-confidence; < 4 → not
  called.
- **SV-burden floor:** a focal event with **< 6** SVs is excluded regardless of oscillation.
- **Multi-state caveat:** canonical events oscillate between two CN states, but real events may span
  multiple chromosomes / additional CN states (Cortés-Ciriano 2020).
- **Empty / null / single-segment** input → explicit validation.

## Assumptions and scope

- **ASSUMPTION — oscillation = adjacent-segment CN-state reversal count.** Korbel & Campbell / Magrangeas
  count "oscillating copy number changes" without fixing one bracketing; transition counting (with the
  ≤ 3-state bound) is the minimal faithful realisation of the explicitly "first-pass" screen, and the
  state-transition count is the directly source-supported quantity.
- **ASSUMPTION — clustering summarised by CV > 1 vs the exponential null.** The sources fix no single
  goodness-of-fit statistic; the exponential's CV = 1 anchors the transparent CV > 1 flag.
- **Scope:** criteria C–F (heterozygosity toggling, haplotype prevalence, fragment-order randomness,
  derivative-chromosome walk) are recorded as the full hallmark set but are not the computed signal;
  this is a **research-grade** pattern classifier, not a clinical caller.

A [[scientific-rigor|research-grade]] correctness reference — **not for clinical or diagnostic use.**
No source contradictions: Korbel & Campbell 2013 / PMC3861665 (six criteria, ≥ 10 first-pass screen,
exponential-null clustering, two-state hallmark) and Cortés-Ciriano 2020 (≥ 7 / 4–6 confidence tiers,
> 60% oscillating fraction, ≥ 6 SV floor) corroborate one another throughout.
