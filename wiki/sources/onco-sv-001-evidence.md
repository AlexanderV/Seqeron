---
type: source
title: "Evidence: ONCO-SV-001 (somatic complex-rearrangement classification — chromothripsis inference)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-SV-001-Evidence.md
sources:
  - docs/Evidence/ONCO-SV-001-Evidence.md
source_commit: 1d2674a92b8ed1d1afe16f362ae4457c56435ff8
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-SV-001

The validation-evidence artifact for test unit **ONCO-SV-001** — **somatic complex-rearrangement
classification (chromothripsis inference)**: given a per-segment integer copy-number profile (plus
breakpoint positions / clustered SVs) it classifies a region as **Chromothripsis** vs **NotComplex**
using the Korbel & Campbell 2013 hallmark criteria operationalised into oscillation counting, a
two-state hallmark, an SV-burden floor, and a breakpoint-clustering test. The **thirty-fourth
ingested unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in its
own concept, [[chromothripsis-inference]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (four, mutually consistent, no contradictions):**
  - **Korbel & Campbell 2013, *Cell* 152:1226–1236** (rank 1, PMID 23498933; criteria extracted via the
    open-access review PMC3861665) — the **six hallmark criteria** for inferring chromothripsis:
    (A) clustering of breakpoints; (B) regularity of **oscillating copy-number states**; (C) interspersed
    retained/lost heterozygosity; (D) haplotype-specific rearrangement prevalence; (E) randomness of
    fragment order and joins; (F) ability to walk the derivative chromosome (head↔tail alternation).
    Criterion A's **clustering test**: the null of random breakpoints predicts **exponentially distributed
    inter-breakpoint distances** — clustering = deviation toward many short distances. Criterion B: an
    **alternating profile between (canonically) two CN states**, not progressively rising amplification.
    **Clustering alone is necessary but not sufficient.**
  - **Maher & Wilson "Chromothripsis and beyond" review, PMC3861665** (rank 1) — the **first-pass
    operational threshold** (Magrangeas 2011): "10, 20, or 50 oscillating copy number changes"; the
    **lowest commonly cited first-pass cutoff is 10 oscillating CN changes**. Re-states the **two-state
    hallmark** and the **exponential null** for breakpoint spacing.
  - **Cortés-Ciriano et al. 2020, *Nat Genet* 52:331–341** (rank 1, PCAWG pan-cancer, PMC7058534) —
    quantitative segment/SV thresholds: **high-confidence = oscillation between two states across ≥ 7
    adjacent segments**; **low-confidence = 4–6 adjacent segments** (below 4 → not called); canonical
    event = **> 60% of CN segments oscillate between two states**; **minimum burden = 6 clustered
    intrachromosomal SVs** (focal events with < 6 SVs excluded). Multi-state caveat: real events may
    involve multiple chromosomes / additional CN states beyond the canonical two.
  - **Magrangeas et al. 2011, *Blood* 118:675–678** (rank 1, cited via PMC3861665) — the origin of the
    **≥ 10 oscillating-CN-change** first-pass screen (multiple myeloma).

- **Documented corner cases / failure modes:** clustering necessary-but-not-sufficient (clustered
  breakpoints with **progressively rising CN / > 2 ascending states** ⇒ NOT chromothripsis — suggests
  BFB / progressive amplification); **too few segments** (below the 10-change screen / < 4 adjacent
  segments) → NotComplex; **confidence tiers** 4–6 → low, ≥ 7 → high, < 4 → not called; **SV burden
  < 6** → not eligible.

- **Datasets (worked, deterministic):**
  - **Screening thresholds:** minimum oscillating CN changes (first-pass) = **10**; canonical states = 2;
    breakpoint-clustering null = exponential inter-breakpoint distances.
  - **Cortés-Ciriano tiers:** high-confidence ≥ 7 adjacent oscillating segments; low-confidence 4–6;
    canonical fraction oscillating > 60%; minimum clustered intrachromosomal SVs = 6.
  - **Hand-constructed two-state oscillating CN profiles (per-segment integer CN → oscillation count):**
    `2,1,2,1,2,1,2,1,2,1,2` (11 segments) → **10** transitions, states {1,2} → **Chromothripsis**
    (10 ≥ 10, 2 states); `2,1,2,1,2,1` (6 segments) → 5 transitions → **NotComplex** (5 < 10);
    `2,3,4,5,6,7` (monotone rising, no down-up reversals) → 0 oscillations, 6 states → **NotComplex**
    (progressive amplification, not two-state oscillation).

- **Coverage recommendations (8 items):** MUST-test two-state profile with exactly 10 transitions →
  Chromothripsis; 5 transitions → NotComplex; monotone-rising (no reversals, > 3 states) → NotComplex
  even with many segments; oscillation counting matches per-segment state-transition count; SV burden
  < 6 → not eligible. SHOULD-test the ≥7-vs-4–6 confidence tier and the breakpoint-clustering CV>1
  flag (regular spacing CV≈0 not flagged). COULD-test empty/null/single-segment validation.

## Deviations and assumptions

- **ASSUMPTION — oscillation = adjacent-segment CN-state reversal count.** Korbel & Campbell / Magrangeas
  count "oscillating copy number changes"; the unit operationalises an *oscillation* as a segment whose
  CN state differs from its predecessor (a state transition), additionally requiring the profile to
  alternate between a bounded number of states (≤ 3, canonically 2) per criterion B. Justified: the
  screen is explicitly a "first-pass" count of CN changes; transition counting is the minimal faithful
  realisation. The exact bracketing varies by tool; the state-transition count is the directly
  source-supported quantity.
- **ASSUMPTION — clustering summarised as coefficient of variation of inter-breakpoint gaps vs the
  exponential null.** Korbel & Campbell fix no single goodness-of-fit statistic (only that random
  breakpoints give exponentially distributed distances). The unit exposes the inter-breakpoint distances
  and flags clustering when **CV > 1** (the exponential has CV = 1; over-dispersion toward short gaps
  with a few large gaps gives CV > 1). A transparent, source-anchored summary — **not a clinical
  caller**.

A [[scientific-rigor|research-grade]] correctness reference — **not for clinical or diagnostic use.**
No source contradictions: Korbel & Campbell / PMC3861665 (six criteria, ≥10 screen, exponential null,
two-state hallmark) and Cortés-Ciriano 2020 (≥7 / 4–6 tiers, > 60% fraction, ≥6 SV floor) corroborate
one another; the two-state hallmark and the "necessary-but-not-sufficient" clustering rule are shared
across every source.
