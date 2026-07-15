---
type: source
title: "Evidence: POP-FST-001 (Fst, F-statistics, pairwise Fst)"
tags: [validation, population-genetics]
doc_path: docs/Evidence/POP-FST-001-Evidence.md
sources:
  - docs/Evidence/POP-FST-001-Evidence.md
source_commit: 6a9852103155b627075f1a105de26fac5b97f70a
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: POP-FST-001

The validation-evidence artifact for test unit **POP-FST-001** — **Fst (fixation index),
F-statistics (Fis, Fit, Fst), and pairwise Fst**: population-differentiation measures computed from
per-population allele frequencies. It is one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the formulae, worked oracles,
invariants, and corner cases are synthesized in the dedicated concept
[[population-differentiation-fst]]. This is a population-genetics `POP-*` unit that **consumes** the
per-population allele frequencies produced by [[allele-genotype-frequencies]] (POP-FREQ-001), sits
in the family anchored by [[ancestry-estimation-admixture]], and shares its heterozygosity
denominator with [[genetic-diversity-statistics]] (POP-DIV-001). See [[test-unit-registry]] for how
units are tracked.

## What this file records

- **Online sources (all "exact match", no deviations):**
  - **Wikipedia "Fixation index"** — Fst as a 0-to-1 measure of population differentiation;
    Wright's (1965) variance definition `Fst = sigma_S^2 / (pBar(1-pBar))`; the alternative
    `Fst = (pBar(1-pBar) - mean(p(1-p))) / (pBar(1-pBar))`; Fst is **not** a metric (fails the
    triangle inequality).
  - **Wikipedia "F-statistics"** — heterozygosity partition `Fis = 1 - Hi/Hs`, `Fit = 1 - Hi/Ht`,
    `Fst = 1 - Hs/Ht`, and the exact identity `(1 - Fit) = (1 - Fis)(1 - Fst)`.
  - **Wright (1950, 1965)** — variance-based F-statistics and size-weighting `c_i = n_i/N`.
  - **Weir & Cockerham (1984)** — the θ (ANOVA variance-component, bias-corrected) estimator, cited
    to mark what the implementation deliberately does **not** do.
  - **Holsinger & Weir (2009)** review; **Cavalli-Sforza et al. (1994)** textbook — reference Fst
    values.

- **Implementation choice (documented):** `CalculateFst` uses **Wright's variance-based population
  parameter directly** from known allele frequencies (`pBar`, `sigma_S^2`, `het` hand-computed per
  locus; multi-locus = ratio-of-sums `Sum sigma_S^2 / Sum het`). It applies **no Weir–Cockerham
  finite-sample correction**. `CalculateFStatistics` uses the heterozygosity-based definitions and
  the partition identity holds algebraically.

- **Oracles & datasets:** identical populations → Fst = 0; fixed differences `p1=1, p2=0` → Fst = 1.0
  exactly (pBar=0.5, var=0.25, het=0.25); moderate case pop1=(0.9,0.8)/pop2=(0.1,0.2) → Fst = 1/2;
  strengthened exact values — unequal sizes 0.006274…, F-statistics components Fis=1/19, Fit=1/13,
  Fst=1/39; interpretation-scale exacts 1/2499 (little) and 61/198 (very great); pairwise-matrix
  cells 1/99, 4/21, 3/25; excess-heterozygosity Fis=-2/3, Fit=-2/5, Fst=4/25. Literature reference
  values from Cavalli-Sforza (1994) and Elhaik (2012, HapMap). Test count 25 (was 22 → -1 duplicate
  removed + 4 new: monomorphic→0, both-fixed-same-allele→0, pairwise exact cells,
  excess-heterozygosity negative Fis).

## Deviations and assumptions

**None** as an algorithm deviation — every formula is an exact match to Wright (1965) / Wikipedia.
The one substantive modelling choice, explicitly documented: the unit computes the **population
parameter** from known allele frequencies and **does not** implement the Weir & Cockerham (1984) θ
estimator or its unequal-sample-size bias correction. Contract behaviours: a zero denominator
(`pBar = 0` or `1`, empty populations, both-fixed-same-allele) **returns 0** for the 0/0 case; the
pairwise matrix is symmetric with a zero diagonal. Scope is Wright Fst + heterozygosity F-statistics
only — no significance testing, no bootstrap confidence intervals. No source contradictions; Open
Questions: none.
