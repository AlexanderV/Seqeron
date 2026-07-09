---
type: concept
title: "Driver-gene classification (20/20 rule)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-DRIVER-001-Evidence.md
source_commit: f640eb404dd41ebb270208e6664505bba6c4cb8e
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-driver-001-evidence
      evidence: "Test Unit ID: ONCO-DRIVER-001 ... Algorithm: Driver Mutation Detection (20/20 rule)"
      confidence: high
      status: current
---

# Driver-gene classification (20/20 rule)

Classifying a cancer gene as an **oncogene (OG)**, a **tumor-suppressor gene (TSG)**, or **ambiguous**
from the *pattern* of its recorded somatic mutations — the **20/20 rule** of Vogelstein et al. (2013,
*Science* "Cancer Genome Landscapes"). This is the **twelfth ingested unit of the Oncology family**
(`IdentifyDriverMutations` + `MatchCancerHotspots` + `ScoreDriverPotential`) and a **pure per-gene
mutation-pattern heuristic** — it counts what *kind* of mutations a gene accumulates, not any external
pathogenicity score. Validated under test unit **ONCO-DRIVER-001**; the record is
[[onco-driver-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern. It operates at the **gene** level,
distinct from — and complementary to — the **variant-level** clinical classifiers
[[cancer-variant-tier-classification-amp-asco-cap]] (AMP/ASCO/CAP tiers) and
[[clinical-actionability-oncokb-levels]] (OncoKB therapeutic levels).

## The 20/20 rule (Vogelstein 2013; verbatim via Tokheim 2020)

Two thresholds over the gene's own mutation catalog:

| Call | Criterion | Mode of action |
|------|-----------|----------------|
| **Oncogene** | **> 20%** of mutations are **missense at recurrent positions** | activating (gain of function) |
| **TumorSuppressor** | **> 20%** of mutations are **inactivating / truncating** | loss of function |

- **Recurrent position** = a single protein position observed **≥ 2 times** (Miller et al. 2017). A
  singleton missense is *not* recurrent.
- **Truncating (inactivating) mutation types** — protein-truncating changes: **nonsense**,
  **frameshift** insertions/deletions, **splice donor/acceptor** site mutations, and **gained or lost
  stop codons** (Schroeder et al. 2014 / OncodriveROLE; Tokheim 2020 names nonsense + frameshift).
- The rule is **lenient**: "all well-documented cancer genes far surpass these criteria" — a dual pass
  of both thresholds is atypical.

## Worked archetypes (oracles)

- **IDH1 — oncogene archetype.** 10 missense mutations, **all at codon 132** (R132H): truncating
  fraction **0/10 = 0.00**, recurrent-missense fraction **10/10 = 1.00 ≫ 0.20 → Oncogene**. When
  virtually every missense falls on one codon the recurrent-missense fraction is ~1.0.
- **Dispersed-truncating — TSG archetype.** 5 nonsense + 2 frameshift + 1 splice = **8 truncating** at
  distinct positions, plus 2 missense at distinct positions: truncating fraction **8/10 = 0.80 > 0.20 →
  TumorSuppressor**; recurrent-missense fraction 0.00 (no missense position seen ≥ 2×).
- **Boundary (strict `>`).** Truncating fraction **exactly 0.20 → NOT a TSG**; **0.30 → TumorSuppressor**.

## Decision rule and its deviations

```
og_frac  = (# recurrent-position missense) / (total mutations)
tsg_frac = (# truncating / inactivating)   / (total mutations)
if og_frac > 0.20 and tsg_frac > 0.20  → classify by the LARGER fraction; exact tie → Ambiguous
elif og_frac  > 0.20                   → Oncogene
elif tsg_frac > 0.20                   → TumorSuppressor
else                                   → Ambiguous          # low recurrence / neither criterion met
```

- **Strict `> 0.20` (assumption).** Vogelstein / Tokheim / Miller state ">20%"; OncodriveROLE writes
  "≥20%" for TSGs. The unit uses **strict `>` 0.20 for both**, matching the primary source. A fraction of
  exactly 0.20 is therefore *not* sufficient. (Flagged as the one place sources differ in glyph.)
- **Dual-pass tie-break (assumption).** Sources prescribe no single label when a gene passes both
  thresholds; the unit picks the **dominant fraction** and reports **Ambiguous on an exact tie** — the
  least-surprising deterministic resolution given that well-documented genes "far surpass" one criterion.
- **Ambiguous is correct conservative behaviour.** A gene satisfying **neither** criterion (few, lowly
  recurrent mutations) is left **Ambiguous** — OncodriveROLE documents that the rule "fails to identify
  drivers … mostly the lowly recurrent ones." The 20/20 rule is a **heuristic, not a statistical test**;
  it can misclassify genes with few mutations (random passenger truncations in an OG drift up in
  frequency and can mislead — Tokheim 2020). Statistical successors (20/20+, MutSigCV) address this; they
  are out of this unit's scope.

## Companion operations

- **`MatchCancerHotspots`** — flags a mutation whose **(gene, position)** is present in a
  **caller-supplied hotspot set**, and not otherwise. Recurrence/hotspot membership is a caller input
  (e.g. a cancerhotspots.org-style catalog), not embedded.
- **`ScoreDriverPotential`** — returns the **20/20-rule driver-signal fraction in [0,1]** = **max of the
  two criterion fractions** as a transparent, source-derived score. The checklist names CADD / SIFT /
  PolyPhen, but those are **externally trained models** that cannot be reproduced here; external
  pathogenicity scores are **caller-supplied / not implemented** (documented, not fabricated).
- **`IdentifyDriverMutations`** — returns a **subset** of the input somatic variants (invariant:
  driver ⊆ somatic).

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the **20/20 mutation-pattern
rule** only. **Inputs are caller-supplied:** the library does not curate cancer-gene catalogs, hotspot
databases, or trained pathogenicity models — it counts mutation classes and applies the thresholds.
Gene-level classification here is **orthogonal** to the variant-level clinical-significance tiers
([[cancer-variant-tier-classification-amp-asco-cap]]) and therapeutic actionability levels
([[clinical-actionability-oncokb-levels]]): a driver *gene* call is not a per-variant clinical call.
**Not for clinical or diagnostic use.** No source contradictions — Vogelstein 2013, Tokheim & Karchin
2020, Schroeder 2014 (OncodriveROLE), and Miller 2017 are mutually consistent (the sole `≥` vs `>` glyph
difference is resolved in favour of the primary source's strict ">20%").
</content>
</invoke>
