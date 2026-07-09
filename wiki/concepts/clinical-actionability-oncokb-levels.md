---
type: concept
title: "Clinical actionability assessment (OncoKB therapeutic levels of evidence)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-ACTION-001-Evidence.md
source_commit: f43bbbd5379ef59f0f01dabd2d4a7378e8cb23da
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-action-001-evidence
      evidence: "Test Unit ID: ONCO-ACTION-001 ... Algorithm: Clinical Actionability Assessment (OncoKB Therapeutic Levels of Evidence)"
      confidence: high
      status: current
---

# Clinical actionability assessment (OncoKB therapeutic levels of evidence)

Ranking a somatic variant's **clinical actionability** — how strong the evidence is that a molecular
alteration predicts drug response — using the **OncoKB Therapeutic Levels of Evidence** (Chakravarty et
al. 2017). This is the **first ingested unit of the Oncology family** and a **pure ranking** operation:
given a variant's caller-supplied set of leveled drug associations, it returns the **maximum** level
under a fixed evidence order (plus per-axis maxima). It does **not** look up biomarkers — the
knowledgebase is a caller input, not embedded in the library. Validated under test unit
**ONCO-ACTION-001**; the record is [[onco-action-001-evidence]], [[test-unit-registry]] tracks the unit,
and [[algorithm-validation-evidence]] describes the artifact pattern.

## The seven levels (Chakravarty 2017)

Each level is an ordered evidence tier stratifying treatment implications by FDA labeling, NCCN /
professional guidelines, and scientific literature. Two implication **axes** exist — **sensitivity**
(response) and **resistance** (the `R` levels):

| Level | Category | Meaning |
|-------|----------|---------|
| **1**  | Standard Care (sensitivity) | FDA-recognized biomarker predictive of response to an FDA-approved drug, **in this indication**. |
| **2**  | Standard Care (sensitivity) | NCCN / guideline biomarker predictive of response to an FDA-approved drug, in this indication. |
| **3A** | Investigational (sensitivity) | **Compelling clinical evidence** the biomarker predicts response to a drug, in this indication. |
| **3B** | Investigational (sensitivity) | Standard/investigational biomarker predictive of response in **another indication**. |
| **4**  | Hypothetical (sensitivity) | **Compelling biological evidence** the biomarker predicts response to a drug. |
| **R1** | Standard Care (resistance) | Standard-care biomarker predictive of **resistance** to an FDA-approved drug, in this indication. |
| **R2** | Investigational (resistance) | Compelling clinical evidence the biomarker predicts **resistance** to a drug. |

**Grouping (OncoKB SOP v3):** Levels **1, 2, R1** are the *standard* (FDA/NCCN) implications; **3A, 3B,
4, R2** are *investigational*. **3A ranks above 3B** because the system was refined to deprioritize
standard-care biomarkers used outside the approved indication. Levels 1/2/R1 are **fixed by guideline
inclusion**, so conflicting literature does not move them. The system is consistent with the AMP/ASCO/CAP
Joint Consensus (Li et al. 2017), whose sibling **four-tier clinical-significance** classification is
[[cancer-variant-tier-classification-amp-asco-cap]] (Tier I–IV) — the second ingested Oncology unit.

## The three ordering axes (oncokb-annotator)

The reference implementation (`oncokb-annotator`) defines three maxima; for a variant with several
leveled drug associations the **actionable level is the maximum** under the applicable order:

```
Combined  HIGHEST_LEVEL            :  R1 > 1 > 2 > 3A > 3B > 4 > R2
Sensitive HIGHEST_SENSITIVE_LEVEL  :   1 > 2 > 3A > 3B > 4
Resistance HIGHEST_RESISTANCE_LEVEL:  R1 > R2
```

The combined order **interleaves** the two axes: resistance **R1 outranks even sensitivity Level 1**,
while **R2 sits below Level 4**. A variant can carry **both** a sensitive and a resistance level; the two
axis-specific maxima are reported independently alongside the combined maximum.

## Worked oracles

- `{2, 3A}` sensitive → highest sensitive = **Level 2**.
- `{3A, 3B, 4}` sensitive → highest sensitive = **Level 3A**.
- `{1, R1}` → combined = **R1** (R1 > 1); sensitive = 1; resistance = R1.
- `{1}` → combined / sensitive = **Level 1**.
- `{4, R2}` → combined = **Level 4** (4 > R2).
- `{R1, R2}` → highest resistance = **R1**.
- `{}` (no associations) → **NotActionable** (no level).

## Invariants and edge cases

- **INV:** the reported actionable level is the maximum of the variant's leveled associations under the
  applicable total order; the three axes never disagree on their shared members.
- **`NotActionable`:** a variant with **no leveled drug association** (VUS-like) has no highest level —
  the annotator leaves `HIGHEST_LEVEL` empty; the library models this as a distinct `NotActionable`
  outcome (the name is the library's, "no level" is OncoKB's observable).
- Null input → `ArgumentNullException`; per-variant outputs preserve input order (mirrors
  `AnnotateCancerVariants`).

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the **level-ranking** logic
only. **The knowledgebase is a caller input:** the library does **not** embed or reproduce the OncoKB
curated database (3,000+ alterations across 418 genes) — it ranks caller-supplied
drug–gene–level evidence, it does not annotate biomarkers from scratch. This is a framework boundary,
not an algorithm parameter. **Not for clinical or diagnostic use.** No source contradictions —
Chakravarty 2017, the OncoKB Levels-of-Evidence PDF, the OncoKB SOP v3, and the oncokb-annotator README
are mutually consistent.
</content>
