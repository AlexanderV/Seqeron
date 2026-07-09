---
type: concept
title: "Cancer variant tier classification (AMP/ASCO/CAP 2017 four-tier)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-ANNOT-001-Evidence.md
source_commit: c38dde50786144e5976f703468375f68abe6fd1b
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-annot-001-evidence
      evidence: "Test Unit ID: ONCO-ANNOT-001 ... Algorithm: Cancer-Specific Variant Annotation (AMP/ASCO/CAP 2017 four-tier clinical-significance classification)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:clinical-actionability-oncokb-levels
      source: onco-annot-001-evidence
      evidence: "Li et al. (2017) Tables 4-6 list Somatic database: COSMIC ... as an evidence source; the OncoKB SOP is explicitly consistent with the AMP/ASCO/CAP Joint Consensus (Li et al. 2017). Both classify somatic variants from caller-supplied curated evidence."
      confidence: high
      status: current
---

# Cancer variant tier classification (AMP/ASCO/CAP 2017 four-tier)

Classifying a somatic sequence variant's **clinical significance** into one of **four tiers** under the
**AMP/ASCO/CAP 2017 Joint Consensus** (Li MM et al. 2017, *J Mol Diagn* 19(1):4–23). This is the
**second ingested unit of the Oncology family** (`AnnotateCancerVariants` + `GetCOSMICAnnotation`) and a
**pure decision rule**: given a variant's caller-supplied **evidence level** (A–D or none),
**population MAF**, and **cancer-association flag**, it applies Figure 2 of the guideline verbatim to
return a tier. It does **not** look up biomarkers or reproduce curated databases — the evidence inputs
are caller-supplied. Validated under test unit **ONCO-ANNOT-001**; the record is
[[onco-annot-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern. It is distinct from — but consistent
with — the sibling therapeutic-actionability ranking [[clinical-actionability-oncokb-levels]] (OncoKB
levels): the OncoKB SOP is explicitly consistent with this AMP/ASCO/CAP consensus.

## The four tiers (Li et al. 2017, abstract)

| Tier | Clinical significance | Evidence-level mapping (Figure 2) |
|------|-----------------------|-----------------------------------|
| **I**   | Strong clinical significance    | Level **A or B** evidence |
| **II**  | Potential clinical significance | Level **C or D** evidence |
| **III** | Unknown clinical significance (VUS) | No evidence level; rare + a cancer association |
| **IV**  | Benign / likely benign          | No evidence level; common (MAF ≥ 1%) **or** no cancer association |

## The four evidence levels (Table 3)

Levels apply across three **categories** — Therapeutic, Diagnosis, Prognosis:

- **Level A** — biomarker predicting response/resistance to an **FDA-approved therapy** for the specific
  tumor type, or included in professional guidelines.
- **Level B** — based on **well-powered studies with expert consensus**.
- **Level C** — FDA/guideline therapies for a **different** tumor type (off-label), clinical-trial
  inclusion criteria, or diagnostic/prognostic significance from multiple small studies.
- **Level D** — **preclinical** (plausible therapeutic significance), or diagnostic/prognostic from small
  studies / few case reports without consensus.

## The decision rule (Figure 2, applied verbatim)

```
if evidence_level ∈ {A, B}        → Tier I   (strong)
elif evidence_level ∈ {C, D}      → Tier II  (potential)
elif MAF ≥ 0.01                   → Tier IV  (benign — common polymorphism)
elif not cancer_association       → Tier IV  (benign — no published cancer association)
else                              → Tier III (unknown significance)
```

**Population frequency is the Tier III ↔ Tier IV discriminator** when no clinical evidence level is
present: a **rare** variant **with** a cancer association is Tier III; a **common** variant (MAF ≥ 1%) or
one with **no** cancer association is Tier IV. The **1% (0.01) MAF cutoff** is the guideline's *primary*
benign cutoff — the work group's recommendation in the absence of paired normal tissue; the guideline
notes there is **no universally standardized MAF cutoff**. The cutoff is inclusive: **MAF exactly 0.01 →
Tier IV**.

**Evidence level dominates frequency:** a clinically significant biomarker stays Tier I even at high MAF
(well-studied germline-counterpart variants such as `TP53`, `PTEN` may appear in population databases yet
remain significant). Categorization is by evidence level (Figure 2), not downgraded by frequency.

## Worked oracles (Figure 2 / Tables 3–7)

- Level **A**, any MAF, any assoc → **Tier I** (`BRAF p.V600E`, the guideline's canonical Level-A/B
  biomarker).
- Level **B** → **Tier I**; Level **C** → **Tier II**; Level **D** → **Tier II**.
- No level, MAF **0.25**, assoc **false** → **Tier IV** (common SNP).
- No level, MAF **0.0001**, assoc **true** → **Tier III** (rare VUS with cancer association).
- No level, MAF **< 0.01**, assoc **false** → **Tier IV** (no published cancer association).
- Boundary: MAF **exactly 0.01** → **Tier IV**; MAF **0.0099** + assoc → **Tier III**.

## COSMIC annotation (`GetCOSMICAnnotation`)

**COSMIC** (Tate JG et al. 2019, *Nucleic Acids Res* 47(D1):D941) is an **external curated database**
(v86: ~6M coding mutations across 1.4M tumour samples from 26,000+ publications) — it **cannot be
hardcoded** in the library. `GetCOSMICAnnotation` is a lookup against a **caller-supplied catalog**:
returns the catalog value on a **hit**, and **null on a miss** (`not found`, never a fabricated
annotation). COSMIC presence is a *somatic-database evidence input* supporting a tier (Li Tables 4–6:
Tier I "Most likely present", Tier II "Likely present"), not by itself determining it.

## Invariants and edge cases

- **INV:** `AnnotateCancerVariants` returns **one annotation per input variant, in input order**; mixed
  batches carry per-variant tiers. Empty collection → empty list.
- **`GetCOSMICAnnotation`:** hit → catalog value; miss → **null** (external-catalog boundary).
- Invalid MAF (negative, > 1, NaN) → `ArgumentOutOfRangeException`; null inputs → `ArgumentNullException`
  (API contract per sibling methods).

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the **four-tier decision rule**
only. **Evidence inputs are caller-supplied:** the library does **not** reproduce professional
guidelines, population databases, somatic databases, or literature — the caller performs those lookups
and supplies the evidence level, population MAF, and cancer-association flag; the Figure 2 rule is applied
verbatim to whatever is supplied. The Tier III ↔ IV discriminator (MAF ≥ 1% **or** absent association →
Tier IV) is a direct reading of the Figure 2 boxes + Tables 6/7, not invented. **Not for clinical or
diagnostic use.** No source contradictions — Li et al. 2017 and Tate et al. 2019 are mutually consistent
(the "≥" glyph lost in PDF extraction is fixed by the guideline's 1% primary benign cutoff).
