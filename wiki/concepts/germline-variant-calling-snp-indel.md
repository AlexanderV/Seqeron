---
type: concept
title: "Variant calling (SNP/indel from reference↔query alignment + Ti/Tv)"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/VARIANT-CALL-001-Evidence.md
source_commit: 5b4dd805db54d51bae30445a884e122fc4d97bd5
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: variant-call-001-evidence
      evidence: "Test Unit ID: VARIANT-CALL-001 — Variant Detection (SNP / Insertion / Deletion calling from a reference↔query comparison, with transition/transversion classification)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:variant-effect-annotation-vep
      source: variant-call-001-evidence
      evidence: "Calling produces the variant that annotation interprets: Danecek 2011 — 'a variant is a difference from reference'; VARIANT-ANNOT-001 takes an already-called variant and predicts its consequence, so calling is the upstream step feeding annotation"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:somatic-variant-calling-tumor-normal
      source: variant-call-001-evidence
      evidence: "Both detect variants but by different evidence: this unit calls germline SNP/indel from a reference↔query global alignment (CallVariantsFromAlignment); ONCO-SOMATIC-001 classifies Somatic/Germline from tumor-vs-matched-normal VAF thresholds"
      confidence: high
      status: current
---

# Variant calling (SNP/indel from reference↔query alignment + Ti/Tv)

**Variant calling** is the **detection** step of the variant-analysis family: it *produces* the variants
every downstream unit interprets. This unit calls **SNPs, insertions, and deletions** (Danecek 2011 — a
variant is a *difference from reference*) by aligning a query sequence against a reference and reading out
the differing columns, then classifies each substitution as a **transition** or **transversion** and
reports the **Ti/Tv** ratio. `CallVariantsFromAlignment` runs `SequenceAligner.GlobalAlign(reference,
query)` and walks the aligned columns, emitting one `Variant` per difference. Validated as test unit
**VARIANT-CALL-001** ([[variant-call-001-evidence]]); see [[test-unit-registry]] for how units are tracked
and [[algorithm-validation-evidence]] for the evidence-artifact pattern. Research-grade
([[research-grade-limitations]]), **not for clinical or diagnostic use**.

## Three variant classes from aligned columns

Walking the reference/query alignment column by column yields exactly three difference types
(Danecek 2011 — the VCF variant classes):

| Column | Class | In-memory representation |
|--------|-------|--------------------------|
| both bases present, differ | **SNP** | ref base, alt base, 0-based `Position` |
| gap in reference, base in query | **Insertion** | ref = `"-"` gap sentinel, alt base |
| base in reference, gap in query | **Deletion** | ref base, alt = `"-"` gap sentinel |

Identical sequences yield **zero variants**. Mismatched aligned lengths throw `ArgumentException`; empty
input yields an empty call set.

## Transition / transversion classification

Every **SNP** is classified case-insensitively as a **transition** (a base change *within* a ring class —
purine↔purine `A↔G` or pyrimidine↔pyrimidine `C↔T`) or a **transversion** (*across* ring classes —
`A↔C`, `A↔T`, `G↔C`, `G↔T`). There are **2 possible transitions but 4 possible transversions** per base,
yet transitions occur *more often* — the **transitional bias** (Collins & Jukes 1994: transition rate
1.71×10⁻⁹ > transversion rate 1.22×10⁻⁹; ≈2/3 of SNPs are transitions). The summary statistic is the

```
Ti/Tv = (#transitions) / (#transversions)
```

**Convention:** the mathematically-undefined zero-transversion case (`#Tv = 0`) is mapped to **0**, not
`+∞` or an exception. Insertions/deletions are not SNPs, so they classify as `Other` (Ti/Tv is defined
only for substitutions). A genome-wide Ti/Tv well above 0.5 (the naive 2:4 expectation) is the standard
sanity check for a real vs artefactual SNP call set.

## Where it sits in the variant-analysis family

Calling is the **head** of the germline variant pipeline. Its output feeds
[[variant-effect-annotation-vep]] (VARIANT-ANNOT-001), which takes each *already-called* variant and
predicts its functional consequence (missense / stop_gained / frameshift …) plus Sequence-Ontology
IMPACT — annotation cannot run without a caller upstream. This unit is the **germline, reference↔query**
caller; the oncology sibling [[somatic-variant-calling-tumor-normal]] (ONCO-SOMATIC-001) is
**`alternative_to`** it — same goal (find variants), different evidence: tumor-vs-matched-normal **VAF**
thresholds rather than a pairwise alignment. Serialized output (`ToVcfLines`) targets the **VCFv4.3**
format; VEP-style annotation and the oncology tier/pathogenicity layers consume such calls.

## Scope and assumptions

- **Internal gap-sentinel indel representation** — the in-memory `Variant` uses `"-"` and a **0-based**
  `Position`; the VCF **padding base** and **1-based POS** conventions (VCFv4.3 field 4) apply only to the
  *serialized* `ToVcfLines` output, not the in-memory model. E.g. the VCFv4.3 §1.1 examples `GTC→G`
  (2-base deletion) and `GTC→GTCT` (1-base insertion) both anchor at the preceding `G` — a serialization
  concern, not a detection concern.
- **Indels are not left-aligned / parsimony-normalized** — per Tan 2015 the canonical representation
  requires **left-alignment + parsimony**; this caller reports the indel at the column `GlobalAlign`
  produces, with no normalization pass. This affects reported **position** in low-complexity / repeated
  regions only (the alignment there is non-unique), **not** variant counts or types.
- **Alignment-based, not pileup-based** — this is a single-query-vs-reference comparison, not a
  read-pileup genotype caller; there is no depth model, base-quality weighting, or diploid genotype
  assignment.

Reference sources — **VCFv4.3** (samtools/hts-specs, the REF/ALT/POS + padding grammar), **Danecek 2011**
(the VCF variant classes), **Tan 2015** (left-align + parsimony normalization), **Collins & Jukes 1994**
(transitional bias), and Wikipedia Transition/Transversion (the classification table) — full record in
[[variant-call-001-evidence]]. No source contradictions.
