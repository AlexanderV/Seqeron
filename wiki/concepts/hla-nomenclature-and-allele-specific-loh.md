---
type: concept
title: "HLA allele nomenclature parsing + allele-specific HLA LOH (LOHHLA)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-HLA-001-Evidence.md
  - docs/Evidence/ONCO-IMMUNE-001-Evidence.md
source_commit: a197fb86ceeffb8de5c09005d269f020e46584f5
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-hla-001-evidence
      evidence: "Test Unit ID: ONCO-HLA-001 ... Algorithm: HLA allele nomenclature parsing/validation + allele-specific HLA LOH (LOHHLA) classification"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:allele-specific-copy-number-ascat
      source: onco-hla-001-evidence
      evidence: "LOHHLA infers, for each HLA gene, the allele-specific copy number of both homologous alleles from tumor coverage relative to germline (logR) and B-allele frequencies (BAF) at polymorphic sites — the HLA-locus specialization of the genome-wide allele-specific copy-number / LOH derivation."
      confidence: high
      status: current
---

# HLA allele nomenclature parsing + allele-specific HLA LOH (LOHHLA)

The Oncology family's **immuno-oncology / antigen-presentation** unit (**ONCO-HLA-001**), covering two
disjoint pieces that share the HLA locus as their subject:

1. **HLA allele nomenclature parsing/validation** — parse a WHO IPD-IMGT/HLA allele name into its
   gene + colon-delimited field tuple + expression suffix, and validate it.
2. **Allele-specific HLA loss-of-heterozygosity (LOH) classification** — the **LOHHLA**
   (McGranahan et al. 2017) decision rule that calls HLA LOH from per-allele copy number and an
   allelic-imbalance test.

The literature-traced record is [[onco-hla-001-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern. HLA LOH is a mechanism of
**immune escape** — a tumour that loses an HLA allele can no longer present the neoantigens restricted
to it, so this unit sits in the antigen-presentation / neoantigen area of the Oncology family. The
affinity gate that decides which peptides an HLA allele presents (and thereby which neoantigens a lost
allele would have shown) is [[mhc-peptide-binding-prediction]] (ONCO-MHC-001). Its
antigen-agnostic quantitative sibling in the immuno-oncology area is
[[immune-infiltration-deconvolution]] (ONCO-IMMUNE-001), which quantifies the immune *cell content* of
the microenvironment rather than the antigen-presentation machinery.

## 1. HLA allele nomenclature (WHO / IPD-IMGT/HLA; Marsh et al. 2010)

An allele name has the form `HLA-[Gene]*[Field1]:[Field2][:Field3][:Field4][suffix]`:

- the **gene** is separated from the numeric fields by an asterisk `*`; numeric fields are separated by
  **colons** (the colon-delimited convention introduced by Marsh 2010);
- **Field 1** = the *type* / allele group (often the serological antigen);
- **Field 2** = the specific **HLA protein** (subtype, non-synonymous coding difference);
- **Field 3** = alleles differing only by **synonymous** coding substitutions;
- **Field 4** = alleles differing only by **non-coding** (intron / 5′ or 3′ UTR) polymorphisms;
- an optional trailing **expression suffix** — one of `N` (Null / not expressed), `L` (Low surface),
  `S` (Secreted / soluble only), `C` (Cytoplasm only), `A` (Aberrant), `Q` (Questionable).

**Validity rules:** every valid allele has **at least two fields** (Field1:Field2 — the "four-digit"
minimum) and **at most four**; the `HLA-` prefix is required; each field is numeric; and any trailing
letter must be from the expression-suffix set. So `HLA-A*02:01` (2 fields), `HLA-B*07:02:01` (3),
`HLA-C*07:02:01:03` (4), and `HLA-A*24:02:01:02L` (4 + `L`) parse; while `HLA-A*02` (one field only),
`A*02:01` (no `HLA-` prefix), `HLA-A*02:01:01:01:01` (five fields), and `HLA-A*02:01X` (`X` not a valid
suffix) are rejected.

## 2. Allele-specific HLA LOH — LOHHLA (McGranahan et al. 2017)

LOHHLA infers, for **each HLA gene**, the **allele-specific copy number of both homologous alleles**
from tumor coverage relative to germline (**logR**) and **B-allele frequencies (BAF)** at polymorphic
sites that distinguish the two alleles. The reference implementation (`mskcc/lohhla`,
`LOHHLAscript.R`) reports per-allele copy numbers (`HLA_type1copyNum_withBAF` /
`HLA_type2copyNum_withBAF`, and `_withoutBAF` variants) — these are the classifier's input.

**The call rule (two source-exact thresholds, both must hold):**

```
HLA LOH  ⟺  (min allele copy number < 0.5)   AND   (allelic-imbalance paired t-test p < 0.01)
```

- **Loss threshold:** "A copy number < 0.5 … is classified as subject to loss, and thereby indicative
  of LOH" — the allele with copy number **< 0.5** is the *lost* allele.
- **Allelic-imbalance guard:** "To avoid over-calling LOH, we also calculate a p value relating to
  allelic imbalance … Allelic imbalance is determined if **p < 0.01** using the paired Student's
  t-Test." (`PairedTtest <- t.test(...,paired=TRUE)` in the reference R). LOH is called **only** when
  the imbalance is significant, even if a raw copy number dips below 0.5 — this is the explicit
  "avoid over-calling" guard.

Both thresholds are **strict `<`**: copy number **exactly 0.5** is *not* loss, and p **exactly 0.01**
is *not* significant.

### Worked oracles (from the evidence datasets)

| Allele1 CN | Allele2 CN | Imbalance p | HLA LOH? | Lost allele | Why |
|-----------|-----------|-------------|----------|-------------|-----|
| 1.8  | 0.30 | 0.001  | yes | allele 2 | CN < 0.5 and p < 0.01 |
| 0.10 | 1.50 | 0.0005 | yes | allele 1 | CN < 0.5 and p < 0.01 |
| 1.10 | 0.90 | 0.30   | no  | — | both retained (CN ≥ 0.5) |
| 1.60 | 0.40 | 0.05   | no  | — | over-calling guard: p ≥ 0.01 |
| 1.50 | 0.50 | 0.001  | no  | — | boundary: 0.5 is not < 0.5 |
| 1.70 | 0.40 | 0.01   | no  | — | boundary: 0.01 is not < 0.01 |

## Relation to the copy-number family

This is the **HLA-locus specialization** of allele-specific copy-number / LOH inference. The genome-wide
allele-specific copy-number layer [[allele-specific-copy-number-ascat]] fits per-segment integer
`nA/nB` from logR/BAF by a joint purity/ploidy grid; LOHHLA restricts the same logR + BAF allelic-
contrast idea to the (hard-to-map, polymorphic) HLA genes and outputs a per-allele **continuous** copy
number plus a dedicated **allelic-imbalance t-test** rather than an integer segmentation. The two are
complementary rather than alternatives: HLA LOH needs allele-personalised references that a generic
segmenter does not build.

## Corner cases and assumptions

- **No-LOH despite low copy:** if imbalance is not significant (p ≥ 0.01), no LOH even when a raw copy
  number < 0.5 (the over-calling guard).
- **Both alleles retained:** both CN ≥ 0.5 → heterozygous-retained, no LOH.
- **Homozygous locus:** allele-specific loss cannot be assessed when the two homologs are identical (no
  polymorphic sites distinguish them).
- **ASSUMPTION — lost-allele tie-break (both alleles < 0.5):** McGranahan 2017 does not specify
  behaviour when *both* alleles are < 0.5 with significant imbalance (biologically a homozygous
  deletion, not allele-specific LOH). This unit calls LOH only when **exactly one** allele is < 0.5;
  when both are, it reports a distinct **`HomozygousLoss`** label. Only the label is affected — the two
  source-exact thresholds (0.5, 0.01) are unchanged.
- **Nomenclature:** the two-field minimum, four-field maximum, required `HLA-` prefix, numeric fields,
  and the fixed N/L/S/C/A/Q suffix set are the only validity constraints.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for HLA allele-name parsing and the LOHHLA
LOH call. The nomenclature rules are the WHO/IPD-IMGT/HLA standard (rank 1–2 sources); the LOH
thresholds are McGranahan 2017 verbatim, corroborated by the `mskcc/lohhla` reference R. **Not for
clinical or diagnostic use.** No source contradictions — the WHO nomenclature standard and the LOHHLA
paper/implementation cover disjoint parts of the unit.
