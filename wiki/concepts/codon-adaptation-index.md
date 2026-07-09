---
type: concept
title: "Codon Adaptation Index (CAI)"
tags: [annotation, algorithm]
sources:
  - docs/Evidence/CODON-CAI-001-Evidence.md
  - docs/algorithms/Codon_Optimization/CAI_Calculation.md
source_commit: 9ce49bade5c11e63eebbf8c06dd642662321d5a2
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: codon-cai-001-evidence
      evidence: "ID: CODON-CAI-001 ... Algorithm: CAI Calculation"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:relative-synonymous-codon-usage
      source: codon-cai-001-evidence
      evidence: "w_i = f_i / max(f_j) where i,j in synonymous codons for amino acid — the per-family relative-adaptiveness weight is the same synonymous-family normalization RSCU formalizes"
      confidence: high
      status: current
---

# Codon Adaptation Index (CAI)

A single **directional codon-usage-bias score** for a whole protein-coding gene, defined by
**Sharp & Li (1987)** as "the most widespread technique for analyzing codon usage bias" and used
to **predict a gene's expression level** from its codon choices. Where
[[relative-synonymous-codon-usage|RSCU]] reports a per-codon over/under-representation, CAI
collapses an entire gene to one number in **[0, 1]** measuring how closely its codons match a
**reference set of (ideally highly expressed) genes**. Validated as
[[codon-cai-001-evidence|CODON-CAI-001]]; see [[test-unit-registry]] for how the unit is tracked.

## The measure

CAI is the **geometric mean of the relative adaptiveness `w_i`** of the gene's codons. For codon
`i`, `w_i` is the codon's frequency divided by the frequency of the **most-used synonymous codon**
for the same amino acid — the same synonymous-family normalization that
[[relative-synonymous-codon-usage|RSCU]] formalizes, but referenced to the family **max** rather
than the family mean:

    w_i = f_i / max_j f_j        (i, j synonymous codons of one amino acid)

    CAI = (∏_{i=1}^{L} w_i)^{1/L} = exp( (1/L) · Σ ln w_i )

where `L` is the number of scored codons. The optimal codon of each family has `w = 1`, so
**CAI = 1** when every codon is optimal and **CAI → 0** as rare codons accumulate.

**Geometric mean, not arithmetic:** the product/`exp(mean-ln)` form makes CAI **sensitive to low
values** — a single rare codon (small `w`) pulls the whole score down, which arithmetic averaging
would mask. This is the defining property that separates CAI from a mean-RSCU statistic.

## Scored-codon rules

- **Stop codons excluded** — they encode no amino acid (Sharp & Li 1987).
- **Single-codon amino acids (Met/AUG, Trp/UGG)** have `w ≡ 1` regardless of any bias, because
  their one codon is trivially the family maximum. The **canonical rule (Sharp & Li 1987; Jansen
  et al. 2003, quoted verbatim in the artifact) is to EXCLUDE them** — otherwise a Met/Trp-rich
  gene gets an inflated CAI despite having no codon bias. Seqeron exposes both conventions:
  `CalculateCAI(seq, table)` **includes** Met/Trp (`w=1`, historical behaviour) while
  `CalculateCAI(seq, table, excludeSingleCodonAminoAcids: true)` **excludes** them per the
  canonical definition. Excluding can leave *zero* scored codons (a Met+Trp-only sequence → L=0 →
  CAI = 0, versus 1.0 inclusive).
- **Empty input → CAI = 0** by convention (no codons to evaluate).

## The zero-frequency handling (a key deviation)

The base definition takes `w = 0` for any codon absent from the reference set, forcing `ln(0)` and
CAI = 0. Seqeron's one documented deviation from strict Sharp & Li (1987): when a codon's frequency
is 0 **but its amino acid has other codons present** in the table (`maxFreq > 0`), `w` is **clamped
to `1e-6`** rather than 0 — an incomplete-table protection so one missing entry cannot zero out the
whole gene. If the amino acid is unknown or the family max is 0, `w = NaN` and the codon is skipped
by the caller.

> **Cross-page nuance / follow-up.** [[relative-synonymous-codon-usage]] states that CAI "adds a
> **0.5 pseudocount** for unobserved codons (Sharp & Li 1987)". That 0.5 is Sharp & Li's
> reference-*table*-building convention; the Seqeron **implementation** documented here instead
> uses a **`1e-6` clamp** at score time. Both avoid `log(0)`, but they are different values applied
> at different stages — the RSCU page's phrasing should be read as "a pseudocount-style guard", not
> as "Seqeron uses 0.5 in CAI". Not a source contradiction (the RSCU claim cites the literature,
> this cites the code), but worth reconciling if the two pages are ever merged.

## Worked oracles (E. coli K12, Kazusa species 316407)

- `AUG` → `w = 1.0` → CAI = 1.0 (single optimal codon).
- `CUG-CCG-ACC` (all family-optimal) → CAI = (1·1·1)^{1/3} = 1.0.
- `CUA-CCA-ACA` (all suboptimal; `w` = 0.08 / 0.3585 / 0.2955) → CAI = 0.1980.
- Exclusion mode: `AUGCUACUA` inclusive = 0.18566…, exclusive (drop Met) = **0.08**;
  `AUGUGGCUA` inclusive = 0.43088…, exclusive (drop Met+Trp) = **0.08**; `CUGCUA` (no
  single-codon AA) = 0.28284… under either flag (no effect).

## Place in the codon-usage family

CAI sits one level above [[relative-synonymous-codon-usage|RSCU]] in the codon-usage-bias family:
both count synonymous-family usage, but RSCU normalizes to the family and stays per-codon, while
CAI normalizes to the family **maximum** and reduces a gene to a single geometric-mean expression
proxy. The **[[effective-number-of-codons|ENC / Nc]]** (Wright 1990) is CAI's reference-free
counterpart — same single-number goal, but measuring intrinsic synonymous-codon evenness on a
[20, 61] scale with no reference gene set. **[[codon-optimization]]** is the family's *rewriting*
operation — its MaximizeCAI strategy drives a gene's CAI toward 1 by synonymous substitution.
**[[rare-codon-analysis]]** localizes the low-`w` codons (and their clusters) that pull CAI down.
Other siblings still in `docs/Evidence/` (raw frequency/usage tables) share the same
synonymous-family counting.
See [[algorithm-validation-evidence]] for the shared evidence-artifact pattern.
