---
type: concept
title: "Variant allele frequency (empirical VAF + Wilson binomial CI + purity/ploidy correction)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-VAF-001-Evidence.md
source_commit: 68661290101fe6f70f2c89ccf5f5076fff5940ce
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-vaf-001-evidence
      evidence: "Test Unit ID: ONCO-VAF-001 ... Algorithm: Variant Allele Frequency Analysis (empirical VAF, Wilson score confidence interval, purity/ploidy correction)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:somatic-variant-calling-tumor-normal
      source: onco-vaf-001-evidence
      evidence: "Empirical VAF = altReads/totalReads (GATK AD-based, alt AD / sum AD) is the model-free per-locus primitive; the ONCO-SOMATIC-001 caller consumes this same VAF = altReads/totalReads to compare f_t vs f_n."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:tumor-purity-from-mutation-vaf
      source: onco-vaf-001-evidence
      evidence: "AdjustVAFForPurity inverts the CNAqc expected-VAF formula v = mπ/(2(1−π)+π·n_tot) to m·CCF = VAF·(2(1−π)+π·n_tot)/π — the same CNAqc/Tarabichi generative model the ONCO-PURITY-001 estimator inverts the other way for π."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:cancer-cell-fraction-clonal-clustering
      source: onco-vaf-001-evidence
      evidence: "The AdjustVAFForPurity output m·CCF = VAF·(2(1−π)+π·n_tot)/π is exactly the per-tumour-copy mutant fraction numerator of the CCF closed form CCF = f·(ρ·N_T+2(1−ρ))/(ρ·m); dividing by multiplicity m gives CCF."
      confidence: high
      status: current
---

# Variant allele frequency (empirical VAF + Wilson binomial CI + purity/ploidy correction)

The **foundational VAF primitive** of the Oncology family — the **thirty-sixth ingested ONCO-\*
unit** — bundling the three model-free quantities every VAF-based tumor method sits on: the
**empirical variant allele frequency** `VAF = altReads / totalReads`, its **Wilson score binomial
confidence interval**, and the **CNAqc purity/ploidy correction** `AdjustVAFForPurity`. Validated
under test unit **ONCO-VAF-001**; the literature-traced record is [[onco-vaf-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern. Research-grade ([[scientific-rigor|research-grade]]), **not for clinical
or diagnostic use**.

## 1. Empirical VAF — the model-free primitive (GATK AD-based)

```
VAF = altReads / totalReads          # = alt AD / sum AD (samtools mpileup / GATK AD field)
```

The **empirical** allele fraction is the simple, model-free ratio of alt-supporting reads to total
reads at a locus. It is deliberately **not** Mutect2's `AF` FORMAT field, which is a *model estimate*
(a posterior mean marginalised over allele fractions), not the raw ratio — the empirical VAF matches
the textbook `AD`-based definition (GATK Mutect2 FAQ). This is the exact same `VAF = altReads/totalReads`
quantity the somatic caller [[somatic-variant-calling-tumor-normal]] compares between tumor and normal,
that [[ctdna-detection-and-tumor-fraction]] averages across reporters, and that the purity / CCF layers
invert — this unit is where it is defined and bounded.

Oracles: `25/100 → 0.25`, `50/100 → 0.50`, `5/20 → 0.25`, `0/10 → 0.00`, `10/10 → 1.00`.

## 2. Wilson score binomial confidence interval (Wilson 1927)

Treating the `altReads`-of-`totalReads` observation as a binomial proportion `p̂ = VAF` on `n = totalReads`
trials, the **Wilson score interval** gives a bounded confidence interval — the genuinely distinct
statistical content of this unit, not represented elsewhere in the wiki:

```
center = (p̂ + z²/2n) / (1 + z²/n)
margin = (z / (1 + z²/n)) · √( p̂(1−p̂)/n + z²/(4n²) )
CI     = center ± margin
```

with **z = 1.96** for the default 95% interval (`z₀.₀₅`, source-cited verbatim — see assumptions).
Unlike the **Wald** normal-approximation interval `p̂ ± (z/√n)·√(p̂(1−p̂))` — criticised for *overshoot*
(bounds outside [0,1]) and *zero-width* intervals at extreme proportions — the Wilson interval **stays
within [0,1]** and keeps **non-zero width** even at `p̂ = 0` or `p̂ = 1`. This is the load-bearing
property: it makes the interval usable for the low-VAF subclonal / ctDNA regime where Wald fails.

Oracles (z = 1.96), `(alt, total) → center [lower, upper]`:

| alt | total | center | lower | upper |
|-----|-------|--------|-------|-------|
| 25 | 100 | 0.2592 | 0.1755 | 0.3430 |
| 50 | 100 | 0.5000 | 0.4038 | 0.5962 |
| 5 | 20 | 0.2903 | 0.1119 | 0.4687 |
| 0 | 10 | 0.1388 | **0.0000** | 0.2775 |
| 10 | 10 | 0.8612 | 0.7225 | **1.0000** |

The last two rows are the **no-overshoot** guarantees: at `p̂ = 0` the lower bound is exactly 0, at
`p̂ = 1` the upper bound is exactly 1. Invariant `lower ≤ center ≤ upper`, all in [0,1].

## 3. Purity/ploidy correction — `AdjustVAFForPurity` (CNAqc inversion)

The CNAqc (Antonello et al., *Genome Biology* 2024) expected-VAF generative model for a mutation of
multiplicity `m` on a segment of total copy number `n_tot` in a tumour of purity `π` is
`v = m·π / [2(1−π) + π·n_tot]` (normal ploidy fixed at 2). Solving for the **per-tumour-copy mutant
fraction** `m·CCF` gives the correction this unit computes:

```
adjusted = m·CCF = VAF · (2(1−π) + π·n_tot) / π
```

This is the **same CNAqc/Tarabichi generative model** that [[tumor-purity-from-mutation-vaf]] inverts
for π and that [[cancer-cell-fraction-clonal-clustering]] divides by `m` to recover CCF — here it
returns the copy-number-corrected mutant fraction directly. Oracles:

| VAF | π | n_tot | adjusted (m·CCF) |
|-----|-----|-------|-------------------|
| 0.40 | 0.80 | 2 | 1.00 |
| 0.35 | 0.70 | 2 | 1.00 |
| 0.20 | 0.50 | 2 | 0.80 |
| 0.30 | 0.50 | 4 | 1.80 |

(The first two rows are the diploid-heterozygous clonal band: at `n_tot = 2` the correction reduces to
`2·VAF·(1−π)/π + VAF`, and a clonal het mutation at expected `VAF = π/2` recovers `m·CCF = 1`.)

## Corner cases and failure modes

- **VAF > 1 from alignment artifacts** — misaligned / overlapping reads can push raw alt counts above
  the locus total; the empirical definition requires `altReads ≤ totalReads`, so a count violating this
  is **invalid input** (`ArgumentOutOfRangeException`); negative counts likewise throw.
- **Zero coverage (`totalReads = 0`)** — VAF is undefined (0/0); the division-by-zero must be guarded
  (defined as 0 / no observation).
- **Extreme proportions (`p̂ = 0` or `1`)** — the Wilson interval stays in [0,1] with non-zero width
  (Wald would overshoot / collapse) — the reason Wilson is used, not Wald.
- **Purity `π = 0`** — `AdjustVAFForPurity` divides by π; pure-normal `π = 0` makes the adjusted
  fraction undefined (`ArgumentOutOfRangeException`).

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference. Every formula is source-derived: the
empirical VAF from the GATK AD-based definition, the Wilson interval verbatim from Wilson 1927 (via the
Wikipedia binomial-proportion-CI page, citing the primary), and the purity correction from the CNAqc
*Genome Biology* 2024 expected-VAF formula (Tarabichi 2017 diploid case corroborating). **Two flagged
assumptions**, both source-backed rather than free: (1) **z = 1.96** for the default 95% interval (the
source states `z₀.₀₅ = 1.96`; used verbatim, not the more precise 1.959964, for traceability), and
(2) `AdjustVAFForPurity` fixes **normal copy number = 2** and assumes a heterozygous-eligible
diploid-normal background (CNAqc/Tarabichi autosomal-diploid convention; sex chromosomes / non-diploid
normals out of scope). **No source contradictions** — GATK (empirical VAF), Wilson 1927 (the CI), and
CNAqc / Tarabichi (the purity correction) cover disjoint facets and agree. **Not for clinical or
diagnostic use.**
