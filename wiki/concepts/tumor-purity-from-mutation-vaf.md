---
type: concept
title: "Tumor purity estimation from somatic SNV VAF (CNAqc expected-VAF inversion)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-PURITY-001-Evidence.md
source_commit: fdf583e25989b1d2bcbc999fa056fb16119f8c31
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-purity-001-evidence
      evidence: "Test Unit ID: ONCO-PURITY-001 ... Algorithm: Tumor Purity Estimation (from somatic SNV VAF / allele-specific copy number)"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:allele-specific-copy-number-ascat
      source: onco-purity-001-evidence
      evidence: "Both estimate tumour purity π but differently: this unit inverts the CNAqc expected-VAF formula v = mπc/[2(1−π)+π(n_A+n_B)] from a somatic-SNV VAF (closed-form), whereas ASCAT (ONCO-ASCAT-001) fits ρ jointly with ploidy over a logR/BAF grid. FACETS 2016 confirms the shared 2(1−π)+π·n_tot denominator."
      confidence: high
      status: current
---

# Tumor purity estimation from somatic SNV VAF (CNAqc expected-VAF inversion)

A **standalone tumor-purity estimator** that inverts the somatic-mutation **expected-VAF generative
model** to recover purity **π (= ρ)** from an observed variant allele frequency and the local
copy-number state — the **ABSOLUTE / CNAqc / FACETS** family of purity calls, computed as a
**closed form** from clonal mutations rather than a grid search. Validated under test unit
**ONCO-PURITY-001**; the literature-traced record is [[onco-purity-001-evidence]],
[[test-unit-registry]] tracks the unit, and [[algorithm-validation-evidence]] describes the
evidence-artifact pattern.

## Distinct from the ASCAT joint fit

This is **not** the [[allele-specific-copy-number-ascat|ASCAT joint purity/ploidy grid fit]]
(ONCO-ASCAT-001). ASCAT *fits* ρ and ψ **jointly from raw per-locus logR/BAF tracks** by minimising a
length-weighted integer-closeness objective over a (ρ, ψ) grid. This unit instead reads an observed
somatic-SNV **VAF** plus the local allele-specific copy number and **inverts one closed-form equation**
for π — the same purity scalar, a different substrate (mutation VAFs, not logR/BAF) and a different
method (algebraic inversion, not a grid). It is the mutation-VAF counterpart to the ploidy-side
[[tumor-ploidy-estimation-and-whole-genome-doubling]], and it shares the copy-number-corrected
diploid-heterozygous VAF inversion `π = 2·VAF` with [[ctdna-detection-and-tumor-fraction]]'s tumor
fraction and with the CCF/multiplicity inversion carried in [[allele-specific-copy-number-ascat]] §4
and [[cancer-cell-fraction-clonal-clustering]].

## 1. The expected-VAF generative model (CNAqc)

For a somatic mutation present in **m** copies of the tumour genome (multiplicity), at cancer-cell
fraction (clonality) **c**, on a segment with allele-specific copy numbers **n_A:n_B** (total
n_tot = n_A + n_B), CNAqc (Antonello et al., *Genome Biology* 2024) gives the expected VAF:

```
v = m·π·c / [ 2(1−π) + π·(n_A + n_B) ]
```

The denominator is the mean allele-copy count per cell in the tumour/normal mixture: `2(1−π)` is the
**healthy diploid normal** contribution (2 copies weighted by 1−π) and `π·n_tot` is the tumour
contribution. **FACETS** (Shen & Seshan, *NAR* 2016) independently confirms this structure — it mixes a
normal `(1,1)` genotype with an aberrant `(m,p)` genotype at cellular fraction Φ via `m* = mΦ + (1−Φ)`,
the same `2(1−π) + π·n_tot` denominator. **ABSOLUTE** (Carter et al., *Nat Biotechnol* 2012) is the
inverse-direction sibling: it converts allelic fractions into per-cancer-cell multiplicity by
correcting for purity and local copy number.

## 2. Inversion for purity

Solving the model for π (clonal `c` and multiplicity `m` known, n_tot from the segment):

```
π = 2v / ( m·c + 2v − v·n_tot )
```

**`EstimatePurityFromVAF` — the copy-neutral diploid heterozygous special case.** For a **clonal**
(c = 1), **heterozygous** (m = 1) somatic SNV on a **copy-neutral diploid** (n_tot = 2) locus this
reduces to the textbook closed form:

```
π = 2 · VAF                    # since v = π/2  at m=1, c=1, n_tot=2
```

This is the robust special case CNAqc's worked band uses, and it is stated explicitly in the API
contract. Invariant `0 ≤ π ≤ 1` (so a diploid-model VAF > 0.5 implying π > 1 is rejected).

**`EstimatePurity` — the general allele-specific inversion.** Given the segment copy number n_tot and
multiplicity m, invert the full formula (`π = 2v/(m·c + 2v − v·n_tot)`) — this recovers purity on
amplified/LOH segments where the diploid shortcut does not hold.

## 3. Worked oracles

| Case | Params | Expected VAF | Recovered π |
|------|--------|--------------|-------------|
| Clonal het diploid | n_tot=2, m=1, c=1 | 0.30 | **0.60** (= 2·VAF) |
| Boundary | m=1, c=1, n_tot=2 | 0.50 / 0.00 | 1.0 / 0.0 |
| Tolerance band | purity 0.55–0.65 | VAF 0.275–0.325 | — |
| 2:1 amplified, π=1 | n_tot=3, m=1, c=1 | 1/3 ≈ 0.333 | **1.0** |
| 2:1 amplified, π=1 | n_tot=3, m=2, c=1 | 2/3 ≈ 0.667 | **1.0** |

(CNAqc *Genome Biology* 2024 worked examples: a real 60% purity → 30% VAF for a diploid het mutation;
a 2:1 segment shows two clonal peaks at 33% and 66% because multiplicity m may be 1 or 2.)

**Aggregation.** When several clonal heterozygous SNVs are supplied, per-variant estimates are combined
by the **median** (robust to subclonal/outlier VAFs). This is a documented, non-correctness-affecting
aggregation policy over the source-derived single-variant formula, recorded as an assumption.

## Corner cases and failure modes

- **Multiplicity ambiguity on amplified segments:** on a 2:1 (n_tot=3) segment clonal mutations form
  two VAF peaks (1/3 and 2/3) because m may be 1 or 2 — purity cannot be inferred from VAF alone without
  the copy-number state and m. Copy-neutral diploid heterozygous (1:1) loci avoid this (m = 1,
  n_tot = 2), which is why `π = 2·VAF` is the robust closed form there.
- **Subclonal mutations (c < 1):** VAF is depressed (v ∝ c); treating a subclonal VAF as clonal
  **underestimates** purity. Purity must be estimated from **clonal** mutations.
- **Purity < 0.1 (below detection):** very low purity yields VAFs near sequencing noise; the estimator
  should return a small π near 0 without error. High stromal contamination is equivalent to low purity.
- **No informative variants:** with no usable heterozygous SNPs / clonal variants, purity is undefined.
- **Invalid inputs** (VAF outside [0,1]; a diploid-model VAF > 0.5 implying π > 1; empty variant list;
  non-positive copy number) are rejected.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for a closed-form purity estimator. The
formula is fully source-derived (CNAqc verbatim, FACETS + ABSOLUTE corroborating); the only modelling
scope choices are the copy-neutral-diploid default for `EstimatePurityFromVAF` and the median
cross-variant aggregation — neither invents a numeric constant. **Not for clinical or diagnostic use.**
No source contradictions: CNAqc (the expected-VAF formula + worked peaks), FACETS (the mixing-model
denominator), and ABSOLUTE (the inverse purity/copy-number correction) agree.
