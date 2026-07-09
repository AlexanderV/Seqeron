---
type: source
title: "Evidence: META-TAXA-001 (significant taxa detection — Mann–Whitney U / Wilcoxon rank-sum)"
tags: [validation, metagenomics]
doc_path: docs/Evidence/META-TAXA-001-Evidence.md
sources:
  - docs/Evidence/META-TAXA-001-Evidence.md
source_commit: b8447d68c5661a777dae965c591071beb8225772
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: META-TAXA-001

The validation-evidence artifact for test unit **META-TAXA-001** — **significant / differentially-
abundant taxa detection** by the **Mann–Whitney U (Wilcoxon rank-sum)** test with the asymptotic
normal approximation, `MetagenomicsAnalyzer.MannWhitneyU` / `FindSignificantTaxa`. The
community-differential-abundance step downstream of profiling: `FindSignificantTaxa` consumes
per-sample [[taxonomic-profile|taxon→abundance]] maps for two labelled groups. One instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the method is
synthesized in its own concept [[significant-taxa-detection]]; [[test-unit-registry]] tracks the
unit. See `docs/Evidence/META-TAXA-001-Evidence.md`.

## What this file records

- **Online sources (mutually consistent):**
  - **Wikipedia — Mann–Whitney U test** (rank 4, citing primary Mann & Whitney 1947, rank 1) —
    `U1 = R1 − n1(n1+1)/2`, `U1 + U2 = n1·n2`, mean `m_U = n1·n2/2`, `σ_U = sqrt(n1·n2·(n1+n2+1)/12)`,
    the tie-corrected σ, `z = (U − m_U)/σ_U`, midrank tie handling, and the tortoise/hare worked
    example.
  - **SciPy `scipy.stats.mannwhitneyu`** (rank 3, reference implementation) — `U2 = nx·ny − U1`,
    asymptotic μ/σ with tie correction, **continuity correction reduces `|U − m_U|` by 0.5**
    (`use_continuity=True` default for `method='asymptotic'`), and the reference output example.
  - **Xia & Sun (2017), *Genes & Diseases*** (rank 1, PMC6128532) — the domain justification: the
    Wilcoxon rank-sum / Mann–Whitney test is "used to identify the statistically significant
    differences in microbial taxa or OTUs" between two sample groups.
  - **Abramowitz & Stegun 7.1.26 erf approximation** (rank 2, via johndcook.com) — the erf
    polynomial the repository `StatisticsHelper.Erf` / `NormalCDF` uses (constants identical);
    max error |ε| ≤ 1.5×10⁻⁷ ⇒ p-values accurate to ≈1×10⁻⁶.

## Algorithm (from the Evidence file)

Pool `n1 + n2 = n` observations, assign **midranks** (tied values get the average of their
positions), accumulate `Σ(t_k³ − t_k)` over tie blocks. With `R1` = group-1 rank sum:
`U1 = R1 − n1(n1+1)/2`, `U2 = n1·n2 − U1`; `m_U = n1·n2/2`;
`σ_U = sqrt(n1·n2·(n1+n2+1)/12 − n1·n2·Σ(t_k³−t_k)/(12·n·(n−1)))`;
`z = (|U − m_U| − cc)/σ_U` on `U = max(U1,U2)`, `cc ∈ {0, 0.5}`; two-tailed
`p = 2·(1 − Φ(z))` clamped to [0,1], with **σ = 0 → p = 1**. `Φ` = `StatisticsHelper.NormalCDF`
(A&S 7.1.26).

Two methods: `MannWhitneyU(group1, group2, useContinuityCorrection=true)` returns U1/U2/z/p;
`FindSignificantTaxa(profiles, groups, pThreshold=0.05, useContinuityCorrection=true)` runs it per
taxon over the two label-defined groups and returns `SignificantTaxon` records ascending by p-value.

## Source-verified invariants and oracles

**Invariants (INV-01..06):** `U1 + U2 = n1·n2`; `0 ≤ U ≤ n1·n2`; `p ∈ [0,1]`;
`Significant ⇔ PValue < pThreshold`; all-tied → σ→0 → z=0, **p=1** (not significant);
swapping group1/group2 leaves p unchanged (z symmetric via `max(U1,U2)`).

**Datasets / oracles:**
- **SciPy reference:** x=[19,22,16,29,24] (n1=5), y=[20,11,17,12] (n2=4) → `U1=17, U2=3`,
  `m_U=10`, `σ_U=sqrt(200/12)=4.08248290463863`; `z(no-cc)=1.7146428199482247` → `p≈0.0864107`;
  `z(cc)=1.5921683328090657` → `p≈0.11134688653314041` (within erf tolerance).
- **Mann–Whitney tortoise & hare:** n1=n2=6; tortoise rank sum `R_T=32`;
  `U_T = 32 − 21 = 11`, `U_H = 25`, `U_T+U_H = 36 = n1·n2`; `σ_U = sqrt(39) = 6.244997998398398`.

## Documented corner cases / failure modes

- **Ties** (ubiquitous in abundance data): midranks + tie-corrected σ; without correction σ is
  overstated and p is conservative.
- **Identical groups / all-tied:** corrected σ → 0, `z` undefined → report **p = 1**.
- **Small samples:** the normal approximation is asymptotic; an exact/permutation method would
  differ at very small n (SciPy offers one; this impl uses the asymptotic path).
- **Multiple testing:** each taxon tested independently; raw p-values are **not** FDR-corrected here
  (caller applies e.g. Benjamini–Hochberg).
- **Validation:** null group → `ArgumentNullException`; empty group / label ∉ {1,2} / length
  mismatch / missing group → `ArgumentException`; empty `profiles` → empty result; absent taxon → 0.

## Assumptions

Three explicit ASSUMPTIONs, all reference-backed, not invented: (1) **continuity correction on by
default** (SciPy `use_continuity=True`, toggle exposed); (2) **two-tailed** alternative
`p = 2·SF(|z|)` — the default for "is taxon differentially abundant"; (3) **exactly two group
labels**, taxa absent from a profile contribute abundance 0 (ASM-03). No source contradictions; the
only simplifications are the asymptotic-not-exact p-value and the A&S-7.1.26 Φ numerics.
