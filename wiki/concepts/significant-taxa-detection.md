---
type: concept
title: "Significant taxa detection (Mann–Whitney U / Wilcoxon rank-sum differential abundance)"
tags: [metagenomics, algorithm]
sources:
  - docs/Evidence/META-TAXA-001-Evidence.md
  - docs/algorithms/Metagenomics/Significant_Taxa_Detection.md
source_commit: b8447d68c5661a777dae965c591071beb8225772
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: meta-taxa-001-evidence
      evidence: "Test Unit ID: META-TAXA-001, Algorithm: Significant Taxa Detection (Mann–Whitney U / Wilcoxon rank-sum), Methods MetagenomicsAnalyzer.MannWhitneyU / FindSignificantTaxa"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:taxonomic-profile
      source: meta-taxa-001-evidence
      evidence: "FindSignificantTaxa consumes per-sample taxon→abundance maps (IReadOnlyList<IReadOnlyDictionary<string,double>>) — the community-abundance vectors that GenerateTaxonomicProfile produces; a taxon absent from a profile contributes abundance 0 (ASM-03)."
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:pathway-enrichment-ora
      source: meta-taxa-001-evidence
      evidence: "Both are per-feature significance tests registered under metagenomics, but ORA uses the hypergeometric right-tail over membership counts whereas this uses the rank-sum U statistic over per-sample abundance vectors between two groups — different null model, different question (over-representation vs differential abundance)."
      confidence: medium
      status: current
---

# Significant taxa detection (Mann–Whitney U / Wilcoxon rank-sum)

**Significant / differentially-abundant taxa detection** answers "which microbial taxa differ in
abundance between two sample groups (e.g. case vs control)?" For **each taxon**, the per-sample
abundances of the two groups are compared with the **Mann–Whitney U test** (equivalently the
**Wilcoxon rank-sum test**) — the standard non-parametric approach for microbiome differential
abundance (Xia & Sun 2017). It is rank-based (no normality assumption), which suits compositional,
heavily-tied, zero-inflated abundance data. Validated under test unit **META-TAXA-001**; the record
is [[meta-taxa-001-evidence]], [[test-unit-registry]] tracks the unit, and
[[algorithm-validation-evidence]] describes the artifact pattern.

This is the **community differential-abundance** unit of the Metagenomics family, downstream of
profiling: `FindSignificantTaxa` **consumes the per-sample [[taxonomic-profile|taxon→abundance
maps]]** and asks the between-group question that profiling and the within/between-sample diversity
units ([[alpha-diversity]], [[beta-diversity]]) do not. It is **distinct** from the other
per-feature significance tests in the wiki because of its **statistical test**: the hypergeometric
right-tail of [[pathway-enrichment-ora]] scores set *over-representation*, and the Fisher's-exact
two-sample test of [[differentially-methylated-regions]] compares *methylation* — this one uses the
**rank-sum U statistic with a normal approximation** over two abundance vectors. All three share the
repository's [[epigenetic-age-horvath-clock|StatisticsHelper]] Abramowitz–Stegun 7.1.26 erf →
`NormalCDF` numerics (the same A&S 7.1.26 approximation the DMR/normal-CDF paths use).

## Two entry points

Impl `MetagenomicsAnalyzer.cs` (`Seqeron.Genomics.Metagenomics`):

- **`MannWhitneyU(group1, group2, useContinuityCorrection=true)`** — the core rank-sum test over two
  `IReadOnlyList<double>` observation vectors; returns a `MannWhitneyResult` `(U1, U2, Z, PValue)`.
  `O(n log n)` (dominated by the pooled sort, `n = n1+n2`).
- **`FindSignificantTaxa(profiles, groups, pThreshold=0.05, useContinuityCorrection=true)`** — runs
  the test **per taxon** over the two label-defined groups and returns `SignificantTaxon`
  `(Taxon, U, Z, PValue, Significant)` records **sorted ascending by p-value**. `profiles` is a list
  of per-sample `taxon→abundance` maps; `groups` is an index-aligned list of labels ∈ {1,2}.
  `O(t · s log s)` for `t` taxa, `s` samples.

## The test

Pool the `n1` observations of group 1 and `n2` of group 2 (`n = n1+n2`), sort, and assign **ranks**;
tied values receive the **midrank** (average of the positions they occupy). With `R1` = rank sum of
group 1:

```
U1 = R1 − n1(n1+1)/2      U2 = n1·n2 − U1       (so U1 + U2 = n1·n2)
m_U = n1·n2 / 2
σ_U = sqrt( n1·n2·(n1+n2+1)/12  −  n1·n2·Σ(t_k³ − t_k) / (12·n·(n−1)) )
z   = ( |U − m_U| − cc ) / σ_U        U = max(U1, U2),  cc ∈ {0, 0.5}
p   = 2·(1 − Φ(z))                    two-tailed, clamped to [0,1]
```

- **Tie correction** is the right-hand subtraction in `σ_U`, where `t_k` = number of tied
  observations at distinct value `k`; it is **0 when there are no ties**. Without it σ is overstated
  and p-values are conservative.
- **Continuity correction** subtracts `0.5` from `|U − m_U|`, **on by default** to match SciPy's
  `use_continuity=True` for the asymptotic method (a documented reference default, exposed as a
  toggle). It *raises* the p-value (`p_cc > p_no-cc`).
- **z uses the larger U** / `|U − m_U|`, so the result is **symmetric** in the two groups.
- **Φ** is the standard-normal CDF via `StatisticsHelper.NormalCDF` (A&S erf approximation 7.1.26,
  |ε| ≤ 1.5×10⁻⁷) ⇒ p-values match an exact normal CDF only to **≈1×10⁻⁶** (numerical, not
  algorithmic, imprecision).

## Invariants and edge cases

- **INV-01:** `U1 + U2 = n1·n2` (definition of `U2`).
- **INV-02:** `0 ≤ U ≤ n1·n2`.
- **INV-03:** `PValue ∈ [0, 1]`.
- **INV-04:** `Significant ⇔ PValue < pThreshold`.
- **INV-05 (degenerate):** **all observations tied** (across both groups) ⇒ `σ_U → 0` ⇒ `z = 0`,
  **`p = 1`, not significant** — no evidence against H₀.
- **INV-06:** swapping group1/group2 leaves `p` unchanged (z is symmetric).
- **Absence = 0:** a taxon absent from a sample is treated as **abundance 0** in that sample
  (ASM-03), not missing — standard for abundance tables.
- **Validation:** `MannWhitneyU` → `ArgumentNullException` (null group) / `ArgumentException` (empty
  group). `FindSignificantTaxa` → `ArgumentNullException` (null `profiles`/`groups`),
  `ArgumentException` (length mismatch, label ∉ {1,2}, or a missing group); **empty `profiles` →
  empty result**.

Worked oracles (from [[meta-taxa-001-evidence]]): the **SciPy reference** x=[19,22,16,29,24],
y=[20,11,17,12] → `U1=17, U2=3`, `σ_U=sqrt(200/12)=4.0825`, `z(cc)=1.5922 → p≈0.11135`,
`z(no-cc)=1.7146 → p≈0.08641`; the **tortoise/hare** example → `U_T=11, U_H=25, U_T+U_H=36`.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for the per-taxon two-group
rank-sum test. **Two-group comparison only**; **asymptotic** (normal-approximation) p-values, not
exact/permutation — at very small `n` an exact test would differ slightly (ASM-02). **No
multiple-testing / FDR correction** is applied across taxa — each taxon is tested independently and
correcting the raw p-values (e.g. Benjamini–Hochberg) is the caller's responsibility. The rank-based
per-taxon test **ignores the compositional structure** of microbiome data (a known limitation of
non-parametric per-taxon approaches). The formulas match their primary sources verbatim (Mann &
Whitney 1947; SciPy `mannwhitneyu`; A&S 7.1.26 for Φ); no source contradictions.
