# Clonal vs Subclonal Mutation Classification

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-CLONAL-001 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Classifies each somatic mutation in a tumour as **clonal** (present in essentially all cancer cells) or **subclonal** (present only in a subpopulation), from per-variant read evidence and tumour purity. For each variant a posterior distribution over the cancer cell fraction (CCF) is built from the alternate/total read counts and the local copy-number state, and the mutation is called clonal when the posterior probability that the CCF exceeds 0.95 is greater than 0.5 [1]. The classification is **probabilistic** (it uses the full posterior, not just the point estimate) and **deterministic** given the inputs (no sampling). A point-estimate helper classifies already-computed CCF values by the same 0.95 threshold [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A bulk tumour sample is a mixture of normal cells and cancer cells, and the cancer cells may themselves comprise several subclones. The cancer cell fraction (CCF) of a somatic mutation is the proportion of cancer cells that carry it [2]. Trunk (early, clonal) mutations occur in all cancer cells; branch (later, subclonal) mutations occur in a subset. Distinguishing the two underlies tumour-evolution timing and therapy decisions [1].

### 2.2 Core Model

For a mutation observed in `a` of `N` reads at a locus of absolute somatic copy number `q`, in a sample of purity `α`, the expected alternate-allele fraction of a mutation present in a fraction `c` of cancer cells at multiplicity `M` (mutated copies per cancer cell) is [1][2]:

```
f(c) = α·M·c / (2(1 − α) + α·q)
```

Landau et al. (2013) state this for one mutated copy (M = 1): `f(c) = αc / (2(1 − α) + αq)` [1]; DeCiFering (Satas et al. 2021, Eq. 1) gives the multiplicity-general form `c ≈ (1/ρ)·(ρ·N_tot + 2(1 − ρ))/M · v̂`, the algebraic inverse of the relation above with `N_tot = q` [2].

The posterior over `c` assumes a uniform prior and a binomial read model [1]:

```
P(c) ∝ Binomial(a | N, f(c)),   c ∈ [0.01, 1]
```

evaluated on a regular grid of 100 values of `c` and normalised by its sum [1].

**Classification rule** [1]: a mutation is **clonal** if `P(CCF > 0.95) > 0.5`, and **subclonal** otherwise.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | ClonalCount + SubclonalCount = total | every variant is classified into exactly one of two states |
| INV-02 | ClonalFraction = ClonalCount / total (0 if empty) | definition of the clonal-fraction summary |
| INV-03 | CCF point estimate ∈ [0.01, 1] | posterior is supported on the grid c ∈ [0.01, 1] [1] |
| INV-04 | ProbabilityClonal ∈ [0, 1] | normalised posterior probability mass |
| INV-05 | Clonal ⇔ P(CCF > 0.95) > 0.5 | Landau (2013) classification rule [1] |
| INV-06 | IdentifyClonalMutations selects i ⇔ ccf[i] > 0.95 (strict) | Landau (2013) clonal threshold CCF > 0.95 [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| variants | IEnumerable\<ClonalityVariant\> | required | per-variant read evidence and copy-number state | not null |
| ClonalityVariant.AltReads | int | required | alternate reads `a` | 0 ≤ a ≤ TotalReads |
| ClonalityVariant.TotalReads | int | required | total reads `N` | ≥ 1 |
| ClonalityVariant.LocalCopyNumber | int | required | absolute tumour total copy number `q` | ≥ 1 |
| ClonalityVariant.Multiplicity | int | 1 | mutated copies per cancer cell `M` | 1 ≤ M ≤ LocalCopyNumber |
| purity | double | required | tumour purity `α` | ρ ∈ (0, 1] |
| ccfValues | IEnumerable\<double\> | required | pre-computed CCF point estimates (IdentifyClonalMutations) | each ∈ [0, 1] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| ClonalityResult.Calls | IReadOnlyList\<ClonalityCall\> | per-variant calls, in input order |
| ClonalityCall.Ccf | double | posterior-mean CCF estimate ∈ [0.01, 1] |
| ClonalityCall.ProbabilityClonal | double | P(CCF > 0.95) |
| ClonalityCall.Status | ClonalityStatus | Clonal / Subclonal |
| ClonalityResult.ClonalCount / SubclonalCount | int | counts per class |
| ClonalityResult.ClonalFraction | double | ClonalCount / total (0 if empty) |
| IdentifyClonalMutations return | IReadOnlyList\<int\> | 0-based indices of clonal CCF values |

### 3.3 Preconditions and Validation

`variants` / `ccfValues` null → `ArgumentNullException`. `purity` NaN or ∉ (0,1] → `ArgumentOutOfRangeException`. A variant with TotalReads < 1, AltReads ∉ [0, N], LocalCopyNumber < 1, or Multiplicity ∉ [1, q] → `ArgumentException`. A CCF value NaN or ∉ [0,1] → `ArgumentException`. An empty variant set returns empty calls with counts 0 and ClonalFraction 0. Read counts are integer; reads are 0-based-agnostic (counts only).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `purity` and (per variant) read counts, local copy number, multiplicity.
2. For each variant, compute the per-unit-CCF allele fraction `α·M / (2(1−α)+αq)`.
3. Over a 100-point grid c ∈ [0.01, 1], compute `f(c) = min(1, perUnit·c)` and the binomial likelihood `f^a·(1−f)^(N−a)` (in log-space).
4. Normalise the grid weights to a posterior; take the posterior mean (CCF estimate) and the mass above 0.95.
5. Call clonal if that mass > 0.5, else subclonal; accumulate counts and the clonal fraction.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Clonal CCF threshold: 0.95 (strict) [1].
- Clonal posterior-probability threshold: 0.5 [1].
- CCF grid: 100 points, c ∈ [0.01, 1] [1].
- Normal diploid copy number: 2 (germline autosomal contribution 2(1−α)) [1].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| ClassifyClonality | O(n · G) | O(n) | n variants, G = 100 grid points (constant); grid weights on the stack |
| IdentifyClonalMutations | O(m) | O(k) | m CCF values, k clonal indices returned |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.ClassifyClonality(variants, purity)`: posterior-grid clonal/subclonal classification with counts and clonal fraction.
- `OncologyAnalyzer.IdentifyClonalMutations(ccfValues)`: point-estimate clonal selection (CCF > 0.95).

### 5.2 Current Behavior

The binomial likelihood is computed in log-space and the constant binomial coefficient C(N,a) is omitted (it cancels under grid normalisation). For the degenerate case of an all-zero posterior (e.g. `f ≈ 0` with `a > 0`), a flat posterior over the grid is used so the result stays well-defined (subclonal). This is a search/matching-free numerical computation, so the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `f(c) = α·M·c / (2(1−α)+α·q)` and posterior `P(c) ∝ Binomial(a|N,f(c))` on a 100-point grid c ∈ [0.01,1], uniform prior, normalised by its sum [1][2].
- Clonal iff `P(CCF > 0.95) > 0.5`; subclonal otherwise [1].
- Point-estimate clonal threshold CCF > 0.95 [1].

**Intentionally simplified:**

- Multiplicity `M` is supplied by the caller rather than inferred; **consequence:** if the caller omits it, M = 1 (heterozygous SNV) is assumed, which can overestimate CCF for multi-copy mutant loci.

**Not implemented:**

- CCF estimation itself beyond the internal posterior, and CCF clustering / subclone inference; **users should rely on:** ONCO-CCF-001 (`EstimateCCF`, `ClusterCCFValues`) once implemented.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Registry `ploidy` scalar → per-variant local copy number `q` | Assumption | API shape only; numerical rule unchanged | accepted | Landau's model uses per-locus q; mirrors prior ONCO-WGD scalar→segments decision |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty variant set | empty calls, counts 0, ClonalFraction 0 | clonal fraction undefined for empty input [contract] |
| Near-1 point estimate, shallow coverage | may be subclonal | classification uses posterior mass above 0.95, not the point estimate [1] |
| Multiplicity M = 2 (multi-copy mutant) | higher CCF for the same VAF | f(c) scales with M [2] |
| CCF exactly 0.95 (IdentifyClonalMutations) | not clonal | threshold is strict (> 0.95) [1] |

### 6.2 Limitations

Single-sample bulk model; does not jointly fit multiple samples or cluster CCFs into subclones (ONCO-CCF-001/ONCO-HETERO-001). Assumes the supplied local copy number and multiplicity are correct; mis-specified copy-number state biases CCF. The posterior mean is reported as the CCF point estimate; a credible interval is not returned.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var variants = new[]
{
    new OncologyAnalyzer.ClonalityVariant(altReads: 400, totalReads: 1000, localCopyNumber: 2), // clonal
    new OncologyAnalyzer.ClonalityVariant(altReads: 240, totalReads: 1000, localCopyNumber: 2), // subclonal
};
OncologyAnalyzer.ClonalityResult result = OncologyAnalyzer.ClassifyClonality(variants, purity: 0.8);
// result.ClonalCount == 1, result.SubclonalCount == 1, result.ClonalFraction == 0.5
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_ClassifyClonality_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_ClassifyClonality_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [ONCO-CLONAL-001-Evidence.md](../../../docs/Evidence/ONCO-CLONAL-001-Evidence.md)

## 8. References

1. Landau DA, Carter SL, Stojanov P, et al. 2013. Evolution and Impact of Subclonal Mutations in Chronic Lymphocytic Leukemia. *Cell* 152(4):714–726. https://doi.org/10.1016/j.cell.2013.01.019
2. Satas G, Zaccaria S, El-Kebir M, Raphael BJ. 2021. DeCiFering the Elusive Cancer Cell Fraction in Tumor Heterogeneity and Evolution. *Cell Systems* 12(10):1004–1018. https://doi.org/10.1016/j.cels.2021.07.006
