# Ancestry Estimation (Supervised / Projection ADMIXTURE)

| Field | Value |
|-------|-------|
| Algorithm Group | Population Genetics |
| Test Unit ID | POP-ANCESTRY-001 |
| Related Projects | Seqeron.Genomics.Population |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Estimates the *ancestry (admixture) proportions* of individuals — what fraction of each genome derives from each of K reference populations — from biallelic SNP genotypes. This implementation is the **supervised / projection** form of ADMIXTURE: the reference-population allele frequencies are supplied as fixed, known inputs, and only the per-individual ancestry vector q is estimated [1]. Estimation maximizes the binomial admixture log-likelihood (Eq. 2 of [1]) by the FRAPPE expectation-maximization (EM) update (Eq. 4 of [1]). The result is probabilistic / model-based, not exact: q converges to a local maximum of the likelihood.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Admixed individuals carry chromosomal segments inherited from several ancestral populations. Given reference panels that summarize each ancestral population's allele frequencies, an individual's genome can be modeled as a mixture and the mixing weights (ancestry fractions) inferred. ADMIXTURE adopts the likelihood model of STRUCTURE but maximizes the likelihood instead of sampling the posterior [1]. When reference individuals of known ancestry are available, the problem becomes supervised: their allele frequencies are fixed and only the unknown individuals' Q is estimated [2][3, §2.10]; equivalently, "projection" mode loads fixed allele frequencies learned from a reference panel and projects new samples onto them [3, §2.14].

### 2.2 Core Model

Let g_ij ∈ {0,1,2} be the number of copies of allele 1 at SNP j of individual i; q_ik the fraction of individual i's genome from population k; f_kj the allele-1 frequency in population k at SNP j [1]. Under random union of gametes the genotype is binomial with success probability Σ_k q_ik f_kj, giving the log-likelihood (Eq. 2 of [1]):

```
L(Q,F) = Σ_i Σ_j { g_ij · ln( Σ_k q_ik f_kj ) + (2 − g_ij) · ln( Σ_k q_ik (1 − f_kj) ) }
```

With F held fixed, the FRAPPE EM update for the ancestry vector is (Eq. 4 of [1]):

```
q_ik^{n+1} = (1 / 2J) · Σ_j [ g_ij · a^n_ijk + (2 − g_ij) · b^n_ijk ]
a^n_ijk = q^n_ik f_kj / ( Σ_m q^n_im f_mj )
b^n_ijk = q^n_ik (1 − f_kj) / ( Σ_m q^n_im (1 − f_mj) )
```

where J is the number of (informative) SNPs. Convergence is declared once the log-likelihood gain falls below ε (Eq. 5 of [1]; ADMIXTURE default ε = 10⁻⁴).

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Hardy-Weinberg proportions within each ancestral source (random union of gametes) | The binomial form of Eq. 2 no longer holds; estimates biased |
| ASM-02 | Reference allele frequencies f_kj are correct and fixed (supervised/projection) | Wrong panels propagate directly into wrong q |
| ASM-03 | SNPs approximately independent (no strong LD) | EM treats markers as independent factors; LD inflates effective information |
| ASM-04 | Reference panels share the same SNP set/order as the genotype vector (index alignment) | Mismatched indexing scrambles f_kj |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Σ_k q_ik = 1 for every individual | Eq. 4 divides by 2J and Σ_k(g·a + (2−g)·b) = g + (2−g) = 2 per SNP; constraint Σ_k q_ik = 1 [1] |
| INV-02 | 0 ≤ q_ik ≤ 1 | q_ik ≥ 0 in Eq. 4 and they sum to 1 (constraint q_ik ≥ 0 [1]) |
| INV-03 | L(Q^{n+1}) ≥ L(Q^n) | EM is a monotone ascent algorithm; basis of the Eq. 5 stopping rule [1] |
| INV-04 | Identical reference panels keep a uniform q uniform | Eq. 4 numerators/denominators are proportional ⇒ uniform is a fixed point |

### 2.5 Comparison with Related Methods

| Aspect | This (supervised ADMIXTURE EM) | Unsupervised ADMIXTURE |
|--------|-------------------------------|------------------------|
| Allele frequencies F | fixed, given by reference panels | jointly estimated with Q |
| Optimizer | FRAPPE EM (Eq. 4) | block relaxation + quasi-Newton acceleration [1] |
| Label identifiability | pinned by reference labels | up to K! permutations [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| individuals | IEnumerable<(string IndividualId, IReadOnlyList<int> Genotypes)> | required | Individuals; `Genotypes[j]` = allele-1 count at SNP j | each value in {0,1,2}; length = panel SNP count |
| referencePops | IEnumerable<(string PopulationId, IReadOnlyList<double> AlleleFrequencies)> | required | Fixed reference allele-1 frequencies, one per SNP | each in [0,1]; same SNP order across panels |
| maxIterations | int | 100 | EM iteration budget per individual | ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| AncestryProportion.IndividualId | string | Echoes the input id |
| AncestryProportion.Proportions | IReadOnlyDictionary<string,double> | Ancestry fraction keyed by reference-population id; sums to 1 (INV-01) |

### 3.3 Preconditions and Validation

Genotypes are 0-based, indexed parallel to the reference allele-frequency vectors. Empty `individuals` or empty `referencePops` → empty result. An individual whose genotype length ≠ the panel SNP count is skipped. A genotype value outside {0,1,2} is treated as missing: that SNP contributes no likelihood term (Eq. 2) and is excluded from J for that individual. No exceptions are thrown for these documented input classes.

## 4. Algorithm

### 4.1 High-Level Steps

1. Materialize individuals and reference panels; if either is empty, yield nothing. Let K = number of panels, J = panel SNP count.
2. For each individual whose genotype length = J: initialize q uniformly (1/K).
3. Repeat up to `maxIterations`: for each informative SNP compute the admixed frequencies Σ_m q_m f_mj and Σ_m q_m(1−f_mj), accumulate the Eq. 4 terms g·a + (2−g)·b per population, then set q = accumulator / 2J.
4. After each update recompute the Eq. 2 log-likelihood; stop early when the gain < 10⁻⁴ (Eq. 5).
5. Emit q keyed by reference-population id.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- ε = 10⁻⁴ log-likelihood convergence tolerance (`AncestryLogLikelihoodTolerance`), from Eq. 5 of [1].
- A genotype is informative iff it lies in {0,1,2}; J is the count of informative SNPs for the individual.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Estimate one individual | O(iterations · J · K) | O(K) | per-iteration term is O(JK) (one pass over SNPs × populations); [1] gives O(IJK²) for joint Q+F, here F fixed ⇒ no K² panel-update term |
| Full call | O(I · iterations · J · K) | O(I·K) output | I individuals |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PopulationGeneticsAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs)

- `PopulationGeneticsAnalyzer.EstimateAncestry(individuals, referencePops, maxIterations)`: public entry; iterates individuals and emits `AncestryProportion`.
- `EstimateIndividualAncestry(...)` (private): FRAPPE EM (Eq. 4) for one individual with F fixed.
- `AncestryLogLikelihood(...)` (private): Eq. 2 log-likelihood used for the Eq. 5 stopping rule.

### 5.2 Current Behavior

F is never updated — this is strictly the supervised/projection case. The EM divides accumulated responsibilities by 2J (J = informative SNP count), which keeps Σ_k q_ik = 1 exactly without a separate renormalization step. Early stopping uses the Eq. 2 log-likelihood gain. No search/matching over a text is performed, so the repository suffix tree is **not applicable** to this unit (the work is per-locus arithmetic, not substring occurrence finding).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Binomial admixture log-likelihood, Eq. 2 of [1].
- FRAPPE EM ancestry update with a_ijk and b_ijk, Eq. 4 of [1].
- Σ_k q_ik = 1, q_ik ≥ 0, 0 ≤ f_kj ≤ 1 constraints [1].
- Eq. 5 convergence rule with ADMIXTURE default ε = 10⁻⁴ [1].
- Supervised/projection semantics: F fixed, only Q estimated [2][3].

**Intentionally simplified:**

- Optimizer: plain FRAPPE EM only; **consequence:** convergence is slower than ADMIXTURE's accelerated block relaxation [1], but the maximized objective and the q estimate are the same.

**Not implemented:**

- Joint estimation of allele frequencies F (unsupervised ADMIXTURE); **users should rely on:** an external ADMIXTURE/STRUCTURE run to learn F, then pass it here as fixed reference panels.
- Standard-error estimation via block bootstrap [1]; **users should rely on:** the reference ADMIXTURE tool.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | EM budget via `maxIterations` plus Eq. 5 early stop | Assumption | Bounds work; same fixed point as pure Eq. 5 stopping | accepted | ASM/Evidence Assumption 1 |
| 2 | Out-of-range genotype treated as missing; length-mismatch individual skipped | Assumption | Determines which SNPs/individuals contribute | accepted | Evidence Assumption 2 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty individuals or empty reference panels | Empty result | Nothing to estimate |
| Genotype length ≠ panel SNP count | Individual skipped | Index alignment cannot be guaranteed (ASM-04) |
| Genotype value outside {0,1,2} | SNP skipped (missing) | No binomial term in Eq. 2 |
| All SNPs missing for an individual | Returns the uniform prior | No informative term to update q |
| Identical reference panels | q stays uniform | INV-04 |

### 6.2 Limitations

Supervised/projection only — does not estimate F. Treats SNPs as independent (ASM-03), so strong LD overstates information. Quality depends entirely on the supplied reference panels (ASM-02). EM finds a local maximum; with informative panels the surface is effectively unimodal for a single individual, but pathological panels could admit ties. No standard errors are produced.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var individuals = new[] { ("ind1", (IReadOnlyList<int>)new[] { 2, 0 }) };
var refs = new[]
{
    ("A", (IReadOnlyList<double>)new[] { 0.8, 0.2 }),
    ("B", (IReadOnlyList<double>)new[] { 0.2, 0.8 }),
};
var result = PopulationGeneticsAnalyzer.EstimateAncestry(individuals, refs, maxIterations: 1).Single();
// result.Proportions["A"] == 0.8, result.Proportions["B"] == 0.2  (one EM iteration, Eq. 4)
```

**Numerical walk-through:** start q = (0.5,0.5), f_A=(0.8,0.2), f_B=(0.2,0.8), g=(2,0).
SNP1 (g=2): mix1 = 0.5·0.8+0.5·0.2 = 0.5; contribution to A = 2·(0.5·0.8/0.5) = 1.6, to B = 2·(0.5·0.2/0.5) = 0.4.
SNP2 (g=0): mix2 = 0.5·0.8+0.5·0.2 = 0.5; contribution to A = 2·(0.5·0.8/0.5) = 1.6, to B = 0.4.
q_A = (1.6+1.6)/(2·2) = 0.8, q_B = (0.4+0.4)/4 = 0.2.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PopulationGeneticsAnalyzer_EstimateAncestry_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Population/PopulationGeneticsAnalyzer_EstimateAncestry_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [POP-ANCESTRY-001-Evidence.md](../../../docs/Evidence/POP-ANCESTRY-001-Evidence.md)
- Related algorithms: [F_Statistics](../Population_Genetics/F_Statistics.md)

## 8. References

1. Alexander DH, Novembre J, Lange K. 2009. Fast model-based estimation of ancestry in unrelated individuals. Genome Research 19(9):1655–1664. https://doi.org/10.1101/gr.094052.109
2. Alexander DH, Lange K. 2011. Enhancements to the ADMIXTURE algorithm for individual ancestry estimation. BMC Bioinformatics 12:246. https://doi.org/10.1186/1471-2105-12-246
3. Alexander DH, Shringarpure SS, Novembre J, Lange K. ADMIXTURE 1.4 Software Manual (§2.10 Supervised analysis, §2.14 Projection analysis). https://dalexander.github.io/admixture/admixture-manual.pdf
