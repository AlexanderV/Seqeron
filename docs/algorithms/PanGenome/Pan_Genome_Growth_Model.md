# Pan-Genome Growth Model (Heaps' Law)

| Field | Value |
|-------|-------|
| Algorithm Group | PanGenome |
| Test Unit ID | PANGEN-HEAP-001 |
| Related Projects | Seqeron.Genomics.Metagenomics |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

The pan-genome growth model quantifies how many *new* gene clusters are discovered as additional genomes of a species are sequenced, and from that decides whether the pan-genome is **open** (new genes keep accumulating) or **closed** (the gene pool is bounded). It fits Heaps' law `n(N) = K · N^(−α)` to the number of new gene clusters contributed by the N-th genome, averaged over many random genome orderings, following the micropan `heaps()` reference implementation [3] of the model introduced by Tettelin et al. [1][2]. The fit is a deterministic least-squares estimate over box-constrained parameters; the open/closed verdict is the threshold rule α < 1 ⇒ open [2][3].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A species pan-genome is the union of all gene families (clusters / ortholog groups) found across its sequenced strains, partitioned into a **core genome** (present in essentially all strains) and a **dispensable genome** (subset-shared and strain-specific genes) [2]. In *S. agalactiae* the core is ≈80% of any single genome, and sequencing successive strains keeps revealing strain-specific genes (161 new at the 2nd genome, 54 at the 5th, with an asymptote of ≈33 new genes per additional genome) — evidence of an open pan-genome [1].

### 2.2 Core Model

For a binary gene presence/absence matrix (genomes × clusters), order the genomes and record, for the N-th genome (N = 2, 3, …, G), the number of clusters that appear for the first time at position N. Heaps' law models this new-gene curve as a power law [2][3]:

```
n(N) = K · N^(−α)
```

where `K` is the intercept and `α` the decay exponent [3]. Parameters are estimated by minimizing the micropan objective over pooled (N, count) points from many random orderings [3]:

```
J(K, α) = sqrt( Σ_i (y_i − K · x_i^(−α))² ) / |x|
```

with box constraints `K ∈ [0, 10000]`, `α ∈ [0, 2]` and start values `(mean count at N=2, 1)` [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Open ⇔ α < 1; closed ⇔ α > 1 | micropan rule "If alpha>1.0 ... closed, if alpha<1.0 ... open" [3]; Tettelin 2008 [2] |
| INV-02 | New-gene count at position N (N≥2) = clusters present at N but absent in all earlier positions (first appearance) | micropan `(cm==1)[i] & (cm==0)[i−1]` over the cumulative-sum matrix [3] |
| INV-03 | Presence is binary: a cluster counts once per genome regardless of copy number | micropan binarization `pan.matrix[pan.matrix>0] <- 1` [3] |
| INV-04 | Fitted α ∈ [0, 2] and K ∈ [0, 10000] | optim box bounds [3] |
| INV-05 | For points exactly on a single power curve within bounds, the recovered (K, α) equal the analytic solution | unique global minimum of J at J=0 [3] |
| INV-06 | predictor(N) = K · N^(−α) is non-increasing in N for α ≥ 0 | monotonicity of `N^(−α)` [3] |

### 2.5 Comparison with Related Methods

| Aspect | Heaps' law (this model) | Tettelin 2005 exponential decay |
|--------|-------------------------|---------------------------------|
| New-gene form | `K · N^(−α)` (power law) [2][3] | `κ·exp(−n/τ) + tg(θ)` (exp. + asymptote) [1] |
| Openness criterion | α < 1 ⇒ open [2][3] | asymptote tg(θ) > 0 ⇒ open [1] |
| Fit method | least squares over (K, α), box-bounded [3] | nonlinear least squares (3 params) [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `matrix` | `IEnumerable<GenePresenceRow>` | required | Gene presence/absence rows, one per genome | nullable; rows with no clusters → degenerate fit |
| `permutations` | `int` | 100 | Random genome orderings pooled for the fit (micropan n.perm) | ≥ 1 (clamped) |
| `genomes` (overload) | `IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>` | required | Genomes to cluster then fit | nullable/empty → degenerate fit |
| `identityThreshold` (overload) | `double` | 0.9 | CD-HIT identity cutoff for clustering | [0,1] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Intercept` | `double` | Heaps' law K (micropan "Intercept"), in [0, 10000] |
| `Alpha` | `double` | Heaps' law decay exponent α (micropan "alpha"), in [0, 2] |
| `IsOpen` | `bool` | `true` ⇔ α < 1 (open pan-genome) |
| `PredictNewGenes` | `Func<int,double>` | `N ↦ Intercept · N^(−Alpha)` (expected new genes at the N-th genome) |

### 3.3 Preconditions and Validation

Null or empty input, or fewer than 2 genomes, yields a degenerate fit `(Intercept=0, Alpha=0, IsOpen=false, predictor→0)` rather than an exception, because the new-gene curve `N = 2..G` is undefined below 2 genomes [3]. Cluster presence is read from `GenePresenceRow.GenePresence` (true = present); duplicate/over-counted presence is collapsed to binary [3]. Cluster ordering across rows is stabilized by first appearance so the matrix columns are deterministic.

## 4. Algorithm

### 4.1 High-Level Steps

1. Collect the distinct cluster IDs (matrix columns) and build a binary presence vector per genome (any present → 1) [3].
2. For each of `permutations` orderings (the first is the natural input order; the rest Fisher-Yates shuffles), walk genomes in order and count, at each position N ≥ 2, clusters appearing for the first time [3].
3. Pool all (N, new-count) points: `x = rep(2:G, permutations)`, `y` = stacked counts [3].
4. Minimize `J(K,α) = sqrt(Σ(y − K·x^(−α))²)/|x|` over `K ∈ [0,10000]`, `α ∈ [0,2]` from start `(mean y at N=2, 1)` [3].
5. Return K, α, `IsOpen = α < 1`, and the predictor `N ↦ K·N^(−α)` [3].

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- New-gene rule (INV-02): cluster c is new at position N iff `cumsum_c[N] == 1 && cumsum_c[N−1] == 0` [3].
- Open/closed rule (INV-01): α < 1 ⇒ open, α > 1 ⇒ closed; the boundary α = 1 is treated as not-open (strict inequality) [3].
- Box constraints: K ∈ [0, 10000], α ∈ [0, 2]; start `(mean count at N=2, 1)` [3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| New-gene curve (all perms) | O(P · G · C) | O(G · C) | P permutations, G genomes, C clusters |
| Objective evaluation | O(P · G) | O(P · G) | per (K,α) point pool size |
| Bounded minimization | O(I · P · G) | O(P · G) | I = coordinate-descent iterations (deterministic) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [PanGenomeAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs)

- `PanGenomeAnalyzer.FitHeapsLaw(IEnumerable<GenePresenceRow>, int)`: canonical Heaps' fit over a presence/absence matrix.
- `PanGenomeAnalyzer.FitHeapsLaw(IReadOnlyDictionary<...>, double, int)`: convenience overload — clusters genomes, builds the matrix, then delegates to the canonical fit.
- `PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(genomes, clusters)`: builds the binary matrix the fit consumes.

### 5.2 Current Behavior

The new-gene curve is computed exactly per INV-02 over the binary matrix. The micropan `optim(L-BFGS-B)` call is realized as deterministic coordinate descent with geometric step refinement over the identical objective, bounds, and start point; for data on an exact power curve within bounds it converges to the analytic optimum to < 1e-9 (see TestSpec M1/M2). Permutations use a fixed seed and the natural input order for the first ordering, so fixed-order / permutation-invariant fits are exactly reproducible. **Search reuse:** the suffix tree is **not** applicable — this unit does no substring/pattern search; it counts set first-appearances and fits a power curve, so no occurrence-enumeration structure is involved.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Model `n(N) = K · N^(−α)` and objective `J = sqrt(Σ(y−K·x^(−α))²)/|x|` [3].
- First-appearance new-gene counting `(cm==1)[i] & (cm==0)[i−1]`, index starting at N=2 [3].
- Binarization of presence (>0 → 1) [3].
- Box bounds K∈[0,10000], α∈[0,2]; start `(mean y at N=2, 1)` [3].
- Open/closed rule α < 1 ⇒ open [2][3].

**Intentionally simplified:**

- Optimizer: micropan uses R `optim(L-BFGS-B)`; this port uses deterministic coordinate descent over the identical objective/bounds/start; **consequence:** none for in-bounds power-curve data (analytic optimum recovered < 1e-9); for noisy real data the minimizer may differ by sub-tolerance amounts from R's gradient solver.

**Not implemented:**

- Tettelin 2005 exponential-decay-with-asymptote core/new-gene model (`κ·exp(−n/τ)+Ω`); **users should rely on:** this Heaps' law (power-law) estimate, which is the openness criterion of Tettelin 2008.
- Confidence intervals / percentile bands over permutations; **users should rely on:** repeated calls with larger `permutations` if dispersion is needed.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Coordinate descent vs L-BFGS-B | Assumption | Optimizer method differs; objective/bounds/start identical | accepted | Verified to match analytic optimum < 1e-9 (INV-05) |
| 2 | Fixed-seed permutations, natural first order | Assumption | RNG stream differs from R `sample()` | accepted | Makes fixed-order/invariant fits reproducible; averaging principle unchanged |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| 0 or 1 genome | Degenerate fit (0,0,false), predictor→0 | curve N=2..G empty [3] |
| null / empty matrix | Degenerate fit, no exception | nothing to fit |
| Duplicate presence of a cluster | Counts once (binary) | binarization [3] |
| Constant new-gene curve | α = 0, K = mean count, open | best power fit of a constant [3] |
| Curve on exact power law | (K, α) recovered exactly | unique J=0 minimum (INV-05) |

### 6.2 Limitations

The fit is meaningful only with several genomes (Tettelin 2008 emphasizes many strains [2]); with few genomes α is poorly determined. The α∈[0,2] bound caps very steep decays at 2 (matching micropan). Permutation averaging here uses a fixed seed; absolute α on noisy real matrices is permutation-sensitive and should be read with the openness verdict, not as a high-precision constant.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var genomes = /* IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> */;
var clusters = PanGenomeAnalyzer.ClusterGenes(genomes).ToList();
var matrix = PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(genomes, clusters);
var fit = PanGenomeAnalyzer.FitHeapsLaw(matrix, permutations: 100);
// fit.IsOpen == fit.Alpha < 1.0; fit.PredictNewGenes(50) = expected new genes at genome 50
```

**Numerical walk-through:** A fixed-order matrix where genome 2 introduces 8 clusters absent in genome 1 and genome 3 introduces 4 clusters absent in genomes 1–2 gives the new-gene curve `x=[2,3], y=[8,4]`. Solving `K·2^(−α)=8`, `K·3^(−α)=4`: `8/4=(3/2)^α ⇒ α = ln2/ln(3/2) = 1.7095113`, `K = 8·2^α = 26.1640014`. Since α > 1, the pan-genome is **closed**.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [PanGenomeAnalyzer_FitHeapsLaw_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/PanGenomeAnalyzer_FitHeapsLaw_Tests.cs) — covers `INV-01`–`INV-06`
- Evidence: [PANGEN-HEAP-001-Evidence.md](../../../docs/Evidence/PANGEN-HEAP-001-Evidence.md)

## 8. References

1. Tettelin H, Masignani V, Cieslewicz MJ, Donati C, Medini D, et al. 2005. Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: Implications for the microbial "pan-genome". PNAS 102(39):13950–13955. https://doi.org/10.1073/pnas.0506758102
2. Tettelin H, Riley D, Cattuto C, Medini D. 2008. Comparative genomics: the bacterial pan-genome. Current Opinion in Microbiology 11(5):472–477. https://doi.org/10.1016/j.mib.2008.09.006
3. Snipen L, Liland KH. micropan: Microbial Pan-Genome Analysis — `heaps()` (R/powerlaw.R). https://raw.githubusercontent.com/larssnip/micropan/master/R/powerlaw.R
