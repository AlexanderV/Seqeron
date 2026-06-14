# Mutational Signature Fitting (NNLS Refitting and Cosine Similarity)

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology / Mutational Signatures |
| Test Unit ID | ONCO-SIG-002 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Given an observed mutational catalog (e.g. a 96-channel SBS spectrum from ONCO-SIG-001) and a set of
caller-supplied reference signatures, *signature refitting* estimates how much each reference signature
contributes to the catalog. The contribution (exposure) vector is the non-negative least-squares (NNLS)
solution that best reconstructs the catalog as a non-negative linear combination of the signatures [1][3].
The module also computes the cosine similarity between two vectors — used both to compare catalogs/signatures
and to score how faithfully the fit reconstructs the observed catalog [1][4]. The fit is exact (a
deterministic convex optimisation), not heuristic. Reference signature profiles are **not** hardcoded; the
caller supplies them (COSMIC SBS profiles cannot be reproduced authoritatively here), so the unit is a
*Framework* that performs the fit on user-supplied data.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A mutational catalog is a non-negative vector d of mutation-type counts (96 channels for SBS). A reference
signature is a non-negative vector over the same channels describing one mutational process. Refitting
"considers that the signatures are known and the goal is to estimate their contributions given a mutational
catalogue," i.e. projecting the catalog onto the non-negative cone spanned by the signatures [1].

### 2.2 Core Model

**Cosine similarity** between two vectors A and B over n components [1][4]:

```
sim(A, B) = Σᵢ Aᵢ Bᵢ / ( √(Σᵢ Aᵢ²) · √(Σᵢ Bᵢ²) )
```

the dot product divided by the product of the Euclidean norms; the result lies in [0, 1] for non-negative
inputs (0 = independent/orthogonal, 1 = identical direction) [1].

**Signature fitting (NNLS)** [1][3]:

```
minₓ ‖ S · x − d ‖₂²   subject to   x ≥ 0
```

where S is the signature matrix (column j = signature j), d is the observed catalog, and x is the fitted
exposure vector. The reconstructed catalog is S·x (R = S·W in deconstructSigs notation [2]); coefficients
must be non-negative because "negative contributions make no biological sense" [2]. Exposures are normalised
into proportions ("the weights W are normalized between 0 and 1" [2]). The cosine similarity between d and
S·x is the reconstruction-quality measure; ≥ 0.95 indicates a successful reconstruction [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | sim(A,B) ∈ [0,1] for non-negative A,B | dot product ≤ product of norms (Cauchy-Schwarz); both non-negative [1] |
| INV-02 | sim(A,A) = 1 for any non-zero A | numerator = ‖A‖², denominator = ‖A‖² [1] |
| INV-03 | sim(A, k·B) = sim(A,B) for k>0 | cosine of an angle is scale-invariant [4] |
| INV-04 | every fitted exposure ≥ 0 | NNLS constraint x ≥ 0 [3] |
| INV-05 | ‖S·x − d‖² ≤ ‖d‖² | x = 0 is feasible, so the minimiser is no worse [1][3] |
| INV-06 | normalised exposures sum to 1 (Σ>0) else all 0 | division by Σ exposures [2] |
| INV-07 | NNLS = unconstrained LS when the latter is non-negative | no constraint is active [3] |

### 2.5 Comparison with Related Methods

| Aspect | NNLS active-set (this) | deconstructSigs greedy heuristic |
|--------|------------------------|----------------------------------|
| Objective | exact convex minimiser of ‖Sx−d‖², x≥0 [1][3] | iterative SSE reduction with a contribution threshold [2] |
| Determinism | deterministic, reproducible | heuristic forward-selection [2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `a`, `b` | `IReadOnlyList<double>` | required | vectors for cosine similarity | same length, non-empty |
| `catalog` | `IReadOnlyList<double>` | required | observed catalog d | length = signature channel count; non-negative |
| `signatures` | `IReadOnlyList<IReadOnlyList<double>>` | required | reference signatures S (one vector per signature) | ≥1 signature, equal-length, non-empty |
| `exposures` | `IReadOnlyList<double>` | required | per-signature weights for reconstruction | count = signature count |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Exposures` | `IReadOnlyList<double>` | NNLS solution x (per-signature contribution, ≥ 0) |
| `NormalizedExposures` | `IReadOnlyList<double>` | exposures / Σexposures (sum 1 when Σ>0, else all 0) |
| `Reconstruction` | `IReadOnlyList<double>` | S·x |
| `ReconstructionCosineSimilarity` | `double` | cosine similarity of d and S·x |

### 3.3 Preconditions and Validation

Null inputs throw `ArgumentNullException`. Empty / ragged / count-mismatched / dimension-mismatched inputs
throw `ArgumentException`. Vectors are indexed 0-based; channel order is the caller's responsibility (it must
be consistent between `catalog` and each signature). A zero-norm vector in `CosineSimilarity` yields 0.0
(division by zero is undefined; treated as no shared direction).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs and the common channel count.
2. Solve `minₓ ‖S·x − d‖²`, `x ≥ 0` with the Lawson-Hanson active-set algorithm [3].
3. Reconstruct S·x and normalise exposures into proportions [2].
4. Compute the cosine similarity of d and S·x as the reconstruction-quality score [1].

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Lawson-Hanson active set [3]: maintain passive set P (free variables) and active set R (clamped to 0).
Initialise x = 0; while R ≠ ∅ and max over R of the gradient w = Sᵀ(d − Sx) exceeds ε, move the max-gradient
index into P, solve the unconstrained LS on P via the normal equations
`s_P = ((S_P)ᵀ S_P)⁻¹ (S_P)ᵀ d`; while any passive component ≤ 0, take the bounded step
`α = min x_i/(x_i − s_i)` (over i∈P with s_i ≤ 0), update x = x + α(s − x), move ≤ 0 indices back to R, and
re-solve; then set x = s. ε = 1e-12; the normal equations are solved by Gaussian elimination with partial
pivoting.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CosineSimilarity | O(n) | O(1) | n = vector length |
| ReconstructCatalog | O(k·n) | O(n) | k signatures, n channels |
| FitSignatures (NNLS) | O(k³ + k²·n) per outer iteration; ≤ O(k) outer iterations | O(k² + n) | k = #signatures, n = #channels; small k in practice |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.CosineSimilarity(a, b)`: cosine similarity of two vectors.
- `OncologyAnalyzer.FitSignatures(catalog, signatures)`: NNLS refit → `SignatureFitResult`.
- `OncologyAnalyzer.ReconstructCatalog(signatures, exposures)`: S·x.

### 5.2 Current Behavior

The NNLS solver uses dense normal equations (small k); a singular passive-set matrix (collinear signatures)
leaves the affected component at 0 rather than throwing. Reference signatures are supplied per call; nothing
is cached or hardcoded. No substring/pattern search is involved, so the repository suffix tree is **not
applicable** to this unit.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Cosine similarity sim = ΣAB / (√ΣA² · √ΣB²) [1][4].
- NNLS objective minₓ ‖S·x − d‖², x ≥ 0, solved by Lawson-Hanson active set [1][3].
- Reconstruction S·x and proportion normalisation of exposures [2].
- Reconstruction-quality cosine between d and S·x [1].

**Intentionally simplified:**

- deconstructSigs' greedy forward-selection heuristic and its contribution threshold are not reproduced;
  the deterministic NNLS minimiser is used instead. **consequence:** exposures are the exact convex optimum
  rather than the heuristic's thresholded subset, so very small contributions are not pruned to zero.

**Not implemented:**

- de novo NMF signature *extraction*; bootstrap confidence intervals (ONCO-SIG-003) — **users should rely on:**
  caller-supplied reference signatures for this unit; later ONCO-SIG units for exposure CIs.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Zero catalog d = 0 | exposures all 0; reconstruction 0; proportions all 0 | only feasible minimiser of ‖Sx‖², x≥0 [1][3] |
| Zero-norm vector in cosine | returns 0.0 | cosine undefined (÷0); no shared direction (Assumption) |
| Identical vectors | cosine = 1 | INV-02 [1] |
| Orthogonal vectors | cosine = 0 | dot product 0 [1] |
| Unconstrained LS coefficient < 0 | clamped to 0, refit on remaining set | active-set constraint [3] |
| null / empty / dimension mismatch | `ArgumentNullException` / `ArgumentException` | input validation |

### 6.2 Limitations

Reference signatures must be supplied by the caller (not bundled). The NNLS solver targets small signature
counts (dense normal equations); it is not tuned for thousands of signatures. Collinear signatures make the
decomposition non-unique; the solver returns one minimiser. No statistical uncertainty (confidence intervals)
is produced here.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** S = [[1,1],[0,1]] (sig1 = [1,0] over channels, sig2 = [1,1]), d = [0,1].
Unconstrained normal equations give x = [−1, 1]; x₁ < 0 ⇒ clamp sig1 to 0 and refit sig2 alone:
x₂ = ([1,1]·[0,1]) / ([1,1]·[1,1]) = 1/2. NNLS solution = [0, 0.5]; reconstruction = [0.5, 0.5].

**API usage example:**

```csharp
var catalog = new double[] { 3, 5 };
var signatures = new IReadOnlyList<double>[] { new double[] {1,0}, new double[] {0,1} };
var fit = OncologyAnalyzer.FitSignatures(catalog, signatures);
// fit.Exposures == [3, 5]; fit.ReconstructionCosineSimilarity == 1.0
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_FitSignatures_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_FitSignatures_Tests.cs) — covers INV-01..INV-07
- Evidence: [ONCO-SIG-002-Evidence.md](../../../docs/Evidence/ONCO-SIG-002-Evidence.md)
- Related algorithms: [SBS96 Trinucleotide Context Catalog](./SBS96_Trinucleotide_Context_Catalog.md)

## 8. References

1. Blokzijl F, Janssen R, van Boxtel R, Cuppen E. 2018. MutationalPatterns: comprehensive genome-wide analysis of mutational processes. Genome Medicine 10:33. https://pmc.ncbi.nlm.nih.gov/articles/PMC5922316/
2. Rosenthal R, McGranahan N, Herrero J, Taylor BS, Swanton C. 2016. deconstructSigs. Genome Biology 17:31. https://pmc.ncbi.nlm.nih.gov/articles/PMC4762164/
3. Lawson CL, Hanson RJ. 1974. Solving Least Squares Problems, Ch. 23 (active-set NNLS). Prentice-Hall. https://en.wikipedia.org/wiki/Non-negative_least_squares
4. Pan W, Wang X. 2020. iMutSig: a web application to identify the most similar mutational signature. https://pmc.ncbi.nlm.nih.gov/articles/PMC7702159/
