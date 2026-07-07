# De-novo Mutational-Signature Extraction via NMF

| Field | Value |
|-------|-------|
| Algorithm Group | Oncology |
| Test Unit ID | ONCO-SIG-002 |
| Related Projects | Seqeron.Genomics.Oncology |
| Implementation Status | Production |
| Last Reviewed | 2026-06-23 |

## 1. Overview

De-novo mutational-signature extraction discovers the latent mutational signatures operative across a cohort
directly from a mutation-count matrix, without any reference catalog. Given a non-negative matrix V of mutation
counts (channels × samples — e.g. the 96 SBS trinucleotide channels) and a number of signatures k, it factorises
V ≈ W·H with W ≥ 0 (the signatures, channels × k) and H ≥ 0 (their per-sample exposures, k × samples) using
Non-negative Matrix Factorization (NMF) [1][3]. The caller may pick either the squared-Euclidean (Frobenius)
objective or the Kullback-Leibler / Poisson objective (the one SigProfiler uses) [1][5]. The rank k may be
supplied directly or chosen automatically by the SigProfiler/Brunet model-stability procedure — many NMF
restarts summarised by the cophenetic correlation of the consensus matrix and per-signature silhouette stability
plus reconstruction error [6][5]. Extracted signatures can be labelled against a caller-supplied reference set
(e.g. COSMIC) by cosine similarity [5]. This is the complementary operation to signature *refitting*
(`FitSignatures`), which takes reference signatures as input and solves only for exposures. NMF is a non-convex
heuristic that converges to a local optimum dependent on initialisation; every entry point is deterministic for a
fixed seed.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Somatic mutations in a tumour are the cumulative imprint of the mutational processes active during its
evolution. Alexandrov et al. (2013) framed the recovery of these processes as a *blind source separation*
problem: each sample's mutation spectrum is a non-negative mixture of a small number of latent signatures, and
NMF separates the mixture into the signatures (sources) and their per-sample contributions (mixing weights) [3].
Single-base substitutions are classified into 96 channels (6 substitution types × 4 5′ bases × 4 3′ bases) [4].

### 2.2 Core Model

NMF approximates a non-negative matrix V (m channels × n samples) as a product of two non-negative factors of
inner dimension k [1]:

```
V ≈ W · H,   W ≥ 0 (m × k),   H ≥ 0 (k × n)
```

The columns of W are the k signatures; the rows of H are their exposures across the n samples. The factors
minimise the squared Euclidean (Frobenius) objective [1][2]:

```
F(W, H) = ‖V − W·H‖²_F = Σ_{i,μ} (V_{iμ} − (W·H)_{iμ})²
```

Lee & Seung (2001) give multiplicative update rules that leave the objective non-increasing (their Theorem 1) [1]:

```
H_aμ ← H_aμ · (Wᵀ V)_aμ / (Wᵀ W H)_aμ
W_ia ← W_ia · (V Hᵀ)_ia / (W H Hᵀ)_ia
```

Each update is a ratio applied element-wise to the current factor, which preserves non-negativity and, when
V = W·H exactly, becomes a matrix of ones (a fixed point) [2]. The signatures (columns of W) are L1-normalised
so each sums to 1 — a probability distribution over the channels, per the COSMIC / SigProfiler convention — with
the removed scale absorbed into H so that W·H is unchanged [3][4].

**Kullback-Leibler / Poisson objective.** Lee & Seung also give the generalized KL divergence objective and its
update rules (their Theorem 2) [1]; this is the objective SigProfiler optimises for mutational signatures [5]:

```
D(V‖WH) = Σ_{i,μ} ( V_{iμ} log(V_{iμ}/(WH)_{iμ}) − V_{iμ} + (WH)_{iμ} )
H_aμ ← H_aμ · ( Σ_i W_ia V_iμ/(WH)_iμ ) / ( Σ_i W_ia )
W_ia ← W_ia · ( Σ_μ H_aμ V_iμ/(WH)_iμ ) / ( Σ_μ H_aμ )
```

D(V‖WH) is monotonically non-increasing under these updates (Theorem 2) [1]. The caller selects Frobenius vs KL.

**Automatic rank selection (Brunet 2004; Alexandrov 2013).** For each candidate rank k, NMF is run many times
from a fixed, deterministic seed sequence. Each run assigns every sample to the cluster of its largest exposure
(metagene), giving a connectivity matrix C (C_ij = 1 iff samples i, j share the argmax metagene); the **consensus
matrix** is the mean connectivity over runs, and the **cophenetic correlation coefficient** is the Pearson
correlation between the consensus-induced sample distances (1 − consensus) and their cophenetic distances from
average-linkage hierarchical clustering [6]. Cross-run signature **stability** is the per-signature average
silhouette width (cosine distance) of the run signatures clustered to a consensus signature [5][7]. The rank is
chosen as the largest candidate whose average stability ≥ 0.80 with minimum per-signature stability ≥ 0.20
(SigProfiler defaults), trading reproducibility against reconstruction error [5].

**Reference matching.** Each extracted signature is labelled with the caller-supplied reference signature of
maximal cosine similarity (the per-signature reduction of SigProfiler's Hungarian "maximise total cosine"
pairing) [5].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-NMF-01 | The catalog is a non-negative linear mixture of k signatures | If the true rank ≠ k, the factorisation under/over-fits; signatures blend or split |
| ASM-NMF-02 | Noise model = the chosen objective (Frobenius squared-error or KL/Poisson) | The KL objective matches count data with low totals; Frobenius assumes additive Gaussian-like error |
| ASM-NMF-03 | Non-negative random initialisation is adequate to reach a useful local optimum | NMF is non-convex; a poor initialisation can stall at a worse local minimum |
| ASM-NMF-04 | Consensus clustering = argmax-metagene assignment + average-linkage cophenetic distance | The Brunet (2004) cophenetic value depends on this clustering; a different linkage changes the stability curve |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-NMF-01 | W ≥ 0 and H ≥ 0 throughout | Multiplicative updates scale non-negative factors by non-negative ratios [1][2] |
| INV-NMF-02 | Each signature (column of W) sums to 1 after normalisation | Explicit L1 column-normalisation; scale absorbed into H [3][4] |
| INV-NMF-03 | The selected objective (‖V − WH‖²_F or D(V‖WH)) is non-increasing across iterations | Lee & Seung (2001) Theorems 1 & 2 [1] |
| INV-NMF-04 | Exactly factorable V = W₀·H₀ ⇒ residual → 0 at convergence | All-ones multiplicative factor at the V = WH fixed point [2] |
| INV-NMF-05 | Determinism for a fixed seed (incl. the derived per-run seed sequence in rank selection) | Seeded RNG initialisation; deterministic update arithmetic |
| INV-NMF-06 | Rank 1 ⇒ cophenetic correlation = 1.0 | k = 1 ⇒ all samples in one cluster ⇒ all-ones consensus [6] |
| INV-NMF-07 | Cosine match of a positively-scaled/exact copy of a reference is that reference, cosine 1.0 | Cosine is scale-invariant [5] |

NMF is **non-convex**: the updates only guarantee a *local* minimum, not the global one [2]. Recovery of planted
ground truth (up to a permutation of the k signatures and a positive rescaling between W and H) is therefore only
guaranteed for data that was generated by a non-negative factorisation and, in general, only when that
factorisation is essentially unique (e.g. under a separability / pure-pixel condition).

### 2.5 Comparison with Related Methods

| Aspect | De-novo NMF extraction (`ExtractSignatures`) | NNLS refitting (`FitSignatures`) |
|--------|-----------------------------------------------|----------------------------------|
| Signatures | Discovered from the data (output) | Supplied by caller (input) |
| Unknowns | Both W (signatures) and H (exposures) | Only exposures x |
| Convexity | Non-convex (local optimum, init-dependent) | Convex (unique NNLS minimiser) |
| Input | Count matrix V (channels × samples) | One catalog d + reference signatures S |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `countMatrix` | `IReadOnlyList<IReadOnlyList<double>>` | required | Matrix V: one row per channel, each a vector over samples (`V[channel][sample]`) | non-null, rectangular, finite, ≥ 0; ≥ 1 channel and ≥ 1 sample |
| `rank` | `int` | required | Number of signatures k to extract | 1 ≤ k ≤ channel count |
| `objective` | `NmfObjective` | `Frobenius` (extract) / `KullbackLeibler` (`SelectRank`) | Frobenius or KL/Poisson update variant | enum |
| `maxIterations` | `int` | `10_000` | Max multiplicative-update iterations | > 0 |
| `tolerance` | `double` | `1e-10` | Relative-improvement convergence stop | ≥ 0 |
| `seed` | `int` | `42` | RNG seed (base seed for the derived per-run sequence in `SelectRank`) | any |
| `minRank`,`maxRank` | `int` | required (`SelectRank`) | Candidate-rank range | 1 ≤ minRank ≤ maxRank ≤ channels |
| `runs` | `int` | `20` (`SelectRank`) | NMF restarts per candidate rank | ≥ 1 |
| `stabilityThreshold` | `double` | `0.80` (`SelectRank`) | Min acceptable average stability | [0, 1] |
| `minStability` | `double` | `0.20` (`SelectRank`) | Min acceptable per-signature stability | [0, 1] |
| `extractedSignatures`,`referenceSignatures` | `IReadOnlyList<IReadOnlyList<double>>` | required (`MatchToReferenceSignatures`) | Signatures to label / reference set | non-empty, equal channel count |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Signatures` | `IReadOnlyList<IReadOnlyList<double>>` | k signatures, each a channel vector L1-normalised to sum to 1 |
| `Exposures` | `IReadOnlyList<IReadOnlyList<double>>` | k × samples exposure matrix (`H[j][s]`), non-negative |
| `FinalResidual` | `double` | Final value of the selected objective (‖V − W·H‖²_F, or D(V‖WH) for KL) |
| `Iterations` | `int` | Number of multiplicative-update iterations performed |
| `ObjectiveHistory` | `IReadOnlyList<double>` | Per-iteration objective value (non-increasing) |
| `RankSelectionResult.SelectedRank` | `int` | Auto-selected k (largest qualifying; else highest-average-stability) |
| `RankSelectionResult.PerRank` | `IReadOnlyList<RankStability>` | Per-k cophenetic, average/minimum stability, mean reconstruction error |
| `SignatureMatch` | record | Per extracted signature: best reference index + cosine similarity |

### 3.3 Preconditions and Validation

Inputs are 0-based. The count matrix must be non-null, have at least one channel (row) and one sample (column),
be rectangular (all rows share the sample count), and contain only finite, non-negative entries; otherwise
`ArgumentException` (or `ArgumentNullException` for null matrix/rows) is thrown. `rank` must be in `[1, channels]`;
`maxIterations` must be > 0; `tolerance` must be ≥ 0. The method has no notion of a DNA alphabet — channels are
opaque indices, so any channel layout (SBS-96 or otherwise) is accepted.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate the count matrix and parameters; copy V into a dense jagged array.
2. Initialise W (channels × k) and H (k × samples) with seeded non-negative random values floored above zero.
3. Repeat up to `maxIterations`:
   a. Update H: `H ← H ⊙ (WᵀV) ⊘ (WᵀWH)`.
   b. Update W: `W ← W ⊙ (VHᵀ) ⊘ (WHHᵀ)`.
   c. Record the objective ‖V − WH‖²_F; stop early when its relative decrease < `tolerance`.
4. L1-normalise each signature column of W to sum to 1, absorbing the scale into the matching row of H.
5. Return signatures (column-major → per-signature vectors), exposures, final residual, iteration count, history.

For the **KL/Poisson** objective, steps 3a/3b use the Theorem-2 updates and the objective recorded is D(V‖WH).

For **automatic rank selection** (`SelectRank`): for each k in [minRank, maxRank], run the above NMF `runs`
times with derived seeds; build the consensus matrix from the per-run argmax-metagene assignments and compute its
cophenetic correlation; compute per-signature silhouette stability across runs and the mean reconstruction error;
choose the largest k meeting the stability thresholds (else the highest-average-stability k); return all per-k
diagnostics.

For **reference matching** (`MatchToReferenceSignatures`): for each extracted signature, compute its cosine
similarity to every reference and emit the index + cosine of the maximum.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Denominator floor:** the update denominators (WᵀWH, WHHᵀ) and the random initialisation are floored by
  `NmfEpsilon = 1e-12` to avoid 0/0 when a row/column collapses to zero [2].
- **Convergence test:** relative improvement `(prevObjective − objective) / max(prevObjective, 1)`; the objective
  is non-increasing (INV-NMF-03), so the decrease is non-negative.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| One iteration | O(m·k·n) | O(m·n) | dense W·H product dominates; m = channels, k = rank, n = samples |
| Full extraction | O(I·m·k·n) | O(m·n + m·k + k·n) | I = iterations to convergence (≤ maxIterations) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [OncologyAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs)

- `OncologyAnalyzer.ExtractSignatures(countMatrix, rank, [objective], maxIterations, tolerance, seed)`: de-novo
  NMF extraction; the original 5-arg overload (Frobenius) is preserved and now delegates to the objective overload.
- `OncologyAnalyzer.NmfObjective`: `Frobenius` (Theorem 1) or `KullbackLeibler` (Theorem 2 / Poisson).
- `OncologyAnalyzer.SelectRank(countMatrix, minRank, maxRank, …)`: automatic rank selection; returns
  `RankSelectionResult` with per-rank `RankStability` (cophenetic, average/minimum stability, mean error).
- `OncologyAnalyzer.MatchToReferenceSignatures(extracted, references)`: cosine matching → `IReadOnlyList<SignatureMatch>`.
- `OncologyAnalyzer.SignatureExtractionResult`: the extraction result record.

### 5.2 Current Behavior

Extraction runs at a caller-specified rank, or k is chosen automatically by `SelectRank`. The caller picks the
objective: Frobenius (default for `ExtractSignatures`, byte-for-byte unchanged from the original implementation)
or KL/Poisson (default for `SelectRank`, the SigProfiler choice). `FinalResidual`/`ObjectiveHistory` report the
selected objective's value. The signatures returned are L1-normalised (sum to 1); the exposures absorb the
per-signature scale, so reconstructing W·H reproduces V up to the convergence residual. With the default
relative-change tolerance the solver may stop before full stationarity (sufficient for production fits); for
exact-reconstruction / planted-recovery verification the tests drive convergence with `tolerance = 0`.
`SelectRank` uses a deterministic per-(rank, run) derived seed so the whole procedure is reproducible. No
substring search/matching is involved (signature "matching" is numeric cosine similarity, not text search), so
the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- V ≈ W·H non-negative factorisation with W = signatures, H = exposures [1][3].
- Lee & Seung (2001) Theorem 1 Frobenius and Theorem 2 KL/Poisson multiplicative updates for W and H [1][2][5].
- Monotone non-increase of the selected objective (‖V − WH‖²_F or D(V‖WH)) across iterations [1].
- L1 column-normalisation of signatures to a probability distribution over channels [3][4].
- Brunet (2004) consensus matrix + cophenetic correlation coefficient for rank stability [6].
- Per-signature silhouette stability and the SigProfiler 0.80/0.20 stability + reconstruction-error rank rule [5][7].
- Cosine-similarity matching of extracted signatures to a caller-supplied reference set [5].

**Intentionally simplified:**

- Rank selection uses **greedy best-cosine per-signature** matching of run signatures to a reference partition
  and labels each extracted signature with its single closest reference, rather than SigProfiler's global
  Hungarian assignment; **consequence:** identical for well-separated signatures; for near-degenerate references
  the per-signature labels may not form a global one-to-one assignment.

**Not implemented:**

- Embedded COSMIC reference profiles; **users should rely on:** supplying the COSMIC SBS matrix to
  `MatchToReferenceSignatures` (the library deliberately does not hardcode reference catalogs, matching the
  `FitSignatures` convention).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Caller chooses Frobenius vs KL/Poisson objective | Assumption | KL matches SigProfiler; Frobenius differs on low counts | accepted | ASM-NMF-02; both are Lee & Seung Theorems 1/2 |
| 2 | Local-optimum / init dependence | Deviation | Output depends on seed | accepted | NMF is non-convex [2]; fixed seed for reproducibility |
| 3 | Greedy per-signature cosine matching (not global Hungarian) | Deviation | Labels may not be globally one-to-one for near-degenerate references | accepted | per-signature reduction of SigProfiler matching [5] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null matrix / row | `ArgumentNullException` | input contract |
| Empty matrix (0 channels) or 0 samples | `ArgumentException` | a factorisation needs both dimensions |
| Ragged rows | `ArgumentException` | V must be rectangular |
| Negative or non-finite entry | `ArgumentException` | V must be non-negative and finite |
| rank < 1 or rank > channels | `ArgumentException` | k must be a valid inner dimension |
| maxIterations ≤ 0 or tolerance < 0 | `ArgumentException` | iteration/convergence parameters |
| Exactly factorable V = W₀·H₀ | residual → 0; planted signatures recovered up to perm/scale | INV-NMF-04 [2] |
| All-zero V | factors remain near zero; residual 0 | trivially factorable |

### 6.2 Limitations

- Local optimum only (non-convex); a poor seed can yield a worse factorisation.
- Rank selection cost is O(runs × ranks) full NMF fits; raise `runs` for production-grade stability estimates.
- Reference matching is greedy per-signature (not the global Hungarian assignment) and the reference set must be
  caller-supplied — COSMIC profiles are not embedded.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// V: 96 channels x N samples of SBS counts.
var result = OncologyAnalyzer.ExtractSignatures(countMatrix, rank: 3);
foreach (var signature in result.Signatures)
{
    // signature is a 96-channel probability distribution (sums to 1).
}
// result.Exposures[j][s] = activity of signature j in sample s.
```

### 7.2 Applications and Use Cases

- **Discovering novel mutational processes** in a cohort where no reference signature set fits the data well.
- **Cohort decomposition** into a small number of latent signatures plus their per-sample activities.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [OncologyAnalyzer_ExtractSignatures_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_ExtractSignatures_Tests.cs) — Frobenius path, covers `INV-NMF-01`…`INV-NMF-05`
- Tests: [OncologyAnalyzer_SelectRank_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Oncology/OncologyAnalyzer_SelectRank_Tests.cs) — KL objective, rank selection, cosine matching; covers `INV-NMF-03`, `INV-NMF-06`, `INV-NMF-07`
- Evidence: [ONCO-SIG-002-Evidence.md](../../../docs/Evidence/ONCO-SIG-002-Evidence.md)
- Related algorithms: [Mutational_Signature_Fitting](Mutational_Signature_Fitting.md), [Mutational_Signature_Exposure_Bootstrap](Mutational_Signature_Exposure_Bootstrap.md), [SBS96_Trinucleotide_Context_Catalog](SBS96_Trinucleotide_Context_Catalog.md)

## 8. References

1. Lee D.D., Seung H.S. (2001). Algorithms for Non-negative Matrix Factorization. Advances in Neural Information Processing Systems 13 (NIPS 2000). https://papers.nips.cc/paper/1861-algorithms-for-non-negative-matrix-factorization
2. Non-negative matrix factorization. Wikipedia. https://en.wikipedia.org/wiki/Non-negative_matrix_factorization
3. Alexandrov L.B., Nik-Zainal S., Wedge D.C., Campbell P.J., Stratton M.R. (2013). Deciphering Signatures of Mutational Processes Operative in Human Cancer. Cell Reports 3(1):246–259. https://doi.org/10.1016/j.celrep.2012.12.008
4. Alexandrov L.B. et al. (2020). The repertoire of mutational signatures in human cancer. Nature 578:94–101. https://doi.org/10.1038/s41586-020-1943-3 ; COSMIC SBS96: https://cancer.sanger.ac.uk/signatures/sbs/sbs96/
5. Islam S.M.A., Díaz-Gay M., Wu Y., et al. (2022). Uncovering novel mutational signatures by de novo extraction with SigProfilerExtractor. Cell Genomics 2(11):100179. https://doi.org/10.1016/j.xgen.2022.100179 ; https://github.com/AlexandrovLab/SigProfilerExtractor
6. Brunet J-P., Tamayo P., Golub T.R., Mesirov J.P. (2004). Metagenes and molecular pattern discovery using matrix factorization. PNAS 101(12):4164–4169. https://doi.org/10.1073/pnas.0308531101
7. Rousseeuw P.J. (1987). Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. J. Comput. Appl. Math. 20:53–65. https://doi.org/10.1016/0377-0427(87)90125-7
