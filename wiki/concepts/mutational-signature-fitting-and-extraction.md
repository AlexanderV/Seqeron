---
type: concept
title: "Mutational signature fitting (NNLS refit) and de-novo extraction (NMF)"
tags: [oncology, algorithm]
sources:
  - docs/Evidence/ONCO-SIG-002-Evidence.md
  - docs/algorithms/Oncology/Mutational_Signature_Fitting.md
  - docs/algorithms/Oncology/Mutational_Signature_Extraction_NMF.md
source_commit: c559752e2d7524d5b9f60057aa6b5798d97e6835
created: 2026-07-10
updated: 2026-07-14
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: onco-sig-002-evidence
      evidence: "Test Unit ID: ONCO-SIG-002 ... Algorithm: Mutational Signature Fitting / Refitting (NNLS decomposition + cosine similarity)"
      confidence: high
      status: current
    - predicate: depends_on
      object: concept:sbs96-mutational-signature-catalog
      source: onco-sig-002-evidence
      evidence: "NNLS refitting: d = the original 96-mutation count vector for a sample; NMF V has 96 channels (SBS-96) — the fitting decomposes the SBS-96 spectrum built by ONCO-SIG-001."
      confidence: high
      status: current
---

# Mutational signature fitting and de-novo extraction

The Oncology family's **signature-deconvolution** unit (**ONCO-SIG-002**): the downstream step that takes the
96-channel single-base-substitution spectrum built by the catalog step
([[sbs96-mutational-signature-catalog]], ONCO-SIG-001) and decomposes it into **signatures × exposures**. It
answers the two complementary questions the catalog deliberately left open:

- **Refit (supervised):** given a set of **known** COSMIC reference signatures, how much of each is present in
  this sample? → **non-negative least squares (NNLS)**.
- **Extract (unsupervised):** without a reference, what latent signatures generated this cohort's spectra? →
  **non-negative matrix factorization (NMF)**.

The literature-traced record is [[onco-sig-002-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the evidence-artifact pattern.

## Supervised refit — NNLS

Fit the observed 96-count vector `d` against a signature matrix `S` (channels × known signatures) by minimising

`min_x ‖S·x − d‖₂²  subject to  x ≥ 0`

where `x` is the **exposure** (contribution) vector. Non-negativity is biological — a signature cannot
contribute a negative number of mutations. Solved by the **Lawson-Hanson active-set** algorithm
(passive/active set, gradient `w = Aᵀ(y − Ax)`, an inner loop that clamps any coefficient driven ≤ 0 and
refits the remaining ones). When the plain least-squares solution is already non-negative, NNLS returns it
unchanged.

The **reconstruction** is `S·x`; its quality is the cosine similarity to the observed `d`, with **≥ 0.95** the
conventional "successful reconstruction" gate (MutationalPatterns). A cosine below the gate signals a
mutational process not represented in `S`, not an error. Raw exposures (counts) and a **proportion** form
`x / Σx` (deconstructSigs "weights between 0 and 1") are both reported.

**Worked NNLS oracles:**

| S | d | Fitted x | Reconstruction | Note |
|---|---|----------|----------------|------|
| identity `[[1,0],[0,1]]` | `[3,5]` | `[3,5]` | `[3,5]`, residual 0, cos 1 | exact |
| `[[1,1],[0,1]]` | `[0,1]` | `[0, 0.5]` | `[0.5,0.5]` | unconstrained `x₁=−1` clamped to 0, sig 2 refit alone |
| identity | `[0,0]` | `[0,0]` | `[0,0]` | zero catalog → zero exposures |

The clamp-and-refit case is the defining NNLS behaviour: the unconstrained normal equations give `x = [−1,1]`;
`x₁ < 0` so signature 1 is removed and signature 2 is refit alone, `x₂ = (s₂·d)/(s₂·s₂) = 1/2`.

The **uncertainty** on this refit — a per-signature confidence interval rather than a bare point estimate — is
the separate ONCO-SIG-003 unit [[signature-exposure-bootstrap-confidence-intervals]], which resamples the
catalog and re-runs this exact NNLS fit per replicate. The **interpretation** of these exposures — mapping the
normalized proportions to named **mutational processes** (APOBEC, Aging, Tobacco, UV, MMR) under a 6% presence
cutoff — is the downstream ONCO-SIG-004 unit [[mutational-process-classification]].

## Cosine similarity — the shared quality/comparison metric

`CS(A,B) = A·B / (‖A‖·‖B‖) ∈ [0,1]` — the cosine of the angle between two non-negative profile vectors. Used
both as the reconstruction gate above and as the pairwise signature-comparison metric.

- **Identical** vectors → 1 exactly; **orthogonal / disjoint support** → 0; **scale-invariant** (`[3,4]` vs
  `[6,8]` → 1). `[1,1]` vs `[1,0]` → `1/√2 = 0.70710678…`.
- **Zero vector:** cosine is undefined (division by a zero norm). The implementation returns **0.0** for a
  zero-norm pair (no shared direction) — the only place this bites is the degenerate empty-catalog case.

## Unsupervised extraction — NMF (`V ≈ W·H`)

Factor a catalog matrix `V` (channels × samples) into a **signature** matrix `W` (channels × k) and an
**exposure** matrix `H` (k × samples) — Alexandrov's **blind source separation** framing: signatures are the
latent sources `W`, per-sample exposures the mixing weights `H`. Solved by **Lee & Seung multiplicative
updates**, which keep both factors non-negative (each update is a ratio applied to the current value) and drive
the objective monotonically down:

- **Frobenius** objective `‖V − WH‖²` (Theorem 1, the default here): `H ← H·(WᵀV)/(WᵀWH)`,
  `W ← W·(VHᵀ)/(WHHᵀ)`.
- **Generalized KL / Poisson** objective `D(V‖WH)` (Theorem 2, the enhancement layer, the objective SigProfiler
  uses): the ratio updates weighted by `V/(WH)`.

Both are **monotonically non-increasing each iteration** — the verification handle. At an exact factorization
`V = WH` the multiplicative factors are all-ones, so the iterates stop (residual ≈ 0). COSMIC convention:
**L1-normalize each signature column of W to sum to 1** (a probability distribution over the 96 channels) and
absorb the removed scale into `H` — this fixes NMF's positive-rescaling ambiguity.

**Key sharp edges of NMF:**

- **Non-convex** → multiplicative updates reach only a **local** minimum, dependent on initialization; recovery
  of a planted ground truth is only guaranteed on data actually generated by a non-negative factorization and
  with a favourable (fixed-seed) init. Tests use planted-truth `V = W₀H₀`.
- **Permutation + scale ambiguity** → factors are identified only up to a component permutation and a positive
  diagonal rescaling; matching is therefore done column-wise before comparing (cosine ≈ 1 after matching).
- **Division-by-zero guard** → a small ε is added to the `WᵀWH` / `WHHᵀ` denominators.

## Rank selection and reference matching

- **How many signatures (rank k)?** Run NMF repeatedly over a range `[kMin, kMax]` and score each rank by
  **stability** and **reconstruction error**:
  - **Consensus clustering / cophenetic correlation** (Brunet 2004): each run's connectivity matrix has
    `C_ij = 1` iff samples i,j share a cluster (cluster = argmax over that sample's H column); the consensus is
    the mean connectivity across runs; the **cophenetic correlation coefficient** is the Pearson correlation
    between the consensus-induced distances and hierarchical (average-linkage) cophenetic distances. Select the
    first rank where it **begins to fall**. Rank 1 is trivially 1.0 (all-ones consensus) and thus uninformative.
  - **Silhouette stability** (SigProfilerExtractor / Rousseeuw): per-signature stability = the average
    silhouette width (`s = (b−a)/max(a,b)`, cosine distance) of the replicate-factorization cluster; a solution
    is stable if average stability ≥ 0.80 and none < 0.20.
  - **Selection rule** (`SelectRank`): the chosen rank is the **largest candidate k** whose average stability
    ≥ `stabilityThreshold` (0.80) with minimum per-signature stability ≥ `minStability` (0.20); if none
    qualifies, the highest-average-stability k is returned. Cophenetic correlation and mean reconstruction
    error are reported per rank as diagnostics but do not drive the choice. Cost is `O(runs × ranks)` full NMF
    fits; a deterministic per-(rank, run) derived seed makes the whole sweep reproducible.
- **Matching de-novo signatures to COSMIC:** label each extracted signature with its **closest reference by
  cosine similarity** (a greedy per-signature reduction of SigProfiler's global Hungarian assignment). A
  scaled or channel-permuted copy of a reference matches it with cosine 1.0. The library **does not embed
  COSMIC profiles** — the reference set is caller-supplied (matching the `FitSignatures` convention).

## API surface and implementation notes

`OncologyAnalyzer` (`Seqeron.Genomics.Oncology`, `OncologyAnalyzer.cs`) exposes the **refit** layer as:

- `CosineSimilarity(a, b)` → cosine similarity of two equal-length non-empty vectors (returns `0.0` when
  either has zero norm).
- `FitSignatures(catalog, signatures)` → `SignatureFitResult` (`Exposures` = the NNLS solution `x ≥ 0`,
  `NormalizedExposures` = `x / Σx` summing to 1 when Σ>0 else all 0, `Reconstruction` = `S·x`,
  `ReconstructionCosineSimilarity` = cosine of `d` vs `S·x`). Each active-set inner iteration solves the
  passive-set normal equations `s_P = ((S_P)ᵀ S_P)⁻¹ (S_P)ᵀ d` by **Gaussian elimination with partial
  pivoting** (dense — tuned for small signature counts, not thousands), gradient tolerance `ε = 1e-12`; a
  **singular passive-set matrix** (collinear signatures) leaves the affected component at 0 rather than
  throwing. Cost `O(k³ + k²·n)` per outer iteration over `≤ O(k)` outer iterations (`k` signatures, `n`
  channels).
- `ReconstructCatalog(signatures, exposures)` → `S·x`.

and exposes the extraction layer as:

- `ExtractSignatures(countMatrix, rank, [objective], maxIterations=10_000, tolerance=1e-10, seed=42)` →
  `SignatureExtractionResult` (`Signatures`, `Exposures`, `FinalResidual`, `Iterations`, `ObjectiveHistory`).
  The 5-arg Frobenius overload is preserved and delegates to the objective overload; `NmfObjective` is
  `Frobenius` (default here) or `KullbackLeibler`.
- `SelectRank(countMatrix, minRank, maxRank, runs=20, stabilityThreshold=0.80, minStability=0.20, …)` →
  `RankSelectionResult` (`SelectedRank`, per-rank `RankStability` = cophenetic, average/minimum stability,
  mean error). Its objective defaults to **KL/Poisson** (the SigProfiler choice), unlike `ExtractSignatures`.
- `MatchToReferenceSignatures(extracted, references)` → `IReadOnlyList<SignatureMatch>` (best reference index
  + cosine per extracted signature).

`countMatrix` is `V[channel][sample]` — row-per-channel, non-negative, rectangular, finite; channels are
**opaque indices** (no DNA-alphabet notion, so any layout, SBS-96 or otherwise, is accepted). Numeric edges:
the update denominators (`WᵀWH`, `WHHᵀ`) and the seeded random init are floored by `NmfEpsilon = 1e-12`; the
convergence test is relative improvement `(prevObjective − objective) / max(prevObjective, 1) < tolerance`.
With the default tolerance the solver may stop before full stationarity (fine for production fits); planted-
recovery tests drive convergence with `tolerance = 0`. One iteration costs `O(m·k·n)` (dense `W·H`), a full
extraction `O(I·m·k·n)`.

## Relation to the oncology family

This is the **downstream partner** of the catalog: [[sbs96-mutational-signature-catalog]] (ONCO-SIG-001) builds
the 96-channel spectrum `d`/`V` (this unit `depends_on` it), and this unit decomposes it into exposures the
interpretation layers can act on. The fitted exposure profile is a somatic-mutation-process biomarker
orthogonal to the copy-number-scar [[homologous-recombination-deficiency-score]] and the mismatch-repair
[[microsatellite-instability-detection]], and it feeds the clinical-interpretation units
([[cancer-variant-tier-classification-amp-asco-cap]], [[clinical-actionability-oncokb-levels]]) — e.g. a strong
SBS signature exposure is exactly the kind of signal those layers consume. The **NNLS refit** and **NMF
extraction** are two alternative deconvolution strategies (supervised vs unsupervised) over the same `V ≈ WH`
model.

## Scope and limitations

A [[scientific-rigor|research-grade]] correctness reference for signature fitting/extraction. The models are
MutationalPatterns / deconstructSigs (NNLS refit), Lawson-Hanson (active-set NNLS), Lee & Seung (Frobenius +
KL NMF), Alexandrov 2013/2020 (blind-source-separation + COSMIC L1 normalization), and Brunet 2004 /
SigProfilerExtractor (rank selection + cosine matching). The core NMF uses the **Frobenius** objective (KL/
Poisson additionally implemented); a **fixed seed** makes it reproducible; exposure **proportions** are a
presentation form over the source-defined NNLS minimiser. **Not for clinical or diagnostic use.** No source
contradictions.
