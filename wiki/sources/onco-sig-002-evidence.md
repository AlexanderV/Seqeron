---
type: source
title: "Evidence: ONCO-SIG-002 (Mutational signature fitting / refitting + de-novo NMF extraction)"
tags: [validation, oncology]
doc_path: docs/Evidence/ONCO-SIG-002-Evidence.md
sources:
  - docs/Evidence/ONCO-SIG-002-Evidence.md
source_commit: 8cb9903ebbef9fdde3b5e9acf7c72f013a66e74b
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: ONCO-SIG-002

The validation-evidence artifact for test unit **ONCO-SIG-002** — the **downstream mutational-signature
fitting / deconvolution** step that consumes the 96-channel spectrum built by the ONCO-SIG-001 catalog
([[sbs96-mutational-signature-catalog]]) and decomposes it into signature exposures. It is the **thirtieth
ingested unit of the Oncology family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The distinct method is synthesized in its own
concept, [[mutational-signature-fitting-and-extraction]]; [[test-unit-registry]] tracks the unit.

This unit is the **separate downstream concern** that [[sbs96-mutational-signature-catalog]] deferred: given
the SBS-96 count vector, either **refit** it against a set of known COSMIC reference signatures (supervised
NNLS) or **extract** signatures de-novo (unsupervised NMF), then quality-check the reconstruction and match
extracted signatures back to references.

## What this file records

The Evidence file was built in three dated layers (see its Change History):

- **2026-06-14 — NNLS refitting (supervised):**
  - **MutationalPatterns** (Blokzijl 2018, Genome Medicine 10:33; rank 3) — refitting is
    `min_x ‖S·x − d‖₂², x ≥ 0` (S = signature matrix, x = exposure/contribution vector, d = the 96-count
    sample vector), a non-negative least-squares problem solved with an active-set method; and the
    reconstruction quality check is the **cosine similarity** between the observed profile and the
    reconstruction `S·x`, with **≥ 0.95** indicating a successful reconstruction.
  - **deconstructSigs** (Rosenthal 2016, Genome Biology 17:31; rank 3) — reconstruction `S·W`, residual
    `R = T − S·W` minimised by sum-squared error; non-negativity ("negative contributions make no biological
    sense"); weights normalised between 0 and 1 → signature **proportions**.
  - **Lawson & Hanson (1974) active-set NNLS** (via Wikipedia; rank 4 citing the primary) — the verbatim
    active-set algorithm (`argmin_x ‖Ax − y‖₂², x ≥ 0`, passive/active set P/R, gradient `w = Aᵀ(y − Ax)`,
    inner clamping loop); when the unconstrained least-squares solution is already non-negative, NNLS returns
    it unchanged.
  - **iMutSig** (Pan & Wang 2020; rank 3) — cosine similarity `CS(P,C) = P·C / (‖P‖·‖C‖) ∈ [0,1]` as the
    standard pairwise signature-comparison metric.

- **2026-06-23 — de-novo NMF extraction (unsupervised):**
  - **Lee & Seung (2001) NMF** (NIPS 13; rank 1) — the multiplicative-update rules for both the **Frobenius**
    objective `‖V − WH‖²` (Theorem 1: `H ← H·(WᵀV)/(WᵀWH)`, `W ← W·(VHᵀ)/(WHHᵀ)`) and the **generalized KL /
    Poisson** objective `D(V‖WH)` (Theorem 2), each **monotonically non-increasing** each iteration; ratios of
    current values preserve non-negativity; fixed point (all-ones factors) at exact factorization `V = WH`.
  - **Alexandrov 2013** (Cell Reports 3:246; rank 1) — NMF as **blind source separation**: signatures = latent
    sources **W**, per-sample contributions/exposures = mixing weights **H**, so `V ≈ WH`.
  - **Alexandrov 2020 / COSMIC SBS96** (rank 1–2) — 96 channels (6×4×4); each COSMIC signature is an
    **L1-normalized probability distribution** over the 96 types (columns sum to 1).

- **2026-06-23 — rank selection + KL objective + reference matching (enhancement):**
  - **Brunet 2004** (PNAS 101:4164; rank 1 + rank 3 impl docs) — **consensus clustering**: connectivity matrix
    (`C_ij = 1` iff samples i,j share a cluster, cluster = argmax over H columns), consensus = mean
    connectivity across runs, **cophenetic correlation coefficient** (Pearson between consensus-induced and
    hierarchical-cophenetic distances; = 1 for a perfect consensus); rank rule "select the first rank where the
    cophenetic coefficient begins to fall".
  - **SigProfilerExtractor** (Islam 2022; rank 3) — per-signature **stability = silhouette width** of the
    replicate-factorization cluster; solution stable if average stability > 0.80 and none < 0.20; de-novo
    signatures matched to COSMIC by **cosine similarity** (Hungarian assignment globally; a per-signature
    greedy best-cosine reduction is used here).
  - **Rousseeuw 1987** — silhouette `s(i) = (b−a)/max(a,b) ∈ [−1,1]`.

- **Documented corner cases / failure modes:** zero observed catalog `d = 0` → all exposures 0, reconstruction
  0; cosine similarity **undefined for a zero-norm vector** (division by 0) → returns 0.0 and is documented;
  identical vectors → 1, orthogonal/disjoint-support → 0; NMF non-convexity → only a **local minimum**
  (initialization-dependent), so ground-truth recovery is only guaranteed on planted-truth data with a fixed
  seed; division-by-zero guard (ε added to `WᵀWH`/`WHHᵀ` denominators); **permutation/scale ambiguity** (NMF
  factors identified only up to component permutation + positive diagonal rescaling — fixed by L1-normalizing W
  and absorbing scale into H); rank-1 cophenetic = 1 trivially (uninformative); KL with `V_ij = 0` handled via
  `x log x → 0` and floored `(WH)_ij`.

- **Datasets (deterministic worked oracles):** cosine worked values (identical [1,2,3]→1, orthogonal [1,0]/[0,1]
  →0, [1,1]/[1,0]→1/√2 = 0.70710678…, scale-invariance [3,4] vs [6,8]→1); NNLS worked values (identity S,
  d=[3,5]→x=[3,5] residual 0; the **constraint-binding** case S=[[1,1],[0,1]], d=[0,1] where the unconstrained
  x₁=−1 is clamped to 0 and signature 2 refit alone → **x=[0, 0.5]**; zero catalog→[0,0]; normalised exposures
  of [3,5]→[0.375, 0.625]); planted-truth NMF (exactly-factorable `V = W₀H₀` reconstructs with residual ≈ 0;
  planted 96-channel k=2 recovered up to permutation/scaling, columns L1-sum to 1); planted KL k₀=2 recovery +
  auto-rank picks k=2; cosine reference-matching (scaled/permuted reference → cosine 1.0 to itself).

- **Coverage recommendations:** MUST test cosine exact values + scale-invariance; NNLS exact/clamp-and-refit +
  non-negativity; reconstruction `S·x` + cosine = 1 for a representable catalog + proportions sum to 1; NMF
  exactly-factorable → residual ≈ 0, planted recovery up to permutation/scaling, non-negativity + L1
  column-normalization, Frobenius/KL objective monotone non-increasing; auto-rank selects k₀ + reports
  stability + error; cophenetic = 1 at rank 1; cosine matching maps scaled/permuted reference to itself;
  SHOULD test determinism (fixed seed) + input validation (null/empty/ragged/negative, k ≤ 0, k > channels,
  maxIterations ≤ 0, runs < 1, empty reference).

## Deviations and assumptions

- **ASSUMPTION — Frobenius (squared-Euclidean) objective for the core NMF.** SigProfiler's final extraction
  uses a Poisson/KL objective; Lee & Seung give both with proven monotone non-increase, and this
  implementation uses the **Frobenius** objective (Theorem 1) as a faithful, citable NMF, with the KL/Poisson
  variant additionally implemented for the enhancement layer. Both factor the same `V ≈ WH` model — the choice
  invents no constant.
- **ASSUMPTION — deterministic seeded initialisation.** Lee & Seung prescribe no init; a fixed RNG seed
  (mirroring the repo `DefaultBootstrapSeed = 42`) makes results reproducible.
- **ASSUMPTION — exposure normalisation to proportions.** Raw NNLS exposures (counts) are exposed **and** a
  proportion form `exposures / Σexposures` (sums to 1 when Σ > 0, all-zero when Σ = 0); a presentation form
  that does not change the fitted exposures.
- **ASSUMPTION — cosine of a zero vector returns 0.0** (no source defines it; degenerate empty-catalog case
  only).
- **ASSUMPTION — consensus clustering = argmax over H + average-linkage cophenetic; per-signature stability =
  average silhouette (cosine distance); matching = greedy best-cosine per extracted signature** (the
  per-signature reduction of SigProfiler's global Hungarian assignment). Taken from the nimfa / renozao-NMF
  reference implementations of Brunet 2004 (the PNAS full text is gated).

No source contradictions — MutationalPatterns, deconstructSigs, Lawson-Hanson, Lee & Seung, Alexandrov,
Brunet, and SigProfiler are mutually consistent on the fitting model, the NNLS/NMF objectives, and the
cosine/stability quality gates.
