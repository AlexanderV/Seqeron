# Validation Report: ONCO-SIG-002 — De-novo Mutational Signature Extraction via NMF (+ rank selection + COSMIC matching)

- **Validated:** 2026-06-24   **Area:** Oncology
- **Canonical method(s):**
  `OncologyAnalyzer.ExtractSignatures(countMatrix, rank, [NmfObjective], maxIterations, tolerance, seed)`,
  `OncologyAnalyzer.SelectRank(countMatrix, minRank, maxRank, …)`,
  `OncologyAnalyzer.MatchToReferenceSignatures(extracted, references)`,
  records `SignatureExtractionResult`, `RankStability`, `RankSelectionResult`, `SignatureMatch`, enum `NmfObjective`.
  (NNLS refitting `FitSignatures`/`CosineSimilarity`/`ReconstructCatalog` already validated CLEAN 2026-06-16; this
  session re-validates the de-novo NMF additions from commits `2bcfb0ff` + `2a0e3ac7`.)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

This session targets the de-novo extraction layer added by the limitations campaign:
- `2bcfb0ff` — `ExtractSignatures` (Lee & Seung 2001 Frobenius multiplicative updates; L1-normalised signatures).
- `2a0e3ac7` — KL/Poisson objective (Theorem 2) overload, automatic `SelectRank` (Brunet 2004 cophenetic +
  SigProfiler silhouette stability + reconstruction error), and `MatchToReferenceSignatures` (cosine vs COSMIC).

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

| Source | URL | What it confirmed |
|--------|-----|-------------------|
| Non-negative matrix factorization (Wikipedia, cites Lee & Seung) | https://en.wikipedia.org/wiki/Non-negative_matrix_factorization | Frobenius MU rules `H ← H ⊙ (WᵀV) ⊘ (WᵀWH)`, `W ← W ⊙ (VHᵀ) ⊘ (WHHᵀ)`; squared-Euclidean and a KL-divergence objective are the two Lee & Seung variants; iterate to convergence. |
| NMF::cophcor (CRAN, cites Brunet 2004) | https://search.r-project.org/CRAN/refmans/NMF/html/cophcor.html | Cophenetic correlation = Pearson correlation between consensus-induced sample distances and the cophenetic distances of a hierarchical clustering (default **average** linkage); consensus treated as a similarity ⇒ distance; value 1 = perfectly stable clustering. |
| Lee & Seung (2001) Theorems 1 & 2 (Evidence doc §, proof guide arXiv:2501.11341) | per Evidence | Frobenius (Thm 1) and KL (Thm 2) update rules verbatim; both objectives monotonically non-increasing; KL `D(V‖WH)=Σ(V·log(V/WH)−V+WH)`. |
| Brunet 2004 / SigProfiler (Islam 2022) / Alexandrov 2013 (Evidence doc) | per Evidence | Consensus = mean connectivity (sample→argmax metagene); cophenetic rank rule; stability = per-signature silhouette (Rousseeuw 1987 `s=(b−a)/max(a,b)`); accept if avg ≥ 0.80, min ≥ 0.20; de-novo→reference matched by cosine. |

### Formula check (vs source, vs code)

- **Frobenius MU (Theorem 1):** `UpdateH`/`UpdateW` (OncologyAnalyzer.cs:3680, 3711) implement
  `H ← H ⊙ (WᵀV) ⊘ (WᵀWH)` and `W ← W ⊙ (VHᵀ) ⊘ (WHHᵀ)` exactly, with the denominator floored by ε=1e-12 to
  avoid 0/0. Matches Lee & Seung / Wikipedia verbatim. ✅
- **KL/Poisson MU (Theorem 2):** `UpdateHKl`/`UpdateWKl` (3744, 3782) implement
  `H_aμ ← H_aμ·(Σ_i W_ia V_iμ/(WH)_iμ)/(Σ_i W_ia)` and `W_ia ← W_ia·(Σ_μ H_aμ V_iμ/(WH)_iμ)/(Σ_μ H_aμ)` exactly;
  (WH) floored to avoid V/0, column/row sums floored. Matches Theorem 2 verbatim. ✅
- **KL objective** `KullbackLeiblerDivergence` (3820): `Σ(V log(V/WH) − V + WH)` with the `V=0 ⇒ V log V→0` limit
  (adds only `+WH`). Matches the generalized-KL definition. ✅
- **L1 column normalization** `NormalizeSignatureColumns` (3910): each signature column of W divided by its sum,
  scale absorbed into the matching H row so W·H is invariant (COSMIC: each signature is a probability
  distribution over the 96 channels). ✅
- **Consensus / cophenetic** (`AssignSamplesToClusters` 4206 = argmax metagene; `CopheneticCorrelation` 4234 on
  distance `1−consensus`; `AverageLinkageCopheneticDistances` 4293 UPGMA): matches Brunet 2004 / the CRAN/nimfa
  definition; all-ones consensus ⇒ zero-variance distance ⇒ returns 1.0 (Brunet "perfect consensus = 1"). ✅
- **Silhouette stability** `SignatureStability` (4392): clusters the runs×k signatures by greedy one-to-one cosine
  match to run-0, per-point silhouette `(b−a)/max(a,b)` with cosine distance; single-cluster ⇒ s=0 (Rousseeuw
  convention). Average and minimum over signatures. Matches Rousseeuw 1987 + SigProfiler "silhouette of the
  cluster corresponding to that signature." ✅
- **Selection rule** `ChooseRank` (4540): largest k with avg ≥ 0.80 and min ≥ 0.20 (SigProfiler defaults), else
  highest-average-stability (ties → smallest/most-parsimonious). ✅
- **Cosine matching** `MatchToReferenceSignatures` (4609): each extracted signature labelled with the reference of
  maximal cosine — the per-signature reduction of SigProfiler's Hungarian "maximise total cosine." ✅

### Edge-case semantics check

- Rank 1 / <2 samples ⇒ cophenetic = 1.0 (every consensus entry 1, all-ones; Brunet trivial case). ✅
- KL with V=0 entries ⇒ log term vanishes, +(WH) retained; (WH) floored to avoid log(·/0). ✅
- NMF non-convexity ⇒ planted-truth recovery only guaranteed up to permutation/positive scaling on separable
  (pure-channel) factorable data; tests use seeded init + separable W₀,H₀ and run to tight convergence. ✅
- Single run ⇒ no cross-run dispersion ⇒ stability treated as 1.0 (documented). ✅
- Validation: null/empty/ragged/negative/NaN matrix, rank<1 / rank>channels / maxIter≤0 / tol<0, minRank<1 /
  maxRank<minRank / maxRank>channels / runs<1 / threshold∉[0,1], reference null/empty/mismatch — all throw
  `ArgumentNullException` / `ArgumentException`. ✅

### Independent cross-check (numbers — recomputed this session, NOT from repo code)

Independent NumPy reimplementations of the Lee & Seung updates on planted V = W₀·H₀ (separable):

| Check | Independent result | Spec/code expectation |
|-------|--------------------|-----------------------|
| Frobenius residual ‖V−WH‖²_F, exact rank-2 V (W₀=[[5,0],[0,4],[2,1],[1,3]], H₀=[[6,0,3,2],[0,5,1,4]]) | 4.07e-9 | < 1e-6 (M1) ✅ |
| Frobenius objective monotone non-increasing over 200k iters | no violation | INV-NMF-3 ✅ |
| L1-normalised extracted signature column sums | [1.0, 1.0] | sum = 1 (M5) ✅ |
| Planted sig 0 / sig 1 best cosine to extracted | 1.000000 / 1.000000 | ≈1 perm/scale (M3) ✅ |
| KL divergence D(V‖WH) on factorable V | 0.0, monotone (no violation) | INV-KL-1 ✅ |
| KL planted recovery cosine (sig 0 / sig 1) | 0.999237 / 0.995802 | > 0.99 (M-KL2) ✅ |
| cos(5·refA, refA) / cos(5·refA, refB), refA=[.7,.1,.1,.1] | 1.0 / 0.3077 | 1.0 / <0.5 (M-MT1) ✅ |

### Findings / divergences

None material. Documented, sourced assumptions: KL is a caller-selectable objective (Frobenius default preserved
byte-for-byte via overload routing); consensus = argmax-metagene + average-linkage cophenetic (Brunet via
nimfa/NMF-CRAN, PNAS full text gated); stability = average silhouette with cosine distance (SigProfiler /
Rousseeuw); matching = greedy best-cosine per extracted signature. The TestSpec note that Brunet's "first rank
where cophenetic falls" is reported but the *selection* uses the SigProfiler stability threshold is correct and
intentional — cophenetic is exposed in `RankStability` for auditability while `ChooseRank` follows Islam 2022.
Stage A: **PASS**.

## Stage B — Implementation

### Code path reviewed (file:line)

`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- `ExtractSignatures` overloads (3431 Frobenius-default → 3485 objective-explicit); `RunNmf` core loop (3537);
  `ValidateCountMatrix` (3602); `ValidateRankAndStop` (3508); `InitializeNonNegativeFactor` (3659).
- Frobenius `UpdateH`/`UpdateW` (3680/3711); KL `UpdateHKl`/`UpdateWKl` (3744/3782);
  `KullbackLeiblerDivergence` (3820); `FrobeniusResidualSquared` (3884); `MultiplyWh` (3857);
  `NormalizeSignatureColumns` (3910); `TransposeColumnsToSignatures` (3942); `RowsToReadOnly` (3962).
- `SelectRank` (4062) + `DerivedSeed` (4191), `AssignSamplesToClusters` (4206), `CopheneticCorrelation` (4234),
  `AverageLinkageCopheneticDistances` (4293), `SignatureStability` (4392), `ChooseRank` (4540).
- `MatchToReferenceSignatures` (4609) reusing `CosineSimilarity` + `ValidateSignatures`.

### Formula realised correctly? (evidence)

Yes. Every update/objective formula is the cited expression verbatim (see Stage A formula table). The relative-
improvement stop (3575) relies on the proven monotone non-increase (Thm 1/2) and so cannot loop forever; the L1
normalization runs once after convergence and is W·H-invariant. `SelectRank` reuses `RunNmf` per (rank, run) with
a deterministic derived seed, builds the Brunet consensus, computes cophenetic + silhouette stability + mean
error, and selects via the SigProfiler rule. The independent NumPy reruns above reproduce the same residuals,
monotonicity, recovery cosines, and matching values that the actual code produces (all 40 tests pass).

### Cross-verification table recomputed vs code

The Stage-A independent numbers were all reproduced by the actual code through the passing test classes:
`OncologyAnalyzer_ExtractSignatures_Tests` (19: exact reconstruction < 1e-6, product recovers V, planted recovery
cosine, nonnegativity, L1 sum=1, monotonicity, SBS-96, determinism, 10 guards) and
`OncologyAnalyzer_SelectRank_Tests` (21: KL monotonicity + planted recovery > 0.999, KL L1 + determinism,
Frobenius backward-compat, SelectRank picks k0=2, rank-1 cophenetic = 1.0, error vs rank, cosine matching
1.0/cross<0.5, validation guards). Filtered run: **40 passed, 0 failed**.

### Variant/delegate consistency

The Frobenius-default `ExtractSignatures` overload routes to the objective-explicit overload with
`NmfObjective.Frobenius` (3437), so the original behaviour is preserved exactly (test M-KL5 asserts equal
FinalResidual). `SelectRank` and `SignatureStability` reuse `RunNmf` and `CosineSimilarity` — no divergent
re-implementation. `MatchToReferenceSignatures` reuses the validated `CosineSimilarity`.

### Numerical robustness

All denominators floored by ε=1e-12 (MU divisors, (WH) in KL, cophenetic variance, silhouette denom); the
relative-improvement stop guards against infinite loops and is backed by the monotonicity theorem; `Math.Log`
only on a floored positive ratio; nonnegativity preserved by the multiplicative structure (init drawn from
(0,1]+ε so no factor entry starts at 0). No overflow concern on count-scale inputs.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** expectations are planted ground truth and paper formulas — residual < 1e-6
  (Lee & Seung fixed point), recovery cosine ≈ 1 (separability), L1 sum = 1 (COSMIC), monotone objective
  (Thm 1/2), cophenetic = 1 at rank 1 (Brunet), cos(5·ref,ref)=1 (scale invariance). All independently
  reproduced here via NumPy, not read off the implementation. A wrong update rule would break the residual /
  recovery / monotonicity assertions; a wrong matcher would break M-MT1/M-MT2.
- **Deterministic:** seeded init + derived per-(rank,run) seeds; determinism tests assert byte-identical repeats.
- **Coverage:** all M/S/V rows of TestSpec Appendices A & B implemented; 10 + 1x validation-guard set.
- **Honest green:** FULL unfiltered `Seqeron.Genomics.Tests` = **18213 passed, 0 failed, 0 skipped** (1 benchmark
  marked Skipped by design); `dotnet build` 0 warnings / 0 errors.

**Gate result: PASS.**

### Findings / defects

No code defect. No code change was required this session.

## Verdict & follow-ups

- **Stage A: PASS.** Frobenius (Thm 1) and KL/Poisson (Thm 2) multiplicative updates, the KL divergence, L1
  signature normalization, Brunet consensus/cophenetic (average linkage on 1−consensus), Rousseeuw silhouette
  stability with the SigProfiler 0.80/0.20 rule, and cosine reference matching all match Lee & Seung 2001,
  Brunet 2004, Alexandrov 2013/2020 (COSMIC), Islam 2022 and Rousseeuw 1987 verbatim. Assumptions are declared
  and standard.
- **Stage B: PASS.** Code realises every validated formula; independent NumPy cross-checks (Frobenius residual
  4e-9, KL divergence → 0, planted recovery cosine 1.0 / 0.999, cos(5·refA,refA)=1.0 vs cross 0.31) match the
  actual code via 40 passing unit tests; numerically robust; tests assert exact sourced values, not tautologies.
- **End-state: ✅ CLEAN.** No defect; the de-novo NMF extraction, automatic rank selection, and COSMIC cosine
  matching are fully functional. Full suite green (18213/0).
