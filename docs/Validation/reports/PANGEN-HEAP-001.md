# Validation Report: PANGEN-HEAP-001 — Pan-Genome Growth Model (Heaps' law)

- **Validated:** 2026-06-15   **Area:** PanGenome
- **Canonical method(s):** `PanGenomeAnalyzer.FitHeapsLaw(IEnumerable<GenePresenceRow>, int)`;
  `PanGenomeAnalyzer.FitHeapsLaw(IReadOnlyDictionary<...>, double, int)` (delegate);
  `PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(genomes, clusters)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (two test-quality defects found and fully fixed this session)
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened this session (not the repo's citations)

1. **micropan `heaps()` source — `R/powerlaw.R`**
   (https://raw.githubusercontent.com/larssnip/micropan/master/R/powerlaw.R, retrieved 2026-06-15).
   Extracted verbatim:
   - Model: `y.hat <- p[1] * x^(-p[2])`
   - Binarization: `pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1`
   - New-cluster count per ordering: `cm <- apply(pan.matrix[sample(nrow(pan.matrix)),], 2, cumsum)`;
     `nmat[,i] <- rowSums((cm == 1)[2:ng,] & (cm == 0)[1:(ng-1),])`
   - Pooled points: `x <- rep((2:nrow(pan.matrix)), times = n.perm)`, `y <- as.numeric(nmat)`
   - Objective: `J <- sqrt(sum((y - y.hat)^2))/length(x)`
   - Start: `p0 <- c(mean(y[which(x == 2)]), 1)`
   - Optim: `optim(p0, objectFun, ..., method = "L-BFGS-B", lower = c(0, 0), upper = c(10000, 2))`
   - Return names: `c("Intercept", "alpha")`
2. **micropan CRAN refman — `heaps.html`**
   (https://search.r-project.org/CRAN/refmans/micropan/html/heaps.html, retrieved 2026-06-15).
   Open/closed rule verbatim: **"If 'alpha>1.0' the pan-genome is closed, if 'alpha<1.0' it is open."**
   n.perm: "The default value of 100 is certainly a minimum."
3. Tettelin et al. (2005) PNAS and Tettelin et al. (2008) Curr Opin Microbiol are the model's
   primary literature (power-law new-gene-discovery curve; open/closed via the decay exponent).

### Formula check

Every formula in the TestSpec/Evidence/algorithm doc matches the retrieved micropan source
**byte-for-byte**: model `n(N)=K·N^(−α)`, objective `J=sqrt(Σ(y−K·x^(−α))²)/|x|`, box bounds
K∈[0,10000] / α∈[0,2], start `(mean y at x=2, 1)`, first-appearance counting
`(cm==1)[i] & (cm==0)[i−1]` with the index starting at N=2, and binarization `>0→1`.

### Edge-case semantics

- < 2 genomes → empty new-gene curve `x = 2:ng`, no fit → degenerate result. Sourced. ✅
- Index starts at 2 (first genome has no predecessor). Sourced. ✅
- Binary presence; copy number > 1 collapses to 1. Sourced. ✅
- α boundary at 1.0: micropan states strict `>1` closed / `<1` open and is silent on exactly 1.0;
  the spec's convention "α=1 ⇒ not-open (closed)" is a reasonable, documented tie-break. ✅ (minor)

### Independent cross-check (numbers recomputed this session, not from the repo)

| Case | Curve | Recomputed expected | Source |
|------|-------|--------------------|--------|
| M1 (closed) | y=[8,4] | α = ln2/ln(3/2) = **1.709511291351455**, K = 8·2^α = **26.1640013949735**, closed | micropan model + hand calc |
| M2 (open)   | y=[1,1] | α = 0, K = 1, open | derived: best power fit of a constant, J=0 |
| M4 (open)   | y=[1,1] | α = 0, K = 1, open | first-appearance on g1{core,a,b}/g2{+c}/g3{+d} |
| M4b (closed)| y=[2,1] | α = ln2/ln(3/2), K = 2·2^α = **6.541000348743378**, closed | first-appearance, "x present g1+g3 absent g2 not recounted" |

All values reproduced by independent computation; M1 matches the spec to the recorded digits.

### Stage A findings

No biological/mathematical defect. The description faithfully reproduces the authoritative
micropan/Tettelin model. **Stage A = PASS.**

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs`:
- `FitHeapsLaw(matrix, permutations)` → `BuildPermutedNewGeneCurves` (L607) → `FitHeapsLawFromOrderings` (L678).
- New-gene counting (L647–666): seeds `seen[]` from the first genome, counts at each position ≥ 2 a
  cluster that is present and not yet `seen` — exactly micropan `(cm==1)[i] & (cm==0)[i−1]`.
- Binarization (L624–628): presence read as `bool` from `GenePresence`; structurally `>0→1`.
- Objective (`HeapsObjective`, L719): `sqrt(Σ(y−K·x^(−α))²)/|x|` — verbatim.
- Bounds/start (L698–701, L741–742): K∈[0,10000], α∈[0,2], start `(mean y at x=2, 1)` — verbatim.
- Optimizer: deterministic bounded coordinate descent with geometric refinement (documented
  divergence from R's L-BFGS-B; objective/bounds/start identical). Recovers analytic optimum to <1e-9.
- `EmptyHeapsFit` (L708): `(0,0,false,n→0)` for <2 genomes / null / empty — sourced degenerate case.
- Open/closed (L703): `isOpen = α < 1.0` — matches the strict micropan rule.
- `CreatePresenceAbsenceMatrix` (L352): `present = cluster.GeneIds.Any(...)`; `PresentGenes` counts
  distinct present clusters — correct binarization.

(Note: `DeterminePanGenomeType`/`EstimateHeapsDecayExponent` use a different log-log regression for
the `ConstructPanGenome` openness verdict; that path belongs to PANGEN-CONSTRUCT, not this unit's
canonical methods, and was not in scope here.)

### Formula realised correctly?

Yes. Each cited formula/bound/start/counting-rule is realised verbatim. The optimizer-method
divergence is documented and verified non-correctness-affecting on in-bounds power-curve data
(M1/M4b recovered to machine precision).

### Cross-verification recomputed vs code

All recomputed values (M1, M2, M4, M4b, M5/M5b) reproduced by the running code to ≤1e-6 (params)
and 1e-9 (analytic optima). 13 tests in the canonical file pass.

### Variant/delegate consistency

M9 (`FitHeapsLaw_DictionaryOverload_DelegatesToMatrixFit`) confirms the dictionary overload produces
identical (Intercept, Alpha, IsOpen) to the matrix path. ✅

### Test quality audit (HARD gate)

Two defects found in the existing canonical test file — **both test defects, no code defect**:

1. **M5 binarization test was tautological.** It asserted `Count(v=>v)==GenePresence.Count` and
   `All(v => v==true||v==false)` — the latter is always true for a `bool`, and the former held by
   construction. The test could not fail for *any* binarization bug → green-wash by tautology.
   **Fix:** rewrote `CreatePresenceAbsenceMatrix_DuplicatePresence_CountsOnce` to a genuine
   multi-member cluster (two identical-sequence genes → one cluster) asserting `PresentGenes==1`
   (not 2) and exactly one present flag; added `FitHeapsLaw_BinarizesMultiMemberPresence` (M5b)
   that exercises binarization in the FIT path (duplicate shared-cluster members → new count 1).

2. **M4 was redundant with M1.** It reused M1's exact matrix and asserted only `predictor(2)=8`,
   `predictor(3)=4` — algebraically identical to M1's parameter assertions, so it added no
   independent coverage of the first-appearance counting rule, while its comment claimed a
   counterfactual ("if shared genes were recounted the fit would differ") that was never executed.
   **Fix:** rewrote `FitHeapsLaw_CountsNewGenesByFirstAppearance` to a *distinguishing* matrix
   (first-appearance y=[1,1] → α=0/open; a naive recount would give y=[3,4] → different verdict) and
   added `FitHeapsLaw_FirstAppearance_NotImmediatePredecessor` (M4b: a cluster present at g1, absent
   g2, present g3 must NOT be recounted at g3 → exact y=[2,1], α=ln2/ln1.5, K=2·2^α, closed; a
   wrong "compare to immediate predecessor" impl would yield y=[2,2]/open).

Other cases (M1, M2, M3 folded into M1/M2, M6, M7, M8, M9, S1, S2, C1) assert exact sourced values
with tight tolerances; S2 uses `Is.InRange` correctly (bounds are a range invariant, INV-04/05).
No skipped/ignored/weakened assertions remain.

**Gate result: PASS** — after fixes, every Stage-A branch, formula path, and documented
edge/error case is covered by sourced (not code-echoed) assertions; full unfiltered suite green.

### Honest-green evidence

- `dotnet build` — 0 warnings, 0 errors.
- Canonical file: **13 passed, 0 failed** (was 11; +M4b, +M5b; M4 and M5 rewritten).
- FULL unfiltered `Seqeron.Genomics.Tests`: **6545 passed, 0 failed**, 0 skipped (the single
  `MFE_Benchmark_AllScenarios` skip is a pre-existing, unrelated benchmark guard).

### Findings / defects

- Two test-quality defects (above) — **fully fixed this session.** No implementation defect.

## Verdict & follow-ups

- **Stage A: PASS.** Description matches the retrieved micropan source and CRAN docs verbatim;
  M1 analytic values reproduced independently.
- **Stage B: PASS-WITH-NOTES.** Implementation is faithful to the validated description; two
  green-washing test defects (tautological M5, redundant M4) were found and completely fixed with
  exact sourced expectations and a distinguishing counterfactual.
- **End-state: ✅ CLEAN.** No outstanding gap; algorithm fully functional, full suite green.
