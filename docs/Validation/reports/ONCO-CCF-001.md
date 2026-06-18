# Validation Report: ONCO-CCF-001 — Cancer Cell Fraction Estimation and CCF Clustering

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.EstimateCcf(vaf, purity, tumorCopyNumber, multiplicity)`, `OncologyAnalyzer.ClusterCcfValues(ccfValues, clusterCount)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
- **Tarabichi et al. 2021, *Nature Methods* (PMC7867630), Box 1** — fetched. The PMC HTML flattens the fractions
  ambiguously (renders CCF and m with numerator/denominator scrambled), so it was used only for the *symbol
  inventory* and the clonal-cluster rule, not as the primary numeric authority. It confirms: CCF is "the fraction
  of cancer cells … carrying a set of SNVs, i.e. CCF = CP / purity"; CCF is inferred from VAF (f), purity (ρ),
  local copy-number (N_T = total tumour copy number) and multiplicity (m); normal contributes the `2(1−ρ)` term;
  and **"The cluster with the highest CP can be deemed clonal."**
- **Zheng et al. 2022, *Bioinformatics* (PICTograph)** — fetched. **Verbatim generative equation:**
  `VAF = m·CCF·p / (c·p + 2·(1−p))` (c = tumour copy number, p = purity, normal = 2). This is the clean,
  unambiguous source. Algebraically inverting gives **CCF = VAF·(c·p + 2(1−p)) / (m·p)**.
- **CNAqc CCF-computation vignette (Caravagna lab)** — fetched. Real reference-implementation worked outputs
  (sample purity 89 %, diploid loci): `VAF=0.08, m=1 → CCF=0.180`; `VAF=0.883, m=2 → CCF=0.993`;
  `VAF=0.471, m=1 → CCF=1.06`. Used as an independent cross-check (below). Also documents that the raw CCF can
  exceed 1 from sampling noise (1.06).
- **Lloyd 1982 via Wikipedia K-means** (Stage-B clustering algorithm) — assignment = nearest centroid by least
  squared Euclidean distance; update = mean of members; objective = minimise WCSS. Matches the Evidence doc.
- McGranahan 2016 *Science* (n_mut = VAF·(1/p)·[p·CN_t + CN_n·(1−p)], CCF = n_mut/m) and DeCiFer (Cell Systems)
  were paywalled (403); the formula is fully corroborated by Zheng 2022 + CNAqc, so this did not block validation.

### Formula check
The repo formula `CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m)` is **exactly** the inverse of the verbatim Zheng 2022
equation `VAF = m·CCF·p / (c·p + 2(1−p))` with c = N_T. ✅ Symbols, the normal-diploid `2(1−ρ)` term, and the
`ρ·m` denominator all match. The multiplicity definition in the docs (`m = f·(ρ·N_T+2(1−ρ))/ρ`) is the same
formula with CCF set to 1, which is the correct self-consistent definition. ✅

### Edge-case semantics
- **Domain:** VAF ∈ [0,1], purity ∈ (0,1] (denominator → reject 0), N_T ≥ 1, m ∈ [1, N_T]. All sourced from
  the formula domain and the multiplicity definition (1 ≤ m ≤ tumour copy number, Tarabichi Box 1). ✅
- **CCF > 1 from noise:** cap reported value at 1, expose raw — matches CNAqc (1.06) + McGranahan clonal
  definition (a mutation in all cancer cells has CCF = 1). Documented as an explicit ASSUMPTION; defensible. 🟡→OK
- **Clustering:** highest-centroid cluster = clonal (Tarabichi). Deterministic Lloyd k-means with quantile
  seeding is a sourced, well-defined choice (Lloyd 1982); determinism is an engineering decision, clearly flagged.

### Independent cross-check (numbers retrieved this session)
CNAqc real outputs reproduced by the formula with ρ=0.89, N_T=2 (hand-computed in Python this session):

| CNAqc VAF | m | CNAqc CCF | Formula CCF (ρ=0.89, N_T=2) | Match |
|-----------|---|-----------|------------------------------|-------|
| 0.08  | 1 | 0.180 | 0.180 | ✅ exact |
| 0.883 | 2 | 0.993 | 0.992 | ✅ (rounding) |
| 0.471 | 1 | 1.06  | 1.058 | ✅ (rounding) |

This is an *independent reference-implementation* cross-check — the formula reproduces CNAqc's published numbers.

Hand-computed unit cases (this session): A(0.40,0.80,2,1)=1.0, B(0.20,0.80,2,1)=0.5, C(0.50,1.0,4,2)=1.0,
D(0.25,0.50,2,1)=1.0, E(0.471,1.0,2,1)=0.942, F(0.60,0.80,2,1)=1.5 raw → cap 1.0. All match the TestSpec.

### Findings / divergences
None material. The PMC HTML flattening (noted in the Evidence doc) is a presentation artefact only; the formula
is unambiguously confirmed by Zheng 2022 and CNAqc. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:4880-5049`.

### Formula realised correctly?
`EstimateCcf` (4905-4910): `totalDnaPerCell = purity*tumorCopyNumber + 2*(1-purity)`,
`rawCcf = vaf*totalDnaPerCell/(purity*multiplicity)`, `cappedCcf = Min(1, rawCcf)`. This is the validated
formula verbatim. ✅ Validation guards (4882-4903) match the Stage-A domain exactly: VAF∈[0,1] and N_T<1 →
`ArgumentOutOfRangeException`; purity∉(0,1] → `ArgumentOutOfRangeException`; m∉[1,N_T] → `ArgumentException`;
NaN handled.

`ClusterCcfValues` (4927-5000): null/empty/NaN/Infinity/k-range guards; sort carrying original indices;
quantile seed at `(j+0.5)/k`; Lloyd assignment (nearest by squared distance, 5003-5029) + update
(mean of members, 5032-5049); relabel clusters ascending by centroid; clonal index = k−1 (highest). Matches
Lloyd 1982 + Tarabichi. ✅ Ties go to the lower index (strict `<` in assignment). Empty clusters retain their
centroid (documented). Iteration bound n+1 with convergence break — safe.

### Cross-verification table recomputed vs code
All six EstimateCcf cases and the M11 clustering case were reproduced *independently in Python* this session
(not read from code) and match both the source formula and the test expectations. The new
`ClusterCcfValues_HighCcfValuesFirst` case ({0.94,0.92,0.90,0.14,0.12,0.10}, k=2) was independently computed →
centroids {0.12, 0.92}, assignments [1,1,1,0,0,0], clonal 1 — confirms the ascending relabeling path.

### Variant/delegate consistency
`InferSubclones` and the phylogeny helper consume `ClusterCcfValues` but belong to other units (ONCO-PHYLO-001
etc.); the canonical methods here have no `*Fast`/instance variants. No inconsistency.

### Test quality audit (HARD gate)
- **Sourced expectations, not code echoes:** every EstimateCcf expected value traces to the source formula
  (hand-computed) and, for the new test, to CNAqc's *published* numbers (0.180/0.993/1.06) — a wrong
  implementation would fail them. ✅
- **No green-washing:** exact equality within 1e-10 used throughout where an exact value is known; the only
  inequality (`Is.GreaterThan`, S2 monotonicity) is a genuine property test (INV-2), not a softened exact check;
  the CNAqc cross-check uses 5e-3 only because the published values are rounded to 3 sig figs (justified). No
  skipped/ignored/weakened assertions. ✅
- **Cover all the logic:** both public methods; all four EstimateCcf validation branches + cap + raw + zero-VAF;
  clustering null/empty/NaN/**Infinity**/k-range, k=1, k=2 with exact centroids/assignments, determinism under
  shuffle, and **non-trivial ascending relabeling** (high values first). ✅
- **Added this session (coverage hardening):** (1) `EstimateCcf_CnaqcWorkedOutputs_MatchReferenceImplementation`
  — locks the independent CNAqc reference numbers incl. the >1 raw / capped-1.0 case; (2)
  `ClusterCcfValues_InfiniteValue_Throws` — covers the Infinity validation branch (previously only NaN);
  (3) `ClusterCcfValues_HighCcfValuesFirst_ClonalIsHighestCentroidAscending` — exercises the relabel-to-ascending
  path that the original tests covered only trivially.
- **Honest green:** FULL unfiltered suite = **6680 passed, 0 failed** (1 pre-existing unrelated benchmark skipped);
  `dotnet build` 0 errors; changed test file builds warning-free.

### Findings / defects
None. The implementation faithfully realises the validated formula and clustering algorithm; tests are sourced
and now cover all branches. **Stage B: PASS.**

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. End-state: ✅ CLEAN.**
- No defects logged. Test coverage hardened with 3 additional sourced/branch tests (19 → 22).
- Test-quality gate: **PASS** (sourced, no green-washing, full-logic coverage, honest green 6680/0).
