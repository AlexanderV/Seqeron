# Validation Report: META-PATHWAY-001 — Metabolic Pathway Enrichment (ORA, hypergeometric test)

- **Validated:** 2026-06-15   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.FindPathwayEnrichment(queryGenes, pathwayDatabase, backgroundGenes?)`, `MetagenomicsAnalyzer.HypergeometricUpperTail(x, bigN, bigM, n)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES

## Scope note (prompt vs repo)

The session prompt's parenthetical hints described a *pathway completeness / reconstruction* family
(KEGG module completeness, coverage = present steps / total steps, MinPath / HUMAnN pathway abundance).
The repository defines META-PATHWAY-001 unambiguously and consistently across every artifact
(`ALGORITHMS_CHECKLIST_V2.md` §META-PATHWAY-001 → canonical `FindPathwayEnrichment`, TestSpec, Evidence,
algorithm doc, implementation) as **hypergeometric Over-Representation Analysis (ORA)** — a different,
well-established pathway-enrichment method. No completeness/MinPath/HUMAnN code exists under this unit ID.
This validation therefore assesses the algorithm **as the repo defines it** (ORA). The completeness
concept is a *different algorithm*, not a defect in this one; it is out of scope for META-PATHWAY-001.

## Stage A — Description

### Sources opened & what they confirm
- **Boyle et al. (2004), GO::TermFinder, *Bioinformatics* 20(18):3710–3715**
  (https://pmc.ncbi.nlm.nih.gov/articles/PMC3037731/). Confirms ORA p-value = right tail of the
  hypergeometric distribution; symbol roles (N background, M term/pathway size, n query size, k overlap);
  one-sided "k or more" (upper-tail) semantics; hypergeometric (sampling without replacement) chosen over
  binomial. Cited correctly by the repo.
- **PNNL Proteomics Data Analysis in R/Bioconductor §8.2 ORA**
  (https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html). Confirms the exact
  formula `P(X ≥ x) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M,n−i)/C(N,n)`, the reference R call
  `phyper(q = x−1, m = M, n = N−M, k = n, lower.tail = FALSE)`, and the worked example
  N=8000, M/n={400,100}, x=20 → P(X≥20) ≈ 7.88×10⁻⁸. Cited correctly.

### Formula check
The documented formula matches both sources verbatim, including the `q = x−1` off-by-one convention that
makes the test P(X ≥ x) (over-representation) rather than P(X > x). Log-space evaluation via ln Γ is a
standard numerical technique, not a formula change.

### Edge-case semantics check
- x = 0 ⇒ empty upper sum ⇒ P = 1 — sourced (PNNL §8.2 corner case). ✓
- Degenerate population (N, M, or n ≤ 0) ⇒ P = 1 — sourced/derived. ✓
- Infeasible partial tables (i > M or n−i > N−M) ⇒ C = 0 contribute 0 — sampling-without-replacement
  constraint. ✓ This also implies the **untested-in-spec** boundary x > min(M,n) ⇒ P = 0 (whole upper
  tail infeasible), which I cross-checked against SciPy and added a test for (see Stage B).

### Independent cross-check (numbers, this session, via SciPy `scipy.stats.hypergeom`)
`P(X ≥ x) = hypergeom.sf(x−1, M=N, n=M, N=n)`:

| Case | N | M | n | x | SciPy P(X≥x) | Exact | Repo expects |
|------|---|---|---|---|--------------|-------|--------------|
| PNNL | 8000 | 400 | 100 | 20 | 7.884747217e-8 | — | 7.884747232900224e-8 (Δ=1.6e-16, within 1e-15) |
| swap M↔n | 8000 | 100 | 400 | 20 | 7.884747217e-8 | symmetry | same ✓ |
| all-in | 10 | 5 | 5 | 5 | 0.0039682539683 | 1/252 | 1/252 ✓ |
| partial | 4 | 2 | 2 | 1 | 0.8333333333 | 5/6 | 5/6 ✓ |
| ≥1 | 10 | 5 | 5 | 1 | 0.996031746 | 251/252 | 251/252 ✓ |
| end-to-end | 10 | 5 | 3 | 3 | 0.083333333 | 1/12 | 1/12 ✓ |
| infeasible | 10 | 3 | 5 | 4 | 0.0 | 0 | 0.0 (added) |

All expected values trace to an external reference (SciPy / exact rational), not to the implementation.

### Findings / divergences (Stage A)
- PASS-WITH-NOTES for the **scope-hint mismatch** above (repo's ORA definition is internally consistent
  and authoritatively sourced; the prompt's completeness hints describe a different algorithm).
- INV-4 (M↔n symmetry) independently verified on a non-trivial asymmetric-looking case (N=10, M=3↔7,
  n=7↔3, x=2 → 0.81667 both ways).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs`:
`FindPathwayEnrichment` (873–916), `HypergeometricUpperTail` (929–947), `LogChoose` (950–956),
`LogGamma` Lanczos (959–987), `PathwayEnrichment` record (100–106).

### Formula realised correctly? (evidence)
Yes. `HypergeometricUpperTail` sums the upper tail directly `Σ_{i=x}^{min(n,M)} C(M,i)C(N−M,n−i)/C(N,n)`
in log-space (avoiding catastrophic cancellation of the `1 − lower-tail` form for tiny tails — correct
choice for the 7.88e-8 example). `LogChoose` returns −∞ for k∉[0,n], correctly zeroing infeasible
partial tables. Degenerate guard `x≤0 || M≤0 || n≤0 || N≤0 ⇒ 1.0` matches INV-2/INV-3. Result clamped
to [0,1] (INV-1). `FindPathwayEnrichment` de-duplicates query/members (set semantics), intersects
pathway with background (M relative to universe), unions query into background, sorts ascending by
p-value (INV-5), and validates nulls.

### Cross-verification table recomputed vs code
All seven rows above reproduced by the live tests (M1–M8, plus the new infeasible case), each matching
the SciPy/exact reference within the asserted tolerance. Full suite green.

### Variant/delegate consistency
`FindPathwayEnrichment` delegates to `HypergeometricUpperTail`; M8 end-to-end p (1/12) equals the direct
core call `HypergeometricUpperTail(3,10,5,3)` — consistent. No `*Fast`/instance variants exist.

### Test quality audit (gate result)
- **Sourced, not code-echoes:** every numeric expectation traces to SciPy or exact rational arithmetic
  retrieved this session, not to implementation output. PNNL value, 1/252, 5/6, 251/252, 1/12 all exact.
- **No green-washing:** exact `Is.EqualTo(...).Within(small)` used where exact values are known; the one
  comparison assertion (M9 ordering) tests the *ordering invariant* itself with the two exact endpoints
  pinned (1/12 and 1.0) — appropriate, not a weakened equality.
- **Coverage gap found & fixed:** the existing 17 tests covered M1–M14, S1–S2, C1 but **omitted the
  x > min(M,n) infeasible-overlap boundary** (P=0), a Stage-A edge case derivable from the
  sampling-without-replacement constraint and confirmed P=0 by SciPy. Added
  `HypergeometricUpperTail_OverlapExceedsAvailableSuccesses_EqualsZero` (both M- and n-limited orientations).
  This is a test-coverage gap, **not a code defect** — the code already returned 0 correctly.
- **Honest green:** FULL unfiltered suite = **6557 passed, 0 failed**, 1 skipped (unrelated MFE benchmark);
  changed test file builds warning-free (the 4 build warnings are pre-existing NUnit2007 in
  `ApproximateMatcher_EditDistance_Tests.cs`, untouched here).

### Findings / defects (Stage B)
- No code defect. One **test coverage gap** (infeasible-overlap boundary) found and fixed in-session.
- PASS-WITH-NOTES, attributable to (a) the closed coverage gap and (b) the documented background-defaulting
  assumption (union of pathway members ∪ query when no background supplied) — a convenience not prescribed
  by the sources, already disclosed in the doc and tested in S1; not correctness-affecting when a background
  is supplied.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (scope-hint mismatch only; ORA description correct & sourced).
- **Stage B:** PASS-WITH-NOTES (code correct; added 1 missing boundary test).
- **End-state:** ✅ CLEAN — no code defect; coverage gap fully closed; build 0 errors, full suite 0 failed.
- **Test-quality gate:** PASS (sourced expectations, no green-washing, all branches/edges now covered,
  honest full-suite green).
