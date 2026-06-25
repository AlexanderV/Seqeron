# Validation Report: ONCO-PURITY-001 — Tumor Purity Estimation

- **Validated:** 2026-06-24   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.EstimatePurityFromVaf(double)`,
  `OncologyAnalyzer.EstimatePurityFromVAF(IEnumerable<VariantObservation>)`,
  `OncologyAnalyzer.EstimatePurity(IEnumerable<PurityVariant>)`
  (private helper `EstimatePurityFromAlleleSpecificVaf`)
  — `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:451-575`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

## Scope vs the limitations campaign (commit 940995dd)

The prompt flags this unit as "changed" by commit `940995dd` (purity derived upstream from
allele-specific segments, ASCAT-style, with multiplicity). I verified this is **scope-only**, not
a change to ONCO-PURITY-001's own contract:

- `git log -L 451,575:OncologyAnalyzer.cs` shows the three `EstimatePurity*` methods and the
  `EstimatePurityFromAlleleSpecificVaf` helper were last touched by **64660c83** (the original
  ONCO-PURITY-001 implementation). Commit `940995dd` did **not** modify them.
- `940995dd` is purely **additive**: it introduced `SegmentAlleleSpecific`, `FitPurityPloidy`,
  `DeriveMultiplicity` (the ASCAT/CNAqc upstream derivation). The commit message itself states the
  existing `EstimatePurity/EstimatePloidy/EstimateCcf/ClassifyClonality` methods are "byte-for-byte
  unchanged". That upstream derivation is validated separately as **ONCO-ASCAT-001** (report:
  `docs/Validation/reports/ONCO-ASCAT-001.md`, CLEAN).

This report therefore (re)validates ONCO-PURITY-001's **own** method and contract: the closed-form
VAF→purity estimator and its allele-specific generalisation — not ASCAT again. All Stage-A sources
were re-opened live this session; all cross-checks recomputed by hand against the current code.

---

## Stage A — Description

### Sources opened (live, this session)

| Source | What it confirmed |
|--------|-------------------|
| CNAqc vignette (caravagnalab.github.io/CNAqc), WebFetch | Expected-VAF of a clonal mutation, **verbatim**: `v_m(c) = mπc / [2(1−π) + π(n_A+n_B)]`; m = mutation copies (multiplicity), π = purity, c = cancer-cell fraction, n_A/n_B = allele-specific copy numbers; (n_A+n_B) = total copy number. |
| Antonello et al. 2024, *Genome Biology* 25:38 (CNAqc paper, ref 1) | The general clonal/subclonal form reduces, for a single clonal state (c=1), to the vignette relation; diploid-het / triploid worked figure at π=1. |
| All-FIT (Bioinformatics 2020) + IINCaSe (PMC6887630) + DeCiFer (bioRxiv) via WebSearch | Independent confirmation: "tumor purity = heterozygous mutation VAF × 2" and "CCF is twice the VAF" specifically for **diploid copy-number-neutral** loci; the simple 2× relation breaks under CNA/LOH (matches the unit's documented scope/pitfalls). |
| Carter et al. 2012, ABSOLUTE (Nat Biotechnol 30:413) | Expected allelic fraction `F = α·s_q / [α·q + 2(1−α)]` (α=purity, q=local CN, s_q=multiplicity); `2(1−α)` = contaminating-normal diploid term. Algebraically identical to CNAqc (α↔π, s_q↔m, q↔n_tot) — independent peer-reviewed model + the 2(1−π) normal term. |
| FACETS (NAR 2016, ref 4) | Mixing `m* = mΦ + (1−Φ)`: normal diploid (1,1) genotype contributes 2 copies weighted (1−Φ) — independently confirms the `2(1−π) + π·n_tot` denominator structure. |

### Formula check (cite source equation)

- **Forward model** `v = m·π / [2(1−π) + π·n_tot]` (c=1) — matches CNAqc vignette verbatim and
  ABSOLUTE `F = αs_q/[αq+2(1−α)]`.
- **Inversion** `π = 2v / [m + v(2 − n_tot)]` — hand-derived this session:
  `v[2(1−π)+π·n_tot] = mπ` → `2v − 2vπ + vπ·n_tot = mπ` → `2v = π[m + v(2 − n_tot)]` →
  `π = 2v / [m + v(2 − n_tot)]`. Exact algebraic inverse (INV-3). ✔
- **Diploid-het special case** `ρ = 2·VAF` (m=1, n_tot=2): `π = 2v/[1 + v·0] = 2v`. Confirmed by
  three independent sources (All-FIT, IINCaSe, DeCiFer). ✔

### Edge-case semantics check

- VAF=0 ⇒ purity 0; VAF=0.5 (diploid) ⇒ purity 1.0; VAF>0.5 (diploid) ⇒ ρ>1 impossible → reject
  (INV-1 + closed form).
- Empty input ⇒ purity undefined (CNAqc corner case) → `ArgumentException`; null ⇒ guard.
- Allele-specific domain: m≥1, n_tot≥1; (VAF,m,n_tot) yielding a non-positive denominator or
  ρ∉[0,1] is rejected (formula domain — the VAF is unreachable by any purity for that copy state).

### Independent cross-check (exact numbers, all from the sourced formulas)

| Case | Inputs | Hand computation | Value |
|------|--------|------------------|-------|
| Diploid het, π=60% | m=1, n=2 | v = 0.6/2 | VAF 0.30; inverse 2·0.30 = **0.60** |
| Tolerance band 55–65% | m=1, n=2 | 2·{0.275,0.325} | VAF band **0.275–0.325** |
| 2:1, m=1, π=1 | n=3 | v = 1/[0+3] | **1/3** ✔ (CNAqc "peak 33%") |
| 2:1, m=2, π=1 | n=3 | v = 2/3 | **2/3** ✔ (CNAqc "peak 66%") |
| general, π=0.5 | m=1, n=4 | v = 0.5/[1+2] | **1/6**; inverse (1/3)/(2/3) = **0.5** |
| ABSOLUTE diploid het | α=0.6, q=2, s_q=1 | F = 0.6/(1.2+0.8) | **0.30** ✔ |

### Findings / divergences (Stage A)

None. The model is biologically and mathematically correct and faithful to two independent
peer-reviewed models (CNAqc and ABSOLUTE) plus the FACETS mixing structure. The VAF-only
estimator's fixed (m=1, n_tot=2) scope is an explicitly stated, correctly-sourced modeling choice
(documented Assumption 1), not an error — the robust copy-neutral-diploid case the literature
itself singles out.

---

## Stage B — Implementation

### Code path reviewed (file:line)

- `OncologyAnalyzer.cs:480-497` `EstimatePurityFromVaf`: `ρ = HeterozygousDiploidPurityFactor(2.0)·vaf`;
  rejects NaN / vaf∉[0,1]; rejects resulting ρ>1 (VAF>0.5).
- `OncologyAnalyzer.cs:451-468` `EstimatePurityFromVAF`: per-variant VAF via `CalculateVAF`,
  ρ=2·VAF, median; empty→`ArgumentException`, null→guard.
- `OncologyAnalyzer.cs:515-575` `EstimatePurity` + `EstimatePurityFromAlleleSpecificVaf`:
  `π = 2·v / [m + v·(2 − n_tot)]` with `NormalDiploidCopyNumber = 2.0`; guards vaf∉[0,1], m<1,
  n_tot<1, denominator≤0, purity∉[0,1]; median aggregation.
- `Median` (`:578-585`): sorts a **copy**, lower-mid average for even counts. Deterministic,
  non-mutating, order-independent.

### Formula realised correctly?

Yes — verbatim. `HeterozygousDiploidPurityFactor * vaf` = `2·VAF`, and
`NormalDiploidCopyNumber * vaf / (m + vaf*(NormalDiploidCopyNumber − n_tot))` =
`2v / [m + v(2 − n_tot)]` — exactly the sourced closed form and its exact inverse.

### Cross-verification table recomputed vs code (tests run)

| Case | Input | Sourced formula value | Code/test | Match |
|------|-------|----------------------|-----------|-------|
| M1/M2/M3/M5 | VAF 0.30 / 0.50 / 0.00 / 0.275 | 0.60 / 1.0 / 0.0 / 0.55 | same | ✔ |
| M4 | VAFs {0.10,0.15,0.30} | median{0.20,0.30,0.60}=0.30 | 0.30 | ✔ |
| S1 | VAF 0.02 | 0.04 | 0.04 | ✔ |
| M6 | v0.30,m1,n2 | 0.60 (=2·VAF) | 0.60 | ✔ |
| M7 | v1/3,m1,n3 | 1.0 | 1.0 | ✔ |
| M8 | v2/3,m2,n3 | 1.0 | 1.0 | ✔ |
| M9 | v1/6,m1,n4 | (1/3)/(2/3)=0.5 | 0.5 | ✔ |
| S2 | mixed {0.60,1.0,0.5} | median=0.60 | 0.60 | ✔ |
| M10 | VAF 0.6 diploid; (0.6,1,2) | ρ=1.2>1 → throw | throws | ✔ |
| M10b | (0.9,1,4) | denom 1+0.9·(−2)=−0.8≤0 → throw | throws | ✔ |
| M11–M14 | empty/null/n0/m0/vaf∉[0,1] | throw | throws | ✔ |

### Variant/delegate consistency

`EstimatePurity` at (m=1, n_tot=2) reduces to ρ=2·VAF, agreeing with `EstimatePurityFromVaf` /
`EstimatePurityFromVAF` (M6 vs M1). Both collection overloads delegate to the single-variant
helper + `Median`. Consistent.

### Numerical robustness

Closed-form arithmetic; no overflow on stated ranges. Denominator≤0 and purity∉[0,1] explicitly
guarded; NaN VAF rejected. `Median` copies the input (no mutation; order/state-independent — C1).

### Test quality audit (HARD gate)

- **Sourced, not code echoes:** every expected value is hand-derived from the CNAqc/ABSOLUTE
  formulas above (e.g. 0.60, 1.0, 0.5, 0.55), not read off the implementation. A wrong-constant
  implementation would fail M1/M6/M7/M8/M9.
- **No green-washing:** all positive assertions use exact `Is.EqualTo(…).Within(1e-10)` (1e-12 for
  determinism). No `Greater`/`AtLeast`/range where an exact value is known; no skip/ignore.
- **Covers all logic:** all three public methods exercised; both VAF-only overloads; allele-specific
  diploid, 2:1 peaks (both multiplicities), general (m1,n4) inversion; median for both overloads;
  determinism; and error branches null/empty/VAF-out-of-range/VAF>0.5/m<1/n_tot<1/purity>1/
  non-positive-denominator (M10b).
- **Honest green:** filtered `OncologyAnalyzer_EstimatePurity_Tests` = **22 passed, 0 failed**;
  `dotnet build` 0 warnings / 0 errors (SDK 10.0.301, net10.0).

### Findings / defects (Stage B)

None. No code or description change required this session. The unchanged methods continue to
realise the validated closed form and its exact inverse.

---

## Verdict & follow-ups

- **Stage A: PASS** — model confirmed verbatim against CNAqc (vignette + paper) and ABSOLUTE,
  cross-checked against FACETS mixing and three purity-from-VAF tools; all worked numbers reproduced
  by hand.
- **Stage B: PASS** — code realises the closed form `ρ=2·VAF` and the exact inverse
  `π = 2v/[m + v(2−n_tot)]`; every cross-check value matches; tests are exact and sourced.
- **End state: CLEAN.** ONCO-PURITY-001's own contract was unchanged by commit 940995dd (that
  commit's ASCAT upstream is the separately-validated ONCO-ASCAT-001). Documented residual: the
  VAF-only estimator deliberately fixes (m=1, n_tot=2) — by-design, with the allele-specific
  `EstimatePurity` covering other copy states via explicit (m, n_tot). No code changed; no test
  added (coverage already complete, incl. the M10b non-positive-denominator branch).

CodeChanged: no. FullSuite: N/A (no code modified; the unit's 22 tests run green).
