# Validation Report: POP-FST-001 — Fixation index F-statistics (FST, FIS, FIT)

- **Validated:** 2026-06-24   **Area:** Population Genetics
- **Canonical method(s):** `PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2)`, `CalculateFStatistics(pop1Name, pop2Name, data)`, `CalculatePairwiseFst(populations)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

---

## Estimator confirmed

This is **Wright (1965) variance-based FST**, computed directly from known per-locus allele
frequencies — explicitly **NOT** Weir & Cockerham (1984) θ (no ANOVA a/b/c variance components,
no finite-sample bias correction) and **NOT** Nei's GST. W&C (1984) is cited in the Evidence doc
only as a *contrast* reference, not as the implemented method. The spec, evidence doc, source
XML-doc and test-class doc all agree; no overclaim — the code computes exactly what it advertises.

- `CalculateFst` → Wright variance FST: `FST = σ²_S / (p̄(1−p̄))` per locus, multi-locus = ratio of sums.
- `CalculateFStatistics` → heterozygosity-based F-statistics: `FIS = 1 − HI/HS`, `FIT = 1 − HI/HT`, `FST = 1 − HS/HT`.

## Stage A — Description

**Sources opened (Evidence doc §1, re-confirmed):**
- Wikipedia *Fixation index* — `F_ST = σ_S² / σ_T² = σ_S² / (p̄(1−p̄))`; equivalent form
  `(p̄(1−p̄) − Σ c_i p_i(1−p_i)) / (p̄(1−p̄))`; range [0,1]; "zero implies complete panmixia",
  "one implies all variation explained by structure / complete differentiation".
- Wikipedia *F-statistics* — partition identity `(1−F_IS)(1−F_ST) = 1−F_IT`, and the
  heterozygosity definitions `HI` (mean observed het within pops), `HS` (mean expected het within
  subpops), `HT` (expected het of pooled total).

**Formula check.** The two Wikipedia forms are algebraically identical (since
`p̄(1−p̄) − mean p_i(1−p_i) = σ²_S`), so `σ²_S/(p̄(1−p̄))` equals the classic `(H_T − H_S)/H_T`.
The code's `het = p̄(1−p̄)` and `variance = Σ c_i (p_i−p̄)²` with `c_i = n_i/N` realise the first
form exactly. Confirmed equivalent to Wright (1965).

**Edge-case semantics (all sourced).** FST ∈ [0,1]; FST=0 identical subpops; FST=1 fixed-different
alleles; denominator p̄(1−p̄)=0 (both fixed same allele / monomorphic) → return 0 (0/0 undefined,
design contract); empty input → 0; FIS may be negative (excess heterozygosity).

**Independent hand cross-checks (exact, recomputed this session with Python `fractions`):**
- p1=0.5, p2=0.5 → σ²_S=0 → **FST=0** ✓
- p1=1.0, p2=0.0 (equal n) → p̄=0.5, σ²_S=0.25, het=0.25 → **FST=1** ✓
- p1=0.8, p2=0.2 → **FST=9/25=0.36** ✓
- [0.9,0.8] vs [0.1,0.2] (multilocus) → **FST=1/2=0.50** ✓
- 0.5 vs 0.9 → **4/21**; 0.5 vs 0.6 → **1/99** ✓ (pairwise cells)

Stage A findings: none. Description is biologically and mathematically correct and faithfully sourced.

## Stage B — Implementation

**Code path:** `src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs`
- `CalculateFst` (609–646): per matched locus `pBar=(n1 p1+n2 p2)/(n1+n2)`,
  `variance=(n1(p1−pBar)²+n2(p2−pBar)²)/(n1+n2)`, `het=pBar(1−pBar)`; accumulates
  `numerator+=variance`, `denominator+=het`; returns `denominator>0 ? num/den : 0`. Multi-locus =
  ratio-of-sums (correct). **Mismatched-locus guard present** (619–623): throws `ArgumentException`
  when `pop1.Count != pop2.Count`, after the empty short-circuit (`Count==0` → 0).
- `CalculatePairwiseFst` (651–669): fills upper+lower triangle from `CalculateFst`, diagonal left 0
  → symmetric with zero diagonal by construction; inherits the mismatch guard.
- `CalculateFStatistics` (674–713): accumulates HetObs, expected het `2p(1−p)n`, pooled het
  `2 p̄(1−p̄)(n1+n2)`; `hi/hs/ht` over total N; `fis=1−hi/hs`, `fit=1−hi/ht`, `fst=1−hs/ht`, each
  guarded by `>0 ? … : 0`. Single combined-tuple list → no two-collection length-mismatch surface.

**Formula realised correctly:** yes — matches the validated Wright variance FST and
heterozygosity-based F-statistics exactly.

**Mismatched-locus-count throw — VERIFIED.** The fix from commit `994e91a7` is present in source
(guard at 619–623) and locked by tests: `CalculateFst_MismatchedLocusCounts_Throws` (pop1 longer,
asserts message contains "match") and `CalculateFst_MismatchedLocusCounts_Pop2Longer_Throws`
(symmetric). `CalculateFst_EqualLength_RegressionExactValues` locks 0/1/0.16/0.36 so equal-length
behaviour is unchanged. All three pass.

**Cross-verification (exact fractions, re-derived independently this session — not by running code):**

| Test | Inputs | Expected | Status |
|------|--------|----------|--------|
| SingleLocus | 0.8 vs 0.2 | 9/25 = 0.36 | ✓ |
| MultiLocus | [0.9,0.8] vs [0.1,0.2] | 1/2 = 0.50 | ✓ |
| Fixed differences | 1.0 vs 0.0 | 1.0 | ✓ |
| Pairwise cells | 0.5/0.9, 0.5/0.6 | 4/21, 1/99 | ✓ |
| FStatistics components | spec data | Fis 1/19, Fit 1/13, Fst 1/39 | ✓ (prior, unchanged) |
| FStatistics excess het | (60,100,80,100,0.3,0.7) | Fis −2/3, Fit −2/5, Fst 4/25 | ✓ (prior) |

**Edge cases in code:** identical → 0 ✓; fixed-different → 1 ✓; both-fixed-same / monomorphic →
denominator 0 → 0 ✓; empty → 0 ✓; mismatched locus counts → throws ✓; negative FIS for excess het ✓.
Partition identity `(1−FIT)=(1−FIS)(1−FST)` holds exactly (reduces to `(HI/HS)(HS/HT)=HI/HT`).

**Numerical robustness:** all divisions guarded by `>0`; exact-fraction test inputs avoid IEEE-754
drift; with valid frequencies σ²_S ≤ p̄(1−p̄) so the ratio stays in [0,1].

**Test quality:** 28 tests, exact sourced values (not "no-throw" tautologies), deterministic,
covering every Stage-A edge case plus the mismatched-locus contract.

## Verdict & follow-ups

- Stage A PASS, Stage B PASS. No defects. Estimator correctly identified (Wright variance FST).
- Mismatched-locus-count throw verified present in source and locked by two real tests.
- Build: succeeded (0 warnings). `~PopulationGeneticsAnalyzer_FStatistics_Tests` = 28/28 passed.
- **End state: CLEAN.** No code or test changes required this session.
