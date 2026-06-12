# Validation Report: POP-FST-001 — Fixation index F-statistics (FST, FIS, FIT)

- **Validated:** 2026-06-12   **Area:** Population Genetics
- **Canonical method(s):** `PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2)`, `CalculateFStatistics(pop1Name, pop2Name, data)`, `CalculatePairwiseFst(populations)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

---

## Estimator claimed and confirmed

The spec, evidence doc, source XML-doc, and test-class doc all **claim Wright (1965) variance-based FST**, and explicitly state it is **NOT** the Weir & Cockerham (1984) θ estimator (no ANOVA a/b/c variance components, no finite-sample bias correction). The checklist cites W&C (1984) only as a *contrast* reference, not as the implemented method. **No overclaim** — code computes exactly what it advertises.

- `CalculateFst` → Wright variance FST: `FST = σ²_S / (p̄(1−p̄))`.
- `CalculateFStatistics` → heterozygosity-based F-statistics: `FIS = 1 − HI/HS`, `FIT = 1 − HI/HT`, `FST = 1 − HS/HT`.

## Stage A — Description

**Sources opened:**
- Wikipedia *Fixation index* — confirms `F_ST = σ_S² / σ_T² = σ_S² / (p̄(1−p̄))`; alternative form `(p̄(1−p̄) − Σc_i p_i(1−p_i)) / (p̄(1−p̄))`; range [0,1]; "zero implies complete panmixia", "one implies all variation explained by structure / complete differentiation".
- Wikipedia *F-statistics* — confirms the partition identity `(1−F_IS)(1−F_ST) = 1−F_IT`. (Current revision does not spell out the HI/HS/HT closed forms.)
- Standard pop-gen sources (Hartl & Clark; Holsinger & Weir 2009; "Heterozygosity" overview) — confirm `HI` = mean observed het within pops, `HS` = mean expected het within subpops (random mating within each), `HT` = expected het of the pooled total; and `FIS = 1−HI/HS`, `FIT = 1−HI/HT`, `FST = 1−HS/HT`.

**Formula check.** The two Wikipedia forms are algebraically identical, and the code's `het = p̄(1−p̄)` (= HT/2) with `variance = Σ c_i (p_i−p̄)²` (= σ²_S) is exactly the first form. The identity `p̄(1−p̄) − mean p_i(1−p_i) = σ²_S` makes `σ²_S/(p̄(1−p̄))` equal to the classic `(H_T − H_S)/H_T`. Confirmed equivalent.

**Edge-case semantics (all sourced).** FST ∈ [0,1]; FST=0 for identical subpops; FST=1 for fixed-different alleles; denominator p̄(1−p̄)=0 (both fixed for same allele / monomorphic) → return 0 (0/0 undefined, design contract); empty input → 0; FIS may be negative (excess heterozygosity).

**Independent cross-checks (hand-computed, exact):**
- p1=0.5, p2=0.5 → σ²_S=0 → FST=0. ✓
- p1=1.0, p2=0.0 (equal n) → p̄=0.5, σ²_S=0.25, het=0.25 → FST=1.0. ✓
- p1=0.7, p2=0.3 (equal n) → p̄=0.5, σ²_S=0.04, het=0.25 → FST=0.16 (= (0.5−0.42)/0.5). ✓
- p1=0.8, p2=0.2 → FST=0.09/0.25=0.36. ✓

Stage A findings: none. Description is biologically and mathematically correct and faithfully sourced.

## Stage B — Implementation

**Code path:** `src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs`
- `CalculateFst` (lines 583–614): per matched locus `pBar = (n1 p1 + n2 p2)/(n1+n2)`, `variance = (n1(p1−pBar)² + n2(p2−pBar)²)/(n1+n2)`, `het = pBar(1−pBar)`; accumulates `numerator += variance`, `denominator += het`; returns `denominator>0 ? num/den : 0`. Multi-locus = ratio-of-sums (correct).
- `CalculatePairwiseFst` (619–637): fills upper+lower triangle from `CalculateFst`, diagonal left 0 → symmetric with zero diagonal by construction.
- `CalculateFStatistics` (642–681): accumulates `HetObs`, expected het `2p(1−p)n` per pop, pooled het `2 p̄(1−p̄)(n1+n2)`; `hi/hs/ht` divided by total N; `fis=1−hi/hs`, `fit=1−hi/ht`, `fst=1−hs/ht`, each guarded by `>0 ? … : 0`.

**Formula realised correctly:** yes — matches the validated Wright variance FST and heterozygosity-based F-statistics exactly.

**Cross-verification table recomputed vs code (exact fractions, verified independently with Python `fractions`):**

| Test | Inputs | Expected | Status |
|------|--------|----------|--------|
| SingleLocus | 0.8 vs 0.2 | 0.36 | ✓ |
| MultiLocus | [0.9,0.8] vs [0.1,0.2] | 0.50 | ✓ |
| Fixed differences | 1.0 vs 0.0 | 1.0 | ✓ |
| UnequalSampleSizes | (0.5,1000) vs (0.9,10) | 0.0062743… ; equal-size = 4/21 | ✓ |
| MultiLocusModerate | binary fracs | 1/19 | ✓ |
| WrightScale | little / very-great | 1/2499 ; 61/198 | ✓ |
| IslandModel | 0.48/0.35/0.15 vs 0.5 | 1/2499, 9/391, 49/351 | ✓ |
| Pairwise cells | 0.5/0.6/0.9 | 1/99, 4/21, 3/25 | ✓ |
| FStatistics ReturnsAllComponents | spec data | Fis 1/19, Fit 1/13, Fst 1/39 | ✓ |
| FStatistics Excess het | (60,100,80,100,0.3,0.7) | Fis −2/3, Fit −2/5, Fst 4/25 | ✓ |

All F-statistics fractions re-derived independently (not by running the code) and matched.

**Edge cases in code:** identical subpops → 0 ✓; fixed-different → 1 ✓; both-fixed-same / monomorphic → denominator 0 → 0 ✓; empty → 0 ✓; single locus valid ✓; negative FIS for excess het ✓. Partition identity `(1−FIT)=(1−FIS)(1−FST)` holds exactly because it reduces to `(HI/HS)(HS/HT)=HI/HT` (tested twice).

**Numerical robustness:** all divisions guarded by `>0` checks; no overflow on stated ranges; exact-fraction test inputs avoid IEEE-754 drift.

**Test quality:** 25 tests, exact sourced values (not "no-throw"), deterministic, covering every Stage-A edge case.

**Minor observations (non-defects):**
1. `CalculateFst` iterates `Math.Min(pop1.Count, pop2.Count)` — if the two populations are passed unequal locus counts, surplus trailing loci are silently dropped rather than throwing. Harmless for the documented matched-locus contract; no test exercises mismatched lengths.
2. FST result is not explicitly clamped to [0,1]; with valid allele frequencies σ²_S ≤ p̄(1−p̄) so the ratio is mathematically in [0,1] and all tests confirm it. No negative-FST clamping is needed because Wright's parametric (not bias-corrected) estimator cannot go negative here.

## Verdict & follow-ups

- Stage A PASS, Stage B PASS. No defects. Estimator correctly identified and not overclaimed (Wright variance FST, with W&C cited only as contrast).
- Build: succeeded (0 warnings). Tests: `~FStatistics` 25/25 passed; full suite 4484/4484 passed (baseline matched).
- **End state: CLEAN.** No code or test changes required.

---

## Fix applied (2026-06-12)

Minor observation #1 (silent truncation of mismatched locus counts in `CalculateFst`) was
re-classified as a silent-data hazard for this mission-critical library and fixed.

**Defect.** `CalculateFst` looped `Math.Min(pop1.Count, pop2.Count)` over the two populations'
per-locus allele frequencies. When the two populations were passed UNEQUAL locus counts, the
surplus trailing loci of the longer population were silently dropped and a value was still
returned — masking a caller bug instead of signalling it.

**Fix.** Added a defensive contract guard in `CalculateFst`
(`src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs`) that throws
`ArgumentException` (matching the existing guard style in this file's siblings, e.g.
`PhylogeneticAnalyzer` "dimensions must match" checks) when `pop1.Count != pop2.Count`, with a
message stating the per-locus counts must match and reporting both counts. The loop now iterates
the full (equal) length. The empty-input short-circuit (`Count == 0` → 0) and all equal-length
behaviour are unchanged — equal-length results are byte-for-byte identical.

**Sibling-method check.** `CalculatePairwiseFst` delegates to `CalculateFst`, so it inherits the
guard automatically. `CalculateFStatistics` takes a single list of combined per-locus tuples
`(HetObs1, N1, HetObs2, N2, AlleleFreq1, AlleleFreq2)` — there are no two separate collections that
can disagree in length, so it has no equivalent silent-min pattern and was left unchanged. No other
`Math.Min(...Count, ...Count)` locus-truncation pattern exists in the file (remaining `Math.Min`
calls are value clamps for MAF / D′).

**Tests added** (`tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_FStatistics_Tests.cs`,
new region "Mismatched-Length Contract Tests"):
- `CalculateFst_MismatchedLocusCounts_Throws` — pop1 longer → `ArgumentException` whose message
  states the counts must match. (Confirmed FAILING first against the pre-fix code, which silently
  truncated and returned a value.)
- `CalculateFst_MismatchedLocusCounts_Pop2Longer_Throws` — symmetric case, pop2 longer → throws.
- `CalculateFst_EqualLength_RegressionExactValues` — regression lock on worked examples
  0 / 1 / 0.16 / 0.36, confirming equal-length results are unchanged.

**Verification.** `dotnet build` 0 errors; `~FStatistics` 28/28 passed; full suite 4498/4498 passed
(prior 4495 baseline + 3 new tests). Equal-length results unchanged.
