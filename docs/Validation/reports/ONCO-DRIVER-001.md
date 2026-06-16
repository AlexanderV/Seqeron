# Validation Report: ONCO-DRIVER-001 — Driver Mutation Detection (20/20 rule)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.ClassifyGene`, `ScoreDriverPotential`, `MatchCancerHotspots`, `IdentifyDriverMutations`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent of the repo)

| Source | URL | What it confirmed |
|--------|-----|-------------------|
| Vogelstein et al. 2013, *Cancer Genome Landscapes*, Science (PMC full text) | https://pmc.ncbi.nlm.nih.gov/articles/PMC3749880/ | Section "A Ratiometric Method to Identify and Classify Mut-Driver Genes". Verbatim: oncogene = "**>20%** of the recorded mutations in the gene are at recurrent positions and are missense"; TSG = "**>20%** of the recorded mutations in the gene are inactivating". Both use the strict **`>20%`** symbol. |
| Tokheim & Karchin 2020 (20/20+), Bioinformatics (PMC) | https://pmc.ncbi.nlm.nih.gov/articles/PMC7703750/ | Verbatim restatement: "OGs have **>20%** mutations causing missense changes at recurrent positions and TSGs have **>20%** mutations causing inactivating changes." Inactivating = "protein-truncating mutations (i.e. nonsense and frame-shifting mutations)". |
| Schroeder et al. 2014 (OncodriveROLE), Bioinformatics | https://academic.oup.com/bioinformatics/article/30/17/i549/201062 | Rule restated; **truncating types list (verbatim):** "mutations causing a frameshift, a gained or lost stop codon as well as mutations in **splice donor or acceptor sites**." (Writes "≥20%" for TSGs — see divergence note.) |
| Miller et al. 2017, Oncotarget | https://www.oncotarget.com/article/15514/text/ | Hotspot/recurrent position = "**at least two** mutations of the same class … at the same position" → ≥2. IDH1 hotspot = codon **R132H** (132). Truncations = "nonsense mutations and frameshift insertions and deletions." |

### Formula check
- **Oncogene** `f_OG = (# missense at recurrent positions) / N`, threshold strict `> 0.20` — matches Vogelstein verbatim ">20%".
- **TSG** `f_TSG = (# truncating) / N`, threshold strict `> 0.20` — matches Vogelstein verbatim ">20%".
- **Denominator convention:** Vogelstein's "% of the recorded mutations in the gene" → denominator = ALL recorded mutations (not just missense). The implementation uses N = total mutations for both fractions. ✅ matches the primary source. (Miller's looser phrasing "20% of all *missense* mutations" is a secondary restatement; the primary Vogelstein convention governs and is the one implemented.)
- **Truncating set** = nonsense + frameshift + splice donor/acceptor — matches OncodriveROLE's explicit list (splice sites included). ✅
- **Recurrent position** = ≥2 missense at the same protein position — matches Miller 2017. ✅

### Edge-case semantics
- Empty mutations → Ambiguous, both fractions 0 (no evidence). Sourced as conservative behavior (OncodriveROLE: low-recurrence genes left unclassified).
- Exactly 0.20 is NOT a driver (strict `>` from "more than 20%"). ✅
- Single missense at a position is not recurrent (Miller ≥2). ✅
- Dual-pass (both criteria > 0.20): source is **silent**; the unit resolves by dominant fraction, exact tie → Ambiguous. This is a documented Assumption (#1), not a sourced fact, and is correctly flagged as such.

### Independent cross-check (hand-computed, this session)
| Case | f_OG | f_TSG | Expected role |
|------|------|-------|---------------|
| IDH1 10 missense @132 | 1.00 | 0.00 | Oncogene |
| Dispersed TSG (5 nonsense+2 fs+1 splice+2 missense) | 0.00 | 0.80 | TumorSuppressor |
| 2 trunc + 8 distinct missense | 0.00 | 0.20 | Ambiguous (strict >) |
| 3 trunc + 7 distinct missense | 0.00 | 0.30 | TumorSuppressor |
| KRAS 2@codon12 + 8 distinct missense | 0.20 | 0.00 | not Oncogene (strict >) |

All match the TestSpec MUST/SHOULD expectations and trace to Vogelstein's strict `>20%` rule.

### Findings / divergences
- **D1 (minor, documented):** OncodriveROLE writes "≥20%" for TSGs; Vogelstein/Tokheim (primary) write ">20%". The unit follows the primary source (strict `>`). Correct choice; documented in Assumption #2.
- **D2 (doc-vs-doc wording, non-blocking):** INV-02 in the algorithm doc says dual-pass Oncogene needs "f_OG ≥ f_TSG", but the code (and TestSpec Assumption #1) require strict `f_OG > f_TSG` with exact tie → Ambiguous. The authoritative documented behavior (Assumption #1, code) is self-consistent; only the INV-02 "≥" phrasing is loose. No code/behavior impact. Noted, not fixed (does not affect any sourced expected value).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`
- `ClassifyGene` (738–783): fractions over N, strict `>` threshold, dual-pass dominant-fraction tie-break.
- `IsTruncating` (864–867): Nonsense | Frameshift | SpliceSite. ✅ matches OncodriveROLE.
- `CountRecurrentMissense` (873–895): per-position missense count, contributes full count when ≥ `RecurrentPositionMinCount` (=2). ✅ matches Miller.
- `ScoreDriverPotential` (794–798): `max(f_TSG, f_OG)` ∈ [0,1]. Transparent 20/20 signal (Assumption #3). ✅
- `MatchCancerHotspots` (810–816) / `IdentifyDriverMutations` (830–859): set-membership lookup; result is input subset in input order. ✅ INV-01.
- Constants: `DriverGeneFractionThreshold = 0.20`, `RecurrentPositionMinCount = 2`. ✅

### Cross-verification vs code
The full suite recomputes every Stage-A value; all match (see hand-computed table above). Confirmed by running the tests.

### Variant/delegate consistency
`ScoreDriverPotential` is defined as `max` over the same fractions `ClassifyGene` produces — consistent by construction.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** All MUST/SHOULD assertions use exact `Is.EqualTo(...).Within(1e-10)` on fractions and exact enum roles traceable to Vogelstein's strict `>20%`. M3/S3 boundary tests would FAIL a `≥`-implementation (they pin 0.20 → not driver). Good.
- **No green-washing:** No weakened assertions, widened tolerances, or skipped tests introduced; all exact values are sourced.
- **Coverage gaps found and FIXED this session (3 new tests):**
  1. `ClassifyGene_DualPassOncogeneDominant_ClassifiesOncogene` — exercised the previously-untested dual-pass dominant-fraction branch (f_OG=0.60 vs f_TSG=0.30 → Oncogene).
  2. `ClassifyGene_DualPassExactTie_ClassifiesAmbiguous` — exercised the exact-tie branch (f_OG=f_TSG=0.30 → Ambiguous).
  3. `ClassifyGene_OtherConsequence_CountsInDenominatorOnly` — exercised `MutationConsequence.Other` (previously never used in any test); confirms Other dilutes both fractions (denominator only), per Vogelstein "% of all recorded mutations".
  - Also strengthened S3 to additionally assert final role `Ambiguous` (parity with M3).
- **Honest green:** FULL unfiltered suite `dotnet test` = **Failed: 0, Passed: 6624** (was 6621; +3). `dotnet build` 0 errors; no new warnings in the touched test file.

### Findings / defects
- No implementation defect. The code faithfully realizes the validated 20/20 rule (strict `>20%`, correct truncating set incl. splice sites, ≥2 recurrence, all-mutations denominator).
- Two Stage-B test-quality coverage gaps (dual-pass branches; `Other` consequence) were the only issues; both fixed this session with sourced exact-value tests.

## Verdict & follow-ups
- **Stage A:** PASS (rule, thresholds, truncating set, recurrence, IDH1 all confirmed verbatim against Vogelstein 2013 primary text + 3 corroborating papers).
- **Stage B:** PASS (code matches; test gaps closed).
- **End-state:** ✅ CLEAN — no code defect; the only gaps (untested branches) were completely fixed, full suite green.
- **Test-quality gate:** PASS (after adding 3 branch-coverage tests + strengthening S3).
- **Non-blocking note:** INV-02 doc wording "f_OG ≥ f_TSG" vs code's strict `>` for dual-pass Oncogene (tie → Ambiguous). Cosmetic doc inconsistency; behavior is correct and self-consistent with Assumption #1.
