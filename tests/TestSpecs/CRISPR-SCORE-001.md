# Test Specification: CRISPR-SCORE-001

**Test Unit ID:** CRISPR-SCORE-001
**Area:** MolTools
**Algorithm:** CRISPR published scoring models (Doench 2014 Rule Set 1, Doench 2016 Rule Set 2 / Azimuth, MIT/Hsu 2013, CFD 2016)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-18

> Scope note. This unit covers the four **published, literature-grounded** predictive scoring methods on
> `CrisprDesigner` (enhancement C7). They are distinct from the repository's heuristic guide-design /
> off-target surface (`DesignGuideRnas` / `EvaluateGuideRna` → CRISPR-GUIDE-001; `FindOffTargets` /
> `CalculateSpecificityScore` → CRISPR-OFF-001), which remain unchanged. Two on-target scorers (Rule Set 1,
> Rule Set 2) and two off-target scorers (MIT/Hsu, CFD) are grouped here because they share the same evidence
> discipline: every coefficient/matrix is transcribed verbatim from an authoritative source (or, for the
> trained Rule Set 2 model, reconstructed from the published model file), and every expected value is
> re-derived independently of the C# code.

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| S1 | Doench, Hartenian, Graham, et al. "Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation." Nat Biotechnol 32:1262-1267 (2014). | Primary literature | doi:10.1038/nbt.3026 · PMID 25184501 | 2026-06-17 |
| S2 | Doench, Fusi, Sullender, Hegde, et al. "Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9." Nat Biotechnol 34:184-191 (2016). | Primary literature (Rule Set 2 **and** CFD) | doi:10.1038/nbt.3437 · PMID 26780180 | 2026-06-18 |
| S3 | Hsu, Scott, Weinstein, Ran, Konermann, et al. "DNA targeting specificity of RNA-guided Cas9 nucleases." Nat Biotechnol 31:827-832 (2013). | Primary literature | doi:10.1038/nbt.2647 · PMID 23873081 | 2026-06-17 |
| R1 | CRISPOR reference implementation `doenchScore.py` (Rule Set 1) and `crispor.py` `calcHitScore`/`calcMitGuideScore`/`hitScoreM` (MIT/Hsu). Haeussler et al. Genome Biol 17:148 (2016). | Canonical reference code | https://github.com/maximilianh/crisporWebsite | 2026-06-17 |
| R2 | Microsoft Research **Azimuth** trained Rule Set 2 models `saved_models/V3_model_nopos.pickle`, `V3_model_full.pickle` (BSD-3-Clause), and `azimuth/features/featurization.py`. | Canonical trained model + featurization | https://github.com/MicrosoftResearch/Azimuth | 2026-06-18 |
| R3 | CFD matrices `mismatch_score.pkl` (240) / `pam_scores.pkl` (16) + `cfd-score-calculator.py` (Doench lab, redistributed by CRISPOR `CFD_Scoring/`) and the independent `bm2-lab/iGWOS` (`CFD/`, `otscore.py` with doctest oracles). | Two independent canonical distributions | https://github.com/maximilianh/crisporWebsite · https://github.com/bm2-lab/iGWOS | 2026-06-17 |
| R4 | Biopython `Bio.SeqUtils.MeltingTemp.Tm_NN` (DNA_NN3, Allawi & SantaLucia 1997), used by Azimuth for Rule Set 2 melting-temperature features. | Canonical reference code | https://github.com/biopython/biopython | 2026-06-18 |

### 1.2 Key Evidence Points

1. **Rule Set 1 is a published logistic linear model** over a 30-nt context (4 up + 20 protospacer + 3 PAM + 3 down): `intercept + abs(10-gcCount(seq[4:24]))·gcWeight + Σ position-specific single/di-nucleotide weights`, then `1/(1+e^-score)` — coefficients in R1 (S1).
2. **Rule Set 2 / Azimuth is a trained gradient-boosted-tree model** (100 trees, depth 3, learning-rate 0.1, init = training-target mean 0.5023237), NOT a coefficient table; it is only reproducible from the trained model file (S2, R2).
3. **Rule Set 2 featurization** uses order-1/2 position-dependent + position-independent nucleotide features, GC features, an NGGX interaction, and **4 nearest-neighbor melting-temperature features** (R2, R4); the 627-/630-column matrix is assembled in the CPython-2.7 dict iteration order used at training (S2, R2).
4. **MIT/Hsu single-hit score** = `Π over mismatched i of (1-W[i]) × distanceTerm × 1/nmm² × 100` with the published 20-element weight vector `W`; aggregate = `100/(100+Σ hitScores)·100` (S3, R1).
5. **CFD off-target score** = `Π over the 20 protospacer positions of the per-position mismatch percent-activity × the PAM-activity score`, key `r{guide,T→U}:d{complement(offTarget)},{i+1}` (S2, R3).
6. For all four models the **30-mer / 20-mer orientation is fixed by the source**: index 0 = 5' / PAM-distal, the high-index end = 3' / PAM-proximal (seed).

### 1.3 Documented Corner Cases

- **Rule Set 1 / Rule Set 2:** require exactly 30 nt with an N`GG` PAM at offsets 25-26 (the model is defined only on that context; S1, S2).
- **Rule Set 2 short Tm segments:** melting temperature is computed on 5-/8-/5-nt sub-segments and can be strongly negative for AT-rich segments — these are valid model inputs (R2, R4).
- **CFD:** defined only for 20-nt vs 20-nt protospacers with a 2- or 3-nt PAM (only the last two PAM nt are scored); insertions/deletions and non-ACGT are undefined (S2, R3).
- **MIT/Hsu `score2`:** the inter-mismatch-distance term applies only with ≥2 mismatches; `score3 = 1/nmm²` only with ≥1 mismatch (R1).

### 1.4 Known Failure Modes / Pitfalls

1. **Position-axis reversal** (the classic CRISPR scoring bug) — silently swapping PAM-distal/PAM-proximal indexing. Guarded explicitly for MIT/Hsu and CFD (S3, S2, R1, R3).
2. **Rule Set 2 column-order drift** — azimuth concatenates feature blocks in CPython-2.7 dict order; any reordering misaligns the trained tree feature indices (S2, R2).
3. **Rule Set 2 upstream fixture is stale** — `azimuth/tests/1000guides.csv` is internally inconsistent with the shipped pickles for ~38% of rows (its own test warns "can fail due to randomness ... feature reordering"); it is therefore NOT a strict oracle (R2).
4. **CFD off-base key** — `dY` is the **complement** of the off-target base, not the off-target base; mis-keying changes the penalty (R3).
5. **Tm parameterization** — Rule Set 2's melting temperature must use Biopython `Tm_NN` with DNA_NN3, salt method 5, dnac1=dnac2=25 nM, Na=50 mM; any other parameters shift the 4 Tm features (R4).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateOnTargetDoench2014(string context30Mer)` | `CrisprDesigner` | **Canonical** | Rule Set 1; returns 0-100. |
| `CalculateOnTargetRuleSet2(string context30Mer)` | `CrisprDesigner` | **Canonical** | Rule Set 2 / Azimuth, sequence-only (V3_model_nopos); returns ~[0,1]. |
| `CalculateOnTargetRuleSet2(string context30Mer, int aminoAcidCutPosition, double percentPeptide)` | `CrisprDesigner` | **Canonical** | Rule Set 2 / Azimuth, gene-context (V3_model_full). |
| `CalculateMitHitScore(string guide20, string offTarget20)` | `CrisprDesigner` | **Canonical** | MIT/Hsu single-hit; returns 0-100. |
| `CalculateMitSpecificityScore(IEnumerable<double>)` and the genome overload | `CrisprDesigner` | **Canonical** | MIT/Hsu aggregate. |
| `CalculateCfdScore(string sgRna20, string offTarget20, string offTargetPam)` | `CrisprDesigner` | **Canonical** | CFD off-target; returns [0,1]. |
| `AzimuthRuleSet2.{Score, Predict, Featurize…, MeltingTemp}` | `AzimuthRuleSet2` | **Internal** | Tested indirectly via `CalculateOnTargetRuleSet2`. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Rule Set 1 score ∈ [0, 100]. | Yes | sigmoid × 100 (S1) |
| INV-2 | Rule Set 2 score ∈ [0, 1] for all valid guides. | Yes | trained model output range (S2); checked over all 947 oracle guides |
| INV-3 | MIT/Hsu single-hit and aggregate scores ∈ [0, 100]; perfect match / no-hit → 100. | Yes | formula special-cases (R1) |
| INV-4 | CFD score ∈ [0, 1]; perfect match + canonical `GG` PAM → exactly 1.0. | Yes | product of factors ≤ 1 (S2, R3) |
| INV-5 | All four scorers are deterministic and case-insensitive (T↔U handled where the source requires). | Yes | pure functions over fixed tables/model |
| INV-6 | Position axis is not reversed: PAM-distal (index 0) vs PAM-proximal (high index) give different, source-correct values. | Yes | orientation counterfactuals (S3, S2) |
| INV-7 | Rule Set 2 C# output reproduces the verified Python reference (≡ `azimuth.model_comparison.predict`) for both models. | Yes | 947-guide oracle, < 1e-5 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Rule Set 1 reference examples | The two 30-mers shipped in `doenchScore.py` | 71.3089368437 / 1.89838463593 (×100, tol 1e-4) | S1, R1 |
| M2 | Rule Set 1 independent value | All-A 30-mer with PAM, re-derived in Python | 4.4338168085 | S1, R1 |
| M3 | Rule Set 2 nopos == reference (all guides) | C# vs verified Python reference, 947 guides | max \|Δ\| < 1e-5 | S2, R2 |
| M4 | Rule Set 2 full == reference (all guides) | C# vs verified Python reference, 947 guides | max \|Δ\| < 1e-5 | S2, R2 |
| M5 | Rule Set 2 == upstream fixture (agreeing subset) | C# vs `1000guides.csv` on the 585 (nopos) / 637 (full) rows where the verified reference and upstream concur | all < 1e-3 | R2 |
| M6 | MIT/Hsu perfect & known single hits | perfect→100; W[0]→100; pos5(0.395)→60.5; pos19(0.583)→41.7 | exact | S3, R1 |
| M7 | MIT/Hsu two-mismatch (all terms) | mm{5,15} → 0.10406·0.345454·0.25·100 | 0.8987 | S3, R1 |
| M8 | MIT/Hsu aggregate | one hit 60.5 → 62.305; two {60.5,41.7} → 49.456; no hits → 100 | exact | S3, R1 |
| M9 | CFD perfect match + GG | 20-mer vs itself + GG | exactly 1.0 | S2, R3 |
| M10 | CFD published doctest oracles | iGWOS `calcCfdScore` doctests | 0.4635989007074176 / 0.5140384614450001 | R3 |
| M11 | CFD single-mismatch matrix entries | rG:dT,1→0.9; rC:dT,5→0.571428571; rU:dG,7→0.6875; rG:dA,16→0.0; rC:dT,20→0.5 | exact matrix values | S2, R3 |
| M12 | CFD PAM application | perfect + GA/AG/TG → 0.069444/0.259259/0.038961; + AA → 0.0 | exact | S2, R3 |
| M13 | Orientation guards | MIT pos13(0.851)→14.9 vs pos0(0)→100; CFD rC:dT,1=1.0 vs rC:dT,20=0.5 | non-reversed values | S3, S2 (INV-6) |
| M14 | Input validation (all four) | null/empty, wrong length, non-ACGT, missing N`GG` PAM (on-target) / wrong PAM length (CFD) | `ArgumentNullException` / `ArgumentException` | implementation contract |

### 4.2 SHOULD Tests

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Range invariants | Rule Set 1 ∈[0,100]; Rule Set 2 ∈[0,1]; CFD ∈[0,1] | bounds hold | INV-1/2/4 |
| S2 | Case-insensitivity | lowercase == uppercase for all four | equal | INV-5 |
| S3 | Determinism | same input twice → identical | equal | INV-5 |
| S4 | Rule Set 2 full ≠ nopos | gene context changes the score | differs | gene-position features (S2) |
| S5 | CFD mismatch × non-canonical PAM | combined product | product of both factors | S2, R3 |

### 4.3 COULD Tests

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | MIT/Hsu genome overload | per-hit scores combined over a scanned genome | matches aggregate formula | reuses `FindOffTargets` |
| C2 | CFD 2-nt vs 3-nt PAM equivalence | only last 2 PAM nt scored | identical | S2, R3 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/CrisprDesigner_Doench2014_Tests.cs` — Rule Set 1 (9 methods).
- `tests/Seqeron/Seqeron.Genomics.Tests/CrisprDesigner_RuleSet2_Tests.cs` — Rule Set 2 / Azimuth (12 methods).
- `tests/Seqeron/Seqeron.Genomics.Tests/CrisprDesigner_Cfd_Tests.cs` — CFD (30 methods).
- `tests/Seqeron/Seqeron.Genomics.Tests/CrisprDesigner_OffTarget_Tests.cs` — MIT/Hsu block (12 methods, alongside the heuristic off-target tests of CRISPR-OFF-001).
- Oracle data: `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Azimuth/{nopos,full}_oracle.csv` (embedded), generated by `scripts/azimuth/extract_azimuth_model.py`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1, M2 (Rule Set 1) | ✅ Covered | `Doench2014_ReferenceExample1/2`, `Doench2014_AllAdenineWithPam` |
| M3, M4 (Rule Set 2 reference) | ✅ Covered | `RuleSet2_{NoPos,Full}_ReproducesVerifiedReference_AllGuides` |
| M5 (Rule Set 2 upstream subset) | ✅ Covered | `RuleSet2_{NoPos,Full}_MatchesUpstreamFixture_OnAgreeingSubset` (counts locked 585/637) |
| M6, M7, M8 (MIT/Hsu) | ✅ Covered | `CalculateMitHitScore_*`, `CalculateMitSpecificityScore_*` |
| M9-M12 (CFD) | ✅ Covered | `Cfd_PerfectMatch_*`, `Cfd_PublishedReferenceOracle_*`, `Cfd_SingleMismatch_*`, `Cfd_*Pam*` |
| M13 (orientation guards) | ✅ Covered | `CalculateMitHitScore_WeightOrientation_*`, `Cfd_OrientationGuard_*` |
| M14 (validation) | ✅ Covered | per-fixture null/empty/length/non-ACGT/PAM tests |
| S1-S4 (range/case/determinism/full≠nopos) | ✅ Covered | `*_ScoreRange*`, `*_IsInUnitInterval`, `*_LowerCase*`/`*CaseInsensitive`, `*Deterministic*`, `RuleSet2_FullModel_AddsGeneContextSignal_OverNoPos` |
| S5, C1, C2 | ✅ Covered | `Cfd_MismatchTimesNonCanonicalPam_*`, `CalculateMitSpecificityScore_Genome_*`, `Cfd_PerfectMatch_ThreeNtPam_*` |

### 5.3 Consolidation Plan

- **Canonical files:** keep one fixture per model (`CrisprDesigner_Doench2014_Tests`, `CrisprDesigner_RuleSet2_Tests`, `CrisprDesigner_Cfd_Tests`, and the MIT block inside `CrisprDesigner_OffTarget_Tests`). No merge — the per-model split mirrors the distinct sources and keeps fixtures focused.
- **Remove:** nothing.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `CrisprDesigner_Doench2014_Tests.cs` | Rule Set 1 | 9 |
| `CrisprDesigner_RuleSet2_Tests.cs` | Rule Set 2 / Azimuth | 12 |
| `CrisprDesigner_Cfd_Tests.cs` | CFD | 30 |
| `CrisprDesigner_OffTarget_Tests.cs` (MIT block) | MIT/Hsu | 12 |

### 5.5 Phase 7 Work Queue

No ❌ or ⚠ rows in §5.2 — all MUST/SHOULD/COULD cases are implemented.

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| — | (none) | — | — | — |

**Total items:** 0 · **✅ Done:** 0 · **⛔ Blocked:** 0 · **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1-M14, S1-S5, C1-C2 | ✅ | All covered; full suite **6825 passed, 0 failed**; build 0 errors, 0 warnings. |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | The verified Python reference equals `azimuth.model_comparison.predict` for the shipped Rule Set 2 pickles. **Mitigated, not blind:** validated three independent ways — extracted trees reproduced by scikit-learn 1.6's own `Tree.predict` bit-for-bit, featurizer 1e-13-identical to a verbatim port of upstream `featurization.py` with real Biopython, and column order matching documented CPython-2.7 dict iteration. | M3, M4 (Rule Set 2 oracle) |

---

## 7. Open Questions / Decisions

1. **Rule Set 2 upstream fixture is not a strict oracle** (resolved): its ~38% drift from the shipped pickles is documented; the authoritative oracle is the verified reference, with the fixture used only as independent corroboration on the agreeing subset.
2. None outstanding.
