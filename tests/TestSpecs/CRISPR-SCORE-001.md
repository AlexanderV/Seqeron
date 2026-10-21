# Test Specification: CRISPR-SCORE-001

**Test Unit ID:** CRISPR-SCORE-001
**Area:** MolTools
**Algorithm:** CRISPR published scoring models (Doench 2014 Rule Set 1, Doench 2016 Rule Set 2 / Azimuth, MIT/Hsu 2013 and its SaCas9 / Tycko 2018 sibling, CFD 2016)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-10-21 (review-2026-09 batch B08: F9-F11, F14, F15, F35/F36, F38, F50, F66/F67, F70)

> Scope note. This unit covers the four **published, literature-grounded** predictive scoring methods on
> `CrisprDesigner` (enhancement C7). They are distinct from the repository's heuristic guide-design /
> off-target surface (`DesignGuideRnas` / `EvaluateGuideRna` → CRISPR-GUIDE-001; `FindOffTargets` /
> `CalculateSpecificityScore` → CRISPR-OFF-001), which remain unchanged. Two on-target scorers (Rule Set 1,
> Rule Set 2) and two off-target scorers (MIT/Hsu, CFD) are grouped here because they share the same evidence
> discipline: every coefficient/matrix is transcribed verbatim from an authoritative source (or, for the
> trained Rule Set 2 model, reconstructed from the published model file), and every expected value is
> re-derived independently of the C# code.
>
> Update 2026-10-21 (review-2026-09 batch B08). The unit grew beyond the original four methods: the SaCas9 off-target
> model of Tycko et al. 2018 (`CalculateSaCas9HitScore`), the PAM-annotated / alternative-PAM MIT aggregates
> (`CalculateMitSpecificityScore` over `PamAnnotatedOffTarget`, `CalculateMitSpecificityScoreWithAlternativePams` with the
> `SpCas9_NGA` / `SaCas9_NNGRRN` systems and `minAlternativePamHitScore`), the CFD genome aggregates with
> `CfdAggregation`, and the Rule Set 2 argument contracts (`percentPeptide == -1` sentinel, `[0, 100]`,
> `aminoAcidCutPosition >= 0`). These are specified in §4.1 (M15-M28); CRISPR-OFF-001 owns the same code from the
> off-target side and holds the long-form numbers (`tests/TestSpecs/CRISPR-OFF-001.md`).

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
| S4 | Tycko, Wainberg, Marinov, et al. "Pairwise library screen systematically interrogates *Staphylococcus aureus* Cas9 specificity in human cells." Nat Commun 9:2962 (2018). | Primary literature (SaCas9 off-target model) | doi:10.1038/s41467-018-05391-2 | see `docs/Evidence/CRISPR-OFF-001-Evidence.md` |
| R5 | CRISPOR `crispor.py` `calcSaHitScore` / `annotateOfftargets` (alternative-PAM x0.2 weight L1947-1958, `ALTPAMMINSCORE` L318) and the SaCas9 table `bin/src/pairwise-library-screen/models/rstan_mult_21_coeffs.txt` (252 coefficients, mechanically diffed: 0 mismatching cells). | Canonical reference code | https://github.com/maximilianh/crisporWebsite | B08 F14, F50, F66, F67 (dates in `CRISPR-OFF-001-Evidence.md` / `review-2026-09/B08.md`) |
| R6 | Azimuth `azimuth/model_comparison.py` `predict` (sentinel `percent_peptide == -1` at L539/L563; `np.isreal` checks at L526-534) and the regression fixture `azimuth/tests/1000guides.csv` (column `AA cut`, integers 0..499). | Canonical reference code | https://github.com/MicrosoftResearch/Azimuth | B08 F35, F36, F70 |

### 1.2 Key Evidence Points

1. **Rule Set 1 is a published logistic linear model** over a 30-nt context (4 up + 20 protospacer + 3 PAM + 3 down): `intercept + abs(10-gcCount(seq[4:24]))·gcWeight + Σ position-specific single/di-nucleotide weights`, then `1/(1+e^-score)` — coefficients in R1 (S1).
2. **Rule Set 2 / Azimuth is a trained gradient-boosted-tree model** (100 trees, depth 3, learning-rate 0.1, init = training-target mean 0.5023237), NOT a coefficient table; it is only reproducible from the trained model file (S2, R2).
3. **Rule Set 2 featurization** uses order-1/2 position-dependent + position-independent nucleotide features, GC features, an NGGX interaction, and **4 nearest-neighbor melting-temperature features** (R2, R4); the 627-/630-column matrix is assembled in the CPython-2.7 dict iteration order used at training (S2, R2).
4. **MIT/Hsu single-hit score** = `Π over mismatched i of (1-W[i]) × distanceTerm × 1/nmm² × 100` with the published 20-element weight vector `W`; aggregate = `100/(100+Σ hitScores)·100` (S3, R1).
5. **CFD off-target score** = `Π over the 20 protospacer positions of the per-position mismatch percent-activity × the PAM-activity score`, key `r{guide,T→U}:d{complement(offTarget)},{i+1}` (S2, R3).
6. For all four models the **30-mer / 20-mer orientation is fixed by the source**: index 0 = 5' / PAM-distal, the high-index end = 3' / PAM-proximal (seed).
7. **Rule Set 2 gene-position sentinel.** In the reference, `percent_peptide == -1` both selects `V3_model_nopos` and replaces the gene-position columns by `-1`; `predict(seq, aa_cut, -1)` is therefore exactly the sequence-only score (R6, F35). The `[0, 100]` bound on `percentPeptide` and the `>= 0` bound on `aminoAcidCutPosition` are this library's contract, justified against the 947-guide fixture, not reference assertions (F36, F70).
8. **SaCas9 off-target model** = 100 x the product, over the 20 positions after the 5'-most one, of the Tycko `Median` coefficient for (guide base read as U, complement of the off-target base, position remapped from the 5' end); 252 coefficients = 21 positions x 12 mismatch types; gaps are refused, as in `calcSaHitScore` (S4, R5, F14).
9. **Which model scores a hit is chosen by system** (`annotateOfftargets`): SaCas9 PAM -> Tycko; other PAM-after-target systems -> `calcHitScore`; PAM-before-target (Cas12a-class) systems get no MIT/CFD score at all (`mitScore = cfdScore = 0`, sentinels `-1`), which this library renders as a refusal, never as a perfect 100 (F9, R5).
10. **Alternative-PAM hits** enter the MIT aggregate at x0.2 when the searched pattern is exactly `NGG` and the hit's last two PAM bases are not `GG`; the CFD score is not downweighted (`CalculateCfdSpecificityScoreWithAlternativePams`, F80, sums the raw CFD of the same hit set). CRISPOR drops alternative-PAM hits whose MIT hit score is below `ALTPAMMINSCORE = 1.0` (parameter here, default 1.0 since G6) (R5, F50, F66, F67).

### 1.3 Documented Corner Cases

- **Rule Set 1 / Rule Set 2:** require exactly 30 nt with an N`GG` PAM at offsets 25-26 (the model is defined only on that context; S1, S2). `null` is an `ArgumentNullException`; `""` is a wrong value and an `ArgumentException` reporting the length, not `ArgumentNullException` (F36, F38).
- **Rule Set 2 gene context:** `percentPeptide` must be in `[0, 100]` or exactly `-1` (sentinel -> sequence-only model, `aminoAcidCutPosition` ignored); `aminoAcidCutPosition >= 0` otherwise (F35, F36, F70). NaN and +-infinity are `ArgumentOutOfRangeException`.
- **MIT/Hsu aggregate on non-20-nt systems:** SaCas9 (21 nt) is scored with the Tycko model per hit (F14); a non-SaCas9 guide longer than 20 nt would use `calcHitScore`'s last-20 fudge and a 19-nt guide its `A`-pad, but no current system reaches those branches (F9). PAM-before-target systems throw `ArgumentException` (F9, F15).
- **CFD aggregates:** `CfdAggregation.MitStyle` = `100/(100 + Σ CFD) x 100` (in `(0, 100]`); `ParkinsonNormalized` = `100/Σ CFD` over RAW CFD values in `[0, 1]`, which is **not** bounded by 100 (raw {0.5, 0.25} -> 133.33) and throws `ArgumentOutOfRangeException` for a zero sum (F11). The genome overloads validate `maxMismatches` in `[0, 5]` like `FindOffTargets` and refuse SaCas9 and PAM-before systems (F9, F10).
- **Rule Set 2 short Tm segments:** melting temperature is computed on 5-/8-/5-nt sub-segments and can be strongly negative for AT-rich segments — these are valid model inputs (R2, R4).
- **CFD:** defined only for 20-nt vs 20-nt protospacers with a 2- or 3-nt PAM (only the last two PAM nt are scored); insertions/deletions and non-ACGT are undefined (S2, R3).
- **MIT/Hsu `score2`:** the inter-mismatch-distance term applies only with ≥2 mismatches; `score3 = 1/nmm²` only with ≥1 mismatch (R1).

### 1.4 Known Failure Modes / Pitfalls

1. **Position-axis reversal** (the classic CRISPR scoring bug) — silently swapping PAM-distal/PAM-proximal indexing. Guarded explicitly for MIT/Hsu and CFD (S3, S2, R1, R3).
2. **Rule Set 2 column-order drift** — azimuth concatenates feature blocks in CPython-2.7 dict order; any reordering misaligns the trained tree feature indices (S2, R2).
3. **Rule Set 2 upstream fixture is stale** — `azimuth/tests/1000guides.csv` is internally inconsistent with the shipped pickles for ~38% of rows (its own test warns "can fail due to randomness ... feature reordering"); it is therefore NOT a strict oracle (R2).
4. **CFD off-base key** — `dY` is the **complement** of the off-target base, not the off-target base; mis-keying changes the penalty (R3).
5. **Tm parameterization** — Rule Set 2's melting temperature must use Biopython `Tm_NN` with DNA_NN3, salt method 5, dnac1=dnac2=25 nM, Na=50 mM; any other parameters shift the 4 Tm features (R4).
6. **A false 100 on a specificity score** - an aggregate that silently drops every hit (the former non-20-nt filter) reports the guide as perfectly specific; the unsafe direction. Guarded by asserting the off-targets exist before asserting the refusal / the real value (F9).
7. **Sentinel read as a coordinate** - feeding `percentPeptide = -1` to the full model (`Percent Peptide = -1`, `<50% = 1`) answers with a score the reference returns for no input (0.608099 vs the correct 0.510655 for the oracle guide) (F35).
8. **Factor-100 trap in `ParkinsonNormalized`** - the reference rescales per-hit CFD by x100 under the same flag; over RAW [0, 1] inputs the factors cancel to `100/Σ` (F11).
9. **Double counting a site** - every `NNGRRT` site is also `NNGRRN`; the alternative-PAM scan drops hits the main PAM already matches (CRISPOR runs ONE search) (F66).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateOnTargetDoench2014(string context30Mer)` | `CrisprDesigner` | **Canonical** | Rule Set 1; returns 0-100, unrounded. `null` -> `ArgumentNullException`; `""`, wrong length, non-ACGT, non-NGG -> `ArgumentException` (F38). |
| `CalculateOnTargetRuleSet2(string context30Mer)` | `CrisprDesigner` | **Canonical** | Rule Set 2 / Azimuth, sequence-only (V3_model_nopos); returns ~[0,1]. |
| `CalculateOnTargetRuleSet2(string context30Mer, int aminoAcidCutPosition, double percentPeptide)` | `CrisprDesigner` | **Canonical** | Rule Set 2 / Azimuth, gene-context (V3_model_full). `percentPeptide == -1` -> sequence-only model (F35); `percentPeptide` outside `[0, 100]` (NaN/inf included) and `aminoAcidCutPosition < 0` (unless the sentinel) -> `ArgumentOutOfRangeException` (F36, F70). |
| `CalculateMitHitScore(string guide20, string offTarget20)` | `CrisprDesigner` | **Canonical** | MIT/Hsu single-hit; returns 0-100. |
| `CalculateSaCas9HitScore(string saGuide21, string offTarget21)` | `CrisprDesigner` | **Canonical** | SaCas9 single-hit, Tycko 2018 via CRISPOR `calcSaHitScore` (21-nt table only); returns 0-100; gaps / non-ACGT / non-21-nt refused (F14). |
| `CalculateMitSpecificityScore(IEnumerable<double>)` | `CrisprDesigner` | **Canonical** | MIT/Hsu aggregate over single-hit scores. |
| `CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType)` | `CrisprDesigner` | **Canonical** | Genome aggregate; per-hit model chosen by system (SaCas9 -> Tycko, other PAM-after systems -> `calcHitScore`); PAM-before-target systems throw `ArgumentException` (F9, F14). |
| `CalculateMitSpecificityScore(string, IEnumerable<PamAnnotatedOffTarget>, CrisprSystemType)` and `FindOffTargetsWithPam` | `CrisprDesigner` | **Canonical** | Aggregate over a PAM-annotated hit set with CRISPOR's x0.2 alternative-PAM weight (F50). |
| `CalculateMitSpecificityScoreWithAlternativePams(string, DnaSequence, int, CrisprSystemType, double minAlternativePamHitScore = 1.0)` | `CrisprDesigner` | **Canonical** | Scans the main PAM plus the alternative PAMs of `GetAlternativeOffTargetPams` that exist as systems (`SpCas9_NAG`, `SpCas9_NGA`, `SaCas9_NNGRRN`); optional `ALTPAMMINSCORE` cutoff; negative / NaN cutoff throws (F50, F66, F67). |
| `CalculateCfdScore(string sgRna20, string offTarget20, string offTargetPam)` | `CrisprDesigner` | **Canonical** | CFD off-target; returns [0,1]. |
| `CalculateCfdSpecificityScore(IEnumerable<double>, CfdAggregation)` and the genome overloads | `CrisprDesigner` | **Canonical** | CFD aggregates (`MitStyle` default, `ParkinsonNormalized`); genome overloads validate `maxMismatches` and refuse SaCas9 / PAM-before systems (F9-F11). |
| `AzimuthRuleSet2.{Score, Predict, Featurize…, MeltingTemp}` | `AzimuthRuleSet2` | **Internal** | Tested indirectly via `CalculateOnTargetRuleSet2`. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Rule Set 1 score ∈ [0, 100]. | Yes | sigmoid × 100 (S1) |
| INV-2 | Rule Set 2 score ∈ [0, 1] for all valid guides. | Yes | trained model output range (S2); checked over all 947 oracle guides |
| INV-3 | MIT/Hsu single-hit and aggregate scores ∈ [0, 100]; perfect match / no-hit → 100. A system the model does not cover (PAM-before-target) is refused, never reported as 100. | Yes | formula special-cases (R1); `EnsureMitModelApplies` (F9) |
| INV-4 | CFD score ∈ [0, 1]; perfect match + canonical `GG` PAM → exactly 1.0. | Yes | product of factors ≤ 1 (S2, R3) |
| INV-5 | All four scorers are deterministic and case-insensitive (T↔U handled where the source requires). | Yes | pure functions over fixed tables/model |
| INV-6 | Position axis is not reversed: PAM-distal (index 0) vs PAM-proximal (high index) give different, source-correct values. | Yes | orientation counterfactuals (S3, S2) |
| INV-7 | Rule Set 2 C# output reproduces the verified Python reference (≡ `azimuth.model_comparison.predict`) for both models. | Yes | 947-guide oracle, < 1e-5 |
| INV-8 | `CalculateOnTargetRuleSet2(ctx, aaCut, -1)` equals `CalculateOnTargetRuleSet2(ctx)` for every `aaCut` (including negative ones). | Yes | `model_comparison.py` L539/L563 (R6); F35, F70 |
| INV-9 | The SaCas9 hit score is 100 for a perfect match and ignores a mismatch at the 5'-most position. | Yes | `calcSaHitScore` (R5); F14 |
| INV-10 | An alternative-PAM hit contributes exactly 0.2 x its unweighted hit score to the MIT sum (NGG searches only); a hit already matched by the main PAM is never counted twice. | Yes | `annotateOfftargets` L1947-1958 (R5); F50, F66 |

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
| M14 | Input validation (all four) | null/empty, wrong length, non-ACGT, missing N`GG` PAM (on-target) / wrong PAM length (CFD) | `ArgumentNullException` (null) / `ArgumentException` (empty string and the rest) | implementation contract; empty-string typing per F36, F38 |
| M15 | Rule Set 2 `percentPeptide == -1` sentinel | 947 guides x two different `aaCut` values equal the sequence-only score and the oracle `ref_score`; oracle guide `AAAAAAAAAAAAAAAAAAATGCAGCGGGAG`, aaCut 181 → 0.510655 (pre-fix 0.608099 NOT returned), full model with a real percentage still 0.570371 | exact | R6, F35 |
| M16 | Rule Set 2 `percentPeptide` range | `-1.0001`, `100.0001`, NaN, ±∞ throw `ArgumentOutOfRangeException`; `0`, `0.18`, `50`, `100`, `-1` are accepted and land in [0, 1]; model blobs declare the 627 / 630-column shape and a wrong declared shape is `InvalidOperationException` (never `IndexOutOfRangeException`) | exception type / range | library contract (F36) |
| M17 | Rule Set 2 `aminoAcidCutPosition` range | `-1`, `-2`, `int.MinValue` with a real percentage throw `ArgumentOutOfRangeException` (`ParamName = aminoAcidCutPosition`); the same values with the sentinel equal the sequence-only score; `0`, `181`, `499`, `100000` accepted | exception type / equality | R6 (`1000guides.csv` AA cut 0..499), F70 |
| M18 | Rule Set 1 reference-executed rows and contract | 11 `TestCase` rows sampled from the 477-context run (gcCount 0/2/5/7/9/10/11/13/17, lowest/highest score, three designed contexts) match the Python reference to 1e-9; both GC coefficients are reachable at the `gcCount` 10/11 boundary; `null` → `ArgumentNullException`, `""` → `ArgumentException` that is not `ArgumentNullException` (message names the 30-nt context, `ParamName = context30Mer`); a non-NGG context throws | tol 1e-9 / exception type | S1, R1, F38 |
| M19 | SaCas9 single-hit (Tycko 2018) | perfect → 100; two reference oracles; `dT:rC` at index 1 → 54.6734106002513 vs index 20 → 10.569373315278801 (orientation); 5'-most mismatch ignored; non-21-nt / non-ACGT / gap / null refused; lower case = upper case | exact | S4, R5, F14 |
| M20 | SaCas9 genome aggregate uses the Tycko model | guide `ACGTACGTACGTACGTACGTA`, genome `ACGTACGTACGTACGTACGTCAAGAAT`, 3 mismatches → 95.91285612347565 (was a false 100); two-mismatch hit and a 5'-end-only mismatch (scored as a perfect hit) match the reference | exact | R5, F9, F14 |
| M21 | PAM-before-target systems are refused, not scored 100 | `CalculateMitSpecificityScore` genome overload throws `ArgumentException` for the Cas12a-class systems even though real off-targets exist (asserted first); 4 PAM-before systems; the CFD genome overloads refuse SaCas9 and invalid arguments | `ArgumentException` | R5, F9, F10, F15 |
| M22 | CFD genome aggregate validates `maxMismatches` | `-1`, `6`, `99` → `ArgumentOutOfRangeException` for both genome overloads; null guide / genome → `ArgumentNullException` before any system check | exception type | F10 |
| M23 | `CfdAggregation` | raw {0.5, 0.25}: `MitStyle` 99.25558312655087, `ParkinsonNormalized` 133.33333333333334; raw {1.0}: 99.00990099009901 / 100.0; `ParkinsonNormalized` with Σ = 0 → `ArgumentOutOfRangeException`; `MitStyle` of an empty / zero-sum set → 100 | exact | R3, F11 |
| M24 | Alternative-PAM x0.2 weight | a `TAG` hit contributes one fifth (59.95203836930455 vs 52.083333333333336 for an `AGG` hit); the weight is on the last two PAM bases; `SpCas9-NAG` searches never downweight; SaCas9 uses Tycko and no weight; `GetAlternativeOffTargetPams` equals the reference table; invalid inputs throw | exact | R5, F50 |
| M25 | Alternative-PAM scan covers `NGA` and `NNGRRN` | one 1-mismatch `NGA` hit (hit score 60.5) → 89.20606601248885 (was 100); SaCas9 `NNGRRN` hit enters unweighted (95.91285612347565) and an `NNGRRT` twin is counted once (92.14668839146395); the new systems carry reference metadata and the existing enum values 0..6 are unchanged | exact | R5, F66 |
| M26 | `minAlternativePamHitScore` | cutoff 0 → 59.92955631050283, 1.0 → 59.95203836930455, 40 → 62.30529595015576 (never drops the main-PAM `NGG` hit); SaCas9 cutoff also applies to main-PAM hits; negative / NaN cutoff → `ArgumentOutOfRangeException`; default 1.0 (CRISPOR's `ALTPAMMINSCORE`, G6) = the explicit-1.0 result, explicit 0 = no filtering | exact | R5, F67 |
| M27 | CFD aggregate over the alternative PAMs | NAG-only genome: main-PAM overload 100, `CalculateCfdSpecificityScoreWithAlternativePams` 99.81515711670932 (per-hit CFD 0.18518518492592592, no 0.2 factor); NGA-only 99.95042141828442; mixed genome 99.24515917350747 (default cutoff 1.0) / 99.20727233740708 (cutoff 0 or 0.3: the cutoff is on the MIT hit score) / 99.29078014212564 (40); SaCas9 and PAM-before refused | exact | R5, F80 |
| M28 | On-target interval | duplicated perfect site `G+AGG+TTTT+G+TGG+TTTTTT`: no interval -> 100 (all exact matches dropped); `[0,23)` -> MIT 50.0 / CFD 99.00990099009901; `[1,23)` (partly inside) and `[100,123)` -> 33.33333333333333 / 98.0392156862745; reverse-strand duplicate and alternative-PAM exact hits likewise (71.42857142857143); invalid bounds throw | exact | R5, F81 |

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

Counts below were taken on 2026-10-21 by counting `[Test]` / `[TestCase]` attributes (methods / expanded `TestCase` rows); no test run is claimed by this documentation update.

- `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CrisprDesigner_Doench2014_Tests.cs` - Rule Set 1 (11 methods, 23 cases).
- `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CrisprDesigner_RuleSet2_Tests.cs` - Rule Set 2 / Azimuth (22 methods, 44 cases).
- `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CrisprDesigner_Cfd_Tests.cs` - CFD single-hit and aggregates (40 methods, 45 cases).
- `tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CrisprDesigner_OffTarget_Tests.cs` - MIT/Hsu block (12 methods, alongside the heuristic off-target tests of CRISPR-OFF-001) plus the B08 regions: F9 MIT on non-20-nt systems (5), F14 SaCas9 / Tycko (6), F50 alternative-PAM weight (10), round-2 E1 `NGA` / `NNGRRN` / `minAlternativePamHitScore` (5) - 38 methods in total for the unit's MIT / SaCas9 / alternative-PAM surface.
- Oracle data: `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Azimuth/{nopos,full}_oracle.csv` (embedded), generated by `scripts/azimuth/extract_azimuth_model.py`.
- MCP-level forwarding of the same methods (`calculate_on_target_rule_set2`, `calculate_on_target_doench2014`, `crispr_mit_specificity_score` including `include_alternative_pams` / `min_alternative_pam_hit_score`, `crispr_cfd_specificity_score`) is tested in `Seqeron.Mcp.MolTools.Tests`, outside this unit's counts.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1, M2 (Rule Set 1) | ✅ Covered | `Doench2014_ReferenceExample1_MatchesPublishedValue`, `Doench2014_ReferenceExample2_MatchesPublishedValue`, `Doench2014_AllAdenineWithPam_MatchesIndependentComputation` |
| M3, M4 (Rule Set 2 reference) | ✅ Covered | `RuleSet2_{NoPos,Full}_ReproducesVerifiedReference_AllGuides` |
| M5 (Rule Set 2 upstream subset) | ✅ Covered | `RuleSet2_{NoPos,Full}_MatchesUpstreamFixture_OnAgreeingSubset` (counts locked 585/637) |
| M6, M7, M8 (MIT/Hsu) | ✅ Covered | `CalculateMitHitScore_*`, `CalculateMitSpecificityScore_*` (MIT block) |
| M9-M12 (CFD) | ✅ Covered | `Cfd_PerfectMatch_*`, `Cfd_PublishedReferenceOracle_*`, `Cfd_SingleMismatch_*`, `Cfd_*Pam*` |
| M13 (orientation guards) | ✅ Covered | `CalculateMitHitScore_WeightOrientation_PamProximalIsHighPenalty`, `Cfd_OrientationGuard_Position1VsPosition20_NotReversed` |
| M14 (validation) | ✅ Covered | per-fixture null/empty/length/non-ACGT/PAM tests; empty-string typing: `Doench2014_NullOrEmpty_Throws`, `RuleSet2_NullOrEmpty_Throws` |
| S1-S4 (range/case/determinism/full≠nopos) | ✅ Covered | `*_ScoreRange*`, `*_IsInUnitInterval`, `*_LowerCase*`/`*CaseInsensitive`, `*Deterministic*`, `RuleSet2_FullModel_AddsGeneContextSignal_OverNoPos` |
| S5, C1, C2 | ✅ Covered | `Cfd_MismatchTimesNonCanonicalPam_MultipliesBoth`, `CalculateMitSpecificityScore_Genome_SingleOffTarget_MatchesFormula`, `Cfd_PerfectMatch_ThreeNtPam_ReturnsExactlyOne` |
| M15 (RS2 sentinel) | ✅ Covered | `RuleSet2_PercentPeptideSentinel_SelectsTheNoPosModel_ForEveryOracleGuide`, `RuleSet2_PercentPeptideSentinel_ReturnsTheVerifiedNoPosScore_NotTheFullModelOne` |
| M16 (RS2 `percentPeptide` range, blob shape) | ✅ Covered | `RuleSet2_PercentPeptideOutOfRange_Throws` (7), `RuleSet2_PercentPeptideInRange_IsAccepted` (5), `RuleSet2_ModelBlob_DeclaresTheShapeTheFeaturizerAssumes` (2), `RuleSet2_ModelBlob_WrongDeclaredShape_ThrowsInvalidOperation` (2) |
| M17 (RS2 `aminoAcidCutPosition`) | ✅ Covered | `RuleSet2_AminoAcidCutPositionNegative_Throws` (3), `RuleSet2_AminoAcidCutPositionInRange_IsAccepted` (4) |
| M18 (Rule Set 1 contract) | ✅ Covered | `Doench2014_ReferenceCrossCheckRows_MatchTheReferenceExactly` (11), `Doench2014_GcTerm_UsesBothPublishedCoefficientsAtTheTenBoundary`, `Doench2014_NullOrEmpty_Throws`, `Doench2014_NonNggPam_Throws` |
| M19 (SaCas9 hit score) | ✅ Covered | `CalculateSaCas9HitScore_PerfectMatch_Returns100`, `..._ReferenceOracles_Match` (2), `..._CoefficientOrientation_NotReversed`, `..._MismatchAtFivePrimeMostPosition_Ignored`, `..._InvalidInputs_Throw`, `..._LowerCaseInput_TreatedAsUpperCase` |
| M20 (SaCas9 genome aggregate) | ✅ Covered | `CalculateMitSpecificityScore_SaCas9_21ntGuide_RealOffTarget_IsNot100`, `..._SaCas9_21ntGuide_TwoMismatches_MatchesReference`, `..._SaCas9_MismatchAtFivePrimeEnd_ScoredAsPerfectHit` |
| M21 (PAM-before refusal) | ✅ Covered | `CalculateMitSpecificityScore_Cas12aClass_WithRealOffTargets_Throws`, `CalculateMitSpecificityScore_PamBeforeTargetSystems_AlwaysThrow` (4), `CfdSpecificity_Genome_SaCas9_Refused` |
| M22 (CFD `maxMismatches`) | ✅ Covered | `CfdSpecificity_Genome_MaxMismatchesOutOfRange_ThrowsArgumentOutOfRange` (3), `CfdSpecificity_Genome_InvalidArguments_Throw` |
| M23 (`CfdAggregation`) | ✅ Covered | `CfdSpecificity_MitStyle_TwoHits_MatchesReference`, `CfdSpecificity_ParkinsonNormalized_TwoHits_MatchesReference`, `CfdSpecificity_SingleUnitHit_BothAggregationsMatchReference`, `CfdSpecificity_MitStyle_EmptyAndZeroSum_Returns100`, `CfdSpecificity_ParkinsonNormalized_ZeroSum_ThrowsArgumentOutOfRange` |
| M24 (alternative-PAM weight) | ✅ Covered | `CalculateMitSpecificityScore_Annotated_*` (5), `GetAlternativeOffTargetPams_MatchesReferenceTable`, `CalculateMitSpecificityScoreWithAlternativePams_FoldsInTheNagScan`, `..._NoScannableAlternative_EqualsSingleSystem`, `FindOffTargetsWithPam_MatchesFindOffTargets_AndValidatesItsArguments` |
| M25 (`NGA` / `NNGRRN`) | ✅ Covered | `CalculateMitSpecificityScoreWithAlternativePams_NgaHit_IsScannedAndDownweighted`, `..._SaCas9_NnggrrnHit_IsScannedUnweighted_AndNotDoubleCounted`, `NewAlternativePamSystems_HaveReferenceMetadata_AndExistingEnumValuesAreUnchanged` |
| M26 (`minAlternativePamHitScore`) | ✅ Covered | `CalculateMitSpecificityScoreWithAlternativePams_MinAlternativePamHitScore_FiltersAlternativeHitsOnly`, `..._SaCas9_CutoffAppliesToMainPamHitsToo` |
| M27 (CFD alt-PAM aggregate) | ✅ Covered | `CalculateCfdSpecificityScoreWithAlternativePams_NagOnlyGuide_IsNoLongerPerfect`, `..._NgaOnlyGuide_IsScanned`, `..._MixedGenome_AppliesTheMitHitScoreCutoff`, `..._KeepsTheCfdGuards`; MCP `CrisprCfdSpecificityScore_IncludeAlternativePams_DefaultsOff_AndScansNagNga_AndForwardsTheCutoff` |
| M28 (on-target interval) | ✅ Covered | `OnTargetInterval_*` (6 methods, 12 cases) in `CrisprDesigner_OffTarget_Tests.cs`; MCP `SpecificityTools_OnTargetInterval_IsForwarded_AndDefaultStaysAllExactDropped` |

### 5.3 Consolidation Plan

- **Canonical files:** keep one fixture per model (`CrisprDesigner_Doench2014_Tests`, `CrisprDesigner_RuleSet2_Tests`, `CrisprDesigner_Cfd_Tests`, and the MIT block inside `CrisprDesigner_OffTarget_Tests`). No merge — the per-model split mirrors the distinct sources and keeps fixtures focused.
- **Remove:** nothing.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `CrisprDesigner_Doench2014_Tests.cs` | Rule Set 1 | 11 methods (23 cases) |
| `CrisprDesigner_RuleSet2_Tests.cs` | Rule Set 2 / Azimuth | 22 methods (44 cases) |
| `CrisprDesigner_Cfd_Tests.cs` | CFD single-hit and aggregates | 40 methods (45 cases) |
| `CrisprDesigner_OffTarget_Tests.cs` (MIT block + F9 / F14 / F50 / round-2 E1 regions) | MIT/Hsu, SaCas9, alternative-PAM | 12 + 5 + 6 + 10 + 5 = 38 methods |

### 5.5 Phase 7 Work Queue

No ❌ or ⚠ rows in §5.2 — all MUST/SHOULD/COULD cases (M1-M28, S1-S5, C1-C2) are implemented.

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| — | (none) | — | — | — |

**Total items:** 0 · **✅ Done:** 0 · **⛔ Blocked:** 0 · **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1-M14, S1-S5, C1-C2 | ✅ | All covered; full suite **6825 passed, 0 failed**; build 0 errors, 0 warnings (figure as of the original 2026-06-18 implementation). |
| M15-M28 | ✅ | Implemented by review-2026-09 batch B08 with the fixes F9-F11, F14, F35/F36, F38, F50, F66, F67, F70, F80, F81 (per-fix results in `docs/Validation/review-2026-09/B08.md`); this spec update (2026-10-21) only records them and did not re-run the suite. |

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
