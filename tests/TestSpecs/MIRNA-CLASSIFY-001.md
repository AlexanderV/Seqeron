# Test Specification: MIRNA-CLASSIFY-001

**Test Unit ID:** MIRNA-CLASSIFY-001
**Area:** MiRNA
**Algorithm:** Pre-miRNA Structure-Feature Classifier
**Status:** ☑ Complete — independently validated 2026-06-25 (Stage A ✅ / Stage B ✅ / CLEAN)
**Last Updated:** 2026-11-25 (B13 review: F23–F26)

> **Validated.** Independently re-validated under the two-stage protocol; see
> `docs/Validation/reports/MIRNA-CLASSIFY-001.md`. miRBase positives verified against the live database, the
> AMFE/MFEI/GC/%paired features and the logistic score reproduced by hand to full double precision, di-shuffle
> negatives (two independent seeds) score below threshold, and the held-out AUC = 1.0. Hardening tests CL13/CL14
> added.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Bonnet et al. (2004), miRBase (public domain), Zhang (2006) MFEI |

## 2. Canonical Method(s)

`ClassifyPreMiRna` (logistic over MFE/AMFE/MFEI/GC/%paired) — an **in-repo model, not a published classifier** (13 positives, 52 shuffles, 5 features)

`ClassifyPreMiRnaTripletSvm` (published Triplet-SVM, Xue et al. 2005, authors' libsvm model; `MiRnaAnalyzer.TripletSvm.cs`, model data `MiRnaAnalyzer.TripletSvmModel.cs`)

`ClassifyPreMiRnaMiPred` (published MiPred, Jiang et al. 2007, NAR 35:W339, authors' 1000-tree randomForest model; `MiRnaAnalyzer.MiPred.cs`, model data `MiRnaAnalyzer.MiPredModel.cs`; opt-in, additive)

- **Source file:** `MiRnaAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs`

## 3. Contract / Invariants

R: probability ∈ [0,1]; D: deterministic; threshold split positive/negative

## 4. Cross-check / Differential Oracle

- **Reference:** held-out miRBase vs shuffled (AUC)
- **Comparison (B13 2026-11-25):** held-out accuracy 0.90, AUC 0.984375 (scikit-learn on the same rows) — CL9.

## 4a. B13 additions (2026-11-25)

| ID | Test | Oracle |
|----|------|--------|
| CL1+ | `BasePairingPropensity` = 32/72 for hsa-mir-21 | dP = \|bp\|/L (Schultes 1999; Ng & Mishra 2007) |
| CL15 | `TrainPreMiRnaClassifier` reproduces the bundled constants (1e-12) and scikit-learn (1e-8) | sklearn 1.9.1 LogisticRegression(C = 1/(1e-3·45)) |
| CL16 | `FitRidgeLogisticRegression` saturated case: b0 = 0, b = ln 2; ridge score equation | hand derivation |
| CL17 | `DinucleotideShuffle` scripted draws = uShuffle 1.1.2 output (3 cases) | `ushuffle.c` + replay randfunc |
| CL18 | low-acceptance graph `A500 G A500 C` not returned unshuffled | uniform: P(input) = 1/999 |
| CL19 | uniform over the 612 valid sequences of `ACGUGUCAGGAUCC` (χ² < 724.75) | Python Eulerian enumeration + scipy χ² |
| CL20 | `CalculateRandfoldPValue` = 21/200 (20 of 199), hsa-mir-21 1/100; guards | randfold.c counting rule (Python), ViennaRNA 2.8 |

## 4b. WP-9 / A12 — published classifier (Triplet-SVM, 2026-11-25)

| ID | Test | Oracle |
|----|------|--------|
| TS1 | 40 package examples: stem-loop span, 32 features (1e-12), decision value (1e-9), label; 37 positive / 3 negative | package `examples/` (1.txt → predict_format.txt → predict_result.txt) + libsvm 3.36 `svm_predict_values` |
| TS2 | Seqeron MFE folds of 5 pri-miRNAs + 5 shuffles (Random(2024)), min arm 22 and 12 (multi-candidate case) | package Perl scripts 2/3/4 + libsvm on the same structures |
| TS3 | default folder = `CalculateMfeStructure`; all 5 miRBase precursors called real | internal consistency |
| TS4 | `TripletSvmDecisionValue` at 3 vectors outside the training data (1e-12) | libsvm `svm_predict_values` |
| TS5 | DNA/lower case normalised; non-ACGU dropped; < min arm pairs empty; guards | package step 1/2 semantics |
| TS6b | `ComputeTripletFeatures` (named row; same test as TS6 plus the two-arm-pair `null` case) | definition in `3_step_triplet_coding_for_queries.pl`; a structure with two arm pairs returns null (package skips) |
| TS6 | hand-computed triplet coding of `GGGGAAAACCCC`/`((((....))))` | definition in `3_step_triplet_coding_for_queries.pl` |

Not committed: 300 random ViennaRNA 2.8 folds × min arm 22/12/6 (981 coded candidates) compared with the Perl scripts + libsvm: 0 differences (decision values < 1e-6, spans, sequences, structures, labels).

## 4c. WP-R2-1 — published classifier (MiPred, 2026-11-25)

Tests: `MiRnaAnalyzer_MiPred_Tests.cs` (data: `MiPredReferenceData.cs`, `MiPredOracleData.cs`).

| ID | Test | Oracle |
|----|------|--------|
| MP1/MP1b | hand-computed 32 triplet frequencies (`GAUCC`/`((.((`: V2, V11, V29 = 0.2); closing-bracket triplets not counted; DNA = RNA | `microRNAcheck_parallel.pl` counting loop (k = 1..n-2, divisor = LENGTH) |
| MP2 | 32 frequencies of 12 sequences (5 precursors, 4 shuffles, 3 pri-miRNAs) = integer counts / length | independent Python port of the Perl, ViennaRNA 2.8 structures |
| MP3/MP3b | pre-judgement (MFE <= -15, >= 18 pairs, no multiloop, 51..137 nt, second-hairpin regex) on 12 structures + every boundary | Perl port; package regex |
| MP4 | forest votes on the 10 Python-derived 34-vectors (exact) | independent numpy traversal of the forest arrays parsed from `model.RData` |
| MP5 | 60 stored training rows of `model.RData`: votes (exact), in-sample class | numpy traversal; stored labels |
| MP6/MP6b | bundled forest: ntree 1000, 56,324 nodes (28,662 terminal), every tree a consistent full binary tree (9..103 nodes); corrupt tables rejected | parse of `model.RData` (ndbigtree, nrnodes 103) |
| MP7 | 60-tree forest trained locally with R `randomForest`, exported with `getTree`, 560 rows (300 fresh, 100 training, 160 built from split points of the forest, 131 of which land exactly on a split point of a tree path; 4 exact 30:30 ties): votes = `predict(type = "vote", norm.votes = FALSE)`; `<=` rule exercised | R 4.x randomForest 4.7-1.1 |
| MP8 | end-to-end (canonical fold, 999 seeded shuffles): 5 miRBase precursors hairpin-like, p <= 0.01, voted real (781/823/836/802/784 of 1000; Python 780/823/836/801/784) | Python port + numpy shuffles + ViennaRNA |
| MP9 | 2 shuffled let-7a-1 hairpins (s7, s9; s3 excluded: canonical engine scores a 1-nt bulge 1.26 kcal/mol below ViennaRNA) passing the pre-judgement: voted pseudo; p-value agrees with the numpy randfold-shuffle p within 4 sigma | Python port |
| MP10/MP11 | pri-miRNAs: hsa-mir-21 (multiloop), hsa-let-7a-1 (140 nt) not hairpin-like -> not classified; hsa-mir-155 (125 nt) passes and votes real | package pre-judgement |
| MP12-MP14 | determinism / default seed / clean-up (`tr`) / argument contract | package semantics |
| MP15 | other classifiers unchanged | existing values |

Differences from the package (documented, not simplifications): MFE/structure from the canonical Turner-2004 engine (package: RNAfold 1.x); shuffles with `DinucleotideShuffle` (package: randfold's triplet-swap shuffle with Perl `rand`), so the p-value feature is statistical (fixed default seed 1); model licence not stated (numeric parameters only) — owner decision.

## 5. Validation Checklist (to restore ☑)

- [x] Stage A: retrieved Bonnet 2004 (di-shuffle null model), Zhang 2006 (AMFE/MFEI), miRBase (sequences verified live); formulas confirmed.
- [x] Stage B: implementation matches the source; features and logistic score reproduced by hand to full precision; held-out AUC = 1.0.
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0 (Genomics 18771 passed, 0 warnings).
- [x] Flipped `☐ → ☑` in ROOT `ALGORITHMS_CHECKLIST_V2.md` (registry + catalog + Quick-Reference). (docs/checklists/* intentionally untouched.)
