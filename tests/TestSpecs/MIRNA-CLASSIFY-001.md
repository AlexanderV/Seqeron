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

## 5. Validation Checklist (to restore ☑)

- [x] Stage A: retrieved Bonnet 2004 (di-shuffle null model), Zhang 2006 (AMFE/MFEI), miRBase (sequences verified live); formulas confirmed.
- [x] Stage B: implementation matches the source; features and logistic score reproduced by hand to full precision; held-out AUC = 1.0.
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0 (Genomics 18771 passed, 0 warnings).
- [x] Flipped `☐ → ☑` in ROOT `ALGORITHMS_CHECKLIST_V2.md` (registry + catalog + Quick-Reference). (docs/checklists/* intentionally untouched.)
