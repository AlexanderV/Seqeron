# Test Specification: IMMUNE-NUSVR-001

**Test Unit ID:** IMMUNE-NUSVR-001
**Area:** Oncology
**Algorithm:** CIBERSORT ν-SVR Immune Deconvolution (+ bundled ABIS-Seq signature matrix)
**Status:** ☑ Validated — Stage A ✅ PASS / Stage B ✅ PASS / State ✅ CLEAN (2026-06-25)
**Last Updated:** 2026-06-25

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Newman AM et al. (2015) "Robust enumeration of cell subsets…", Nat Methods 12(5):453–457, doi:10.1038/nmeth.3337 — CIBERSORT ν-SVR. |
| 2 | CIBERSORT reference `CIBERSORT.R` (CoreAlg) — verbatim pipeline: z-score `X`/`y`; `svm(type="nu-regression",kernel="linear",nu=c(0.25,0.5,0.75))`; `weights=coefs·SV`; `weights[<0]=0; w=weights/sum`; `which.min(rmse)`. |
| 3 | Schölkopf B, Smola A, Williamson R, Bartlett P (2000) "New Support Vector Algorithms", Neural Comput 12(5):1207–1245 — ν-SVR dual; ν bounds SV fraction / margin errors (Thm 9). |
| 4 | Monaco G et al. (2019) Cell Reports 26(6):1627–1640.e7, doi:10.1016/j.celrep.2019.01.041 (CC BY 4.0) — ABIS-Seq matrix (Table S5; 1296 genes × 17 cell types). |
| 5 | Cancer Immunol Immunother review (2018) doi:10.1007/s00262-018-2150-z — confirms ν∈{0.25,0.5,0.75}, lowest-RMSE selection, non-negativity + sum-to-1. |

## 2. Canonical Method(s)

`DeconvoluteImmuneCellsNuSvr`, `LoadBundledAbisSignatureMatrix`, `LoadSignatureMatrix`

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs`
- **Bundled resource:** `Resources/ABIS_sigmatrixRNAseq.tsv` (embedded; CC BY 4.0)
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ImmuneAnalyzer_ImmuneInfiltration_Tests.cs`

## 3. Contract / Invariants

- R (range): fractions ≥ 0 and (when any signal) Σ = 1; zero-clip + renormalise.
- D (deterministic): identical inputs → identical fractions and ν.
- Planted truth `m = B·f` recovered within documented tolerances.
- ν selected from {0.25, 0.5, 0.75} by lowest standardised-scale RMSE.
- No-overlap / empty matrix → empty-or-zero fractions, sentinel `BestNu=0`.
- All-zero mixture → all-zero fractions (zero-SD standardisation branch; no NaN).

## 4. Cross-check / Differential Oracle

- **Reference:** scikit-learn 1.6.1 `NuSVR(kernel='linear', C=1)` on the same z-standardised problem; planted truth `m=B·f`; published ABIS Table S5 values.
- **Comparison numbers (this session):**
  - Disjoint 3×3 (planted 0.5/0.2/0.3): C# = sklearn = 0.508464 / 0.179557 / 0.311979 (|Δ| < 1e-6; sklearn ν=0.75).
  - Bundled ABIS (planted NK 0.60 / Mono-C 0.40): C# = sklearn = NK 0.650132 / Mono-C 0.349868 (|Δ| < 1e-6).
  - Tolerances: < 2e-3 vs sklearn (achieved < 1e-6); planted-truth recovery within 0.005–0.06.

## 5. Validation Checklist

- [x] Stage A: sources retrieved; CIBERSORT pipeline + ν-SVR dual confirmed against the reference R impl and Schölkopf (2000); ABIS provenance/licence confirmed.
- [x] Stage B: implementation matches; cross-checked vs scikit-learn `NuSVR` (< 1e-6) and planted truth; ABIS loader integrity (dimensions + exact values) confirmed.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0.
- [x] Flip `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` and the 10 `docs/checklists/*.md`.

## 6. Notes / Boundaries

- **LM22 not bundled** (Stanford no-redistribution): caller supplies it via `LoadSignatureMatrix`. **Exact-CIBERSORT/LM22 parity is NOT claimed.** The bundled matrix is **ABIS-Seq** (CC-BY). This is a documented, acceptable boundary — the ν-SVR engine itself is verified vs sklearn and planted truth.
- Added 4 edge-case tests this session (NSVR-S7…S10: all-zero mixture, empty matrix, partial overlap, ν-parameter effect).
