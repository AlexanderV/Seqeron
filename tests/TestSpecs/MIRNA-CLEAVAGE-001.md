# Test Specification: MIRNA-CLEAVAGE-001

**Test Unit ID:** MIRNA-CLEAVAGE-001
**Area:** MiRNA
**Algorithm:** Drosha/Dicer Cleavage-Site Prediction
**Status:** ☑ Validated — Stage A ✅ PASS / Stage B 🟡 PASS-WITH-NOTES → State CLEAN (2026-06-25)
**Last Updated:** 2026-06-25

> Independently re-validated under the two-stage protocol on 2026-06-25; see
> `docs/Validation/reports/MIRNA-CLEAVAGE-001.md`. The Han 2006 +11-bp Drosha ruler and Park 2011 22-nt
> 5'-counting Dicer rule are correct and the 5p mature reproduces miRBase hsa-miR-21-5p exactly. The 3p
> (miRNA*) span is a documented linear-geometry approximation (no hairpin folding); the trained miRDeep2
> classifier is a documented out-of-scope boundary.

---

## 1. Evidence Summary

| # | Source | Confirms |
|---|--------|----------|
| 1 | Han et al. (2006), *Cell* 125:887 (PubMed 16751099) | Drosha cleaves ~11 bp from the basal stem–ssRNA junction. |
| 2 | Park et al. (2011), *Nature* 475:201 (PubMed 21753850) | Dicer 5'-counting rule: cleaves ~22 nt from the 5' end (fixes mature length). |
| 3 | Auyeung et al. (2013), *Cell* 152:844 | CNNC motif ~17 nt (16–18 nt window) 3' of the Drosha basal cut (optional confidence). |
| 4 | Lee et al. (2003) / Han 2006 | RNase III leaves a 2-nt 3' overhang. |
| 5 | miRBase MI0000077 | hsa-miR-21-5p `UAGCUUAUCAGACUGAUGUUGA` (MIMAT0000076); hsa-miR-21-3p `CAACACCAGUCGAUGGGCUGU` (MIMAT0004494). |

## 2. Canonical Method(s)

`MiRnaAnalyzer.PredictDroshaDicerCleavage(string sequence, int basalJunction)` → `DroshaDicerCleavage?`

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/MiRnaAnalyzer_PreMiRna_Tests.cs` (DD1–DD14)

## 3. Contract / Invariants

- `DroshaCut5' = basalJunction + 11` (Han 2006); 0-based into the upper-cased (T→U) sequence.
- `mature = [DroshaCut5', DroshaCut5'+21]` → 22 nt (Park 2011).
- `starEnd = matureEnd + 2`, `starStart = starEnd − 21`; **linear-geometry** 2-nt 3' overhang (no folding).
- Deterministic. null/empty, junction ∉ [0,len), or any cut index ≥ len ⇒ `null`. No alphabet validation (non-ACGU passed through).

## 4. Cross-check / Differential Oracle

- **Reference:** miRBase hsa-mir-21 (MI0000077). With the 11-nt pri-miRNA lower stem reconstructed and junction = 0, the predicted 5p mature equals `UAGCUUAUCAGACUGAUGUUGA` (MIMAT0000076) exactly.
- **Limitation locked:** the linear star (`GCUUAUCAGACUGAUGUUGACU`) is NOT the real 3' arm hsa-miR-21-3p.

## 5. Validation Checklist

- [x] Stage A: Han/Park/Auyeung/Lee + miRBase retrieved this session; constants (11, 22, 2, CNNC 16–18) confirmed verbatim.
- [x] Stage B: implementation matches the rules; 5p cross-checked vs miRBase; 3p linear approximation documented (not a hidden defect).
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0 (Genomics 18775 passed).
- [x] `☐ → ☑` in root `ALGORITHMS_CHECKLIST_V2.md` (validation tracker + registry only; docs/checklists unchanged).
