# Test Specification: RNA-ACCESS-001

**Test Unit ID:** RNA-ACCESS-001
**Area:** RnaStructure
**Algorithm:** McCaskill Unpaired (Accessibility) Probabilities
**Status:** ☑ Validated — Stage A PASS / Stage B PASS / CLEAN (2026-06-25)
**Last Updated:** 2026-06-25

> Independently re-validated under the two-stage protocol. See
> `docs/Validation/reports/RNA-ACCESS-001.md` for the full Stage A/B report, the analytic
> derivations, and the cross-verification table.

---

## 1. Evidence Summary

| # | Source | Confirms |
|---|--------|----------|
| 1 | McCaskill JS (1990) *Biopolymers* 29:1105–1119 (PMID 1695107) | Z = Σ_S exp(−E(S)/RT); P(i,j); p(S)=exp(−E/RT)/Z |
| 2 | Bernhart, Hofacker, Stadler (2006) *Bioinformatics* 22:614 + RNAplfold man page | accessibility = P(length-L window ending at i is wholly unpaired) = Z_open/Z |
| 3 | Turner-2004 NN parameters (NNDB turner04) | loop-energy model E(S); RT = 0.61626805 kcal/mol at 37 °C |

## 2. Canonical Method(s)

`CalculateUnpairedProbabilities`, `CalculateRegionUnpairedProbability`

- **Source file:** `RnaSecondaryStructure.cs` (L2556, L2672; inside DP `FillPartitionDp` L2733)
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_UnpairedProbabilities_Tests.cs`

## 3. Contract / Invariants

- R: 0 ≤ P(i,j) ≤ 1; Σ_j P(i,j) ≤ 1; p_unpaired(i) = 1 − Σ_j P(i,j) ∈ [0,1]; Z > 0.
- ΔG_ensemble = −RT·ln Z ≤ MFE (same Turner-2004 model).
- MON: region accessibility is monotone non-increasing in window length (Bernhart 2006).
- CONS: region(L=1 @ i) = per-base p_unpaired(i); both = Z_forbid(i)/Z.
- D: deterministic.

## 4. Cross-check / Differential Oracle

- **Analytic two-state pins (hand-derived from Turner-2004):**
  - GAAAC: Z = 1.0001565052764922, P(0,4) = 0.00015648078642340854,
    p_unpaired(0) = 0.9998435192135765, ΔG = −9.64416549414892e-05, region(5@4) = 1/Z.
  - CAAAAG: Z = 1.0012902114608, P(0,5) = 0.0012885489601637966,
    p_unpaired(0) = 0.9987114510398362, ΔG = −0.0007946036078507769, region(6@5) = 1/Z.
- **Independent brute-force ensemble enumeration** reproduces the GAAAC ensemble (Z, pu, P) exactly.
- **Comparison:** equality to 1e-12…1e-15.

## 5. Test cases (locked to sourced/analytic values)

| ID | Case | Asserts |
|----|------|---------|
| MCC-001/002 | GAAAC analytic | Z, P(0,4), pu, ΔG exact |
| MCC-002b | CAAAAG analytic (4-nt loop, terminal mismatch) | Z, P(0,5), pu, ΔG exact |
| MCC-003 | invariants over 5 seqs | P∈[0,1], Σ_jP≤1, pu=1−ΣP |
| MCC-004 | EFE ≤ MFE | over 4 seqs |
| MCC-007 | single base | Z=1, pu=1, no pairs, ΔG=0 |
| MCC-008 | non-ACGU (GGGNNNCCC) | finite, deterministic, N pu=1, no N pair |
| empty/null, too-short, T≤0 | edge | Z=1 / throws |
| MCC-005/006 | region bounds + GAAAC=1/Z | ∈[0,1]; =1/Z; =1 when no pair |
| MCC-009 | region monotone in length | non-increasing over L=1..14 |
| MCC-010 | region(L=1) = per-base pu | cross-method consistency |
| region OOB / len 0 | edge | throws |

Full unfiltered `dotnet test Seqeron.sln` — Failed: 0.
