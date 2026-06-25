# Test Specification: CHROM-HOR-001

**Test Unit ID:** CHROM-HOR-001
**Area:** Chromosome
**Algorithm:** Higher-Order Repeat (HOR) Detection
**Status:** ☑ Complete — independently validated (Stage A ✅ PASS / Stage B ✅ PASS) 2026-06-25
**Last Updated:** 2026-06-25

> **Validated 2026-06-25** under the two-stage protocol (report `docs/Validation/reports/CHROM-HOR-001.md`).
> Stage A re-grounded against McNulty & Sullivan 2018 (PMC6121732) and Rosandić 2024 (PMC11050224) retrieved this
> session; Stage B cross-checked with an independent k=4×m=7 HOR harness + the in-tree fixture (period = k×171,
> copy number = ⌊monomers/k⌋, inter-HOR ≥ intra-HOR identity, non-HOR array → no HOR). One test-only gap closed
> (non-ACGT trailing-partial monomer). State **CLEAN**; suprachromosomal-family assignment is a documented
> data-blocked boundary.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | McNulty & Sullivan (2018), Alkan et al. (2007) |

## 2. Canonical Method(s)

`DetectHigherOrderRepeat`

- **Source file:** `ChromosomeAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_HigherOrderRepeat_Tests.cs`

## 3. Contract / Invariants

R: inter-HOR identity ≥ intra-monomer identity; R: HOR period = k×monomer; D: deterministic

## 4. Cross-check / Differential Oracle

- **Reference:** known HOR arrays (D-region)
- **Comparison:** period + copy number agree

## 5. Validation Checklist (completed)

- [x] Stage A: sources retrieved this session (PMC6121732, PMC11050224) and confirmed verbatim — monomer 171 bp; intra-HOR 50–70%; inter-HOR <5% divergence (~95–100%); HOR period n ⇒ n×171 bp; inter ≥ intra ordering.
- [x] Stage B: implementation reviewed (`ChromosomeAnalyzer.cs:751`); cross-checked vs an independent k=4×m=7 HOR harness and a non-HOR control — period, copy number, unit length, monomer count, inter-identity, and the inter≥intra invariant all match the literature definition / hand-computation.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0 (Genomics 18780 passed).
- [x] Flip `☐ → ☑` in ROOT `ALGORITHMS_CHECKLIST_V2.md` (registry row + catalog header + Quick-Reference). (The 10 `docs/checklists/*.md` are intentionally NOT touched in this campaign.)
