# Test Specification: SPLICE-MAXENT3-001

**Test Unit ID:** SPLICE-MAXENT3-001
**Area:** Splicing
**Algorithm:** MaxEntScan score3ss (maximum-entropy 3' splice-acceptor model)
**Status:** ☑ Validated — independent Stage A/B re-validation complete (2026-06-25)
**Last Updated:** 2026-06-25

> Validated under the two-stage protocol. Stage A PASS, Stage B PASS, State CLEAN.
> See `docs/Validation/reports/SPLICE-MAXENT3-001.md`.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Yeo G, Burge CB (2004) *J Comput Biol* 11(2–3):377–394, DOI 10.1089/1066527041410418 — original MaxEntScan score3ss model |
| 2 | maxentpy (kepbod/maxentpy, MIT) `maxent.py::score3` + `score3_matrix.txt` — reference factorisation & tables (embedded as `Data/maxent_score3.txt`) |

## 2. Canonical Method(s)

`ScoreAcceptorMaxEnt(string window)`

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs` (~1003–1227)
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_AcceptorSite_Tests.cs` (region ME1–ME10)

## 3. Contract / Invariants

- Window = exactly **23 nt** (20 intron + 3 exon), AG at 0-based positions 18–19.
- Returns `log2( P_maxent / P_background )` in **bits** (finite, deterministic).
- T≡U; case-insensitive. Non-AG dinucleotide is **scored** (penalised), not rejected.
- Length ≠ 23, non-A/C/G/T(/U) → `ArgumentException`; null → `ArgumentNullException`.

## 4. Cross-check / Differential Oracle

- **Reference:** MaxEntScan score3 (maxentpy port) — documented worked examples:
  `…tAGgga` = 2.8867730651152104 (→2.89); `…cAGtgg` = 8.190965 (→8.19); `…tAGcaa` = −0.080278 (→−0.08).
- **Comparison:** exact score (C# reproduces all to <5e-7).

## 5. Validation Outcome

- [x] Stage A — sources retrieved; factorisation, constants, tables confirmed against maxentpy + Yeo & Burge.
- [x] Stage B — implementation line-for-line equivalent to maxentpy `score3`; differential oracle (8 windows) matches to <5e-7.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0 (Genomics 18762 passed).
- [x] ☑ in `ALGORITHMS_CHECKLIST_V2.md`. (docs/checklists/*.md untouched — separate campaigns.)
