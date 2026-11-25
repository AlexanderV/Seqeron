# Test Specification: SPLICE-MAXENT3-001

**Test Unit ID:** SPLICE-MAXENT3-001
**Area:** Splicing
**Algorithm:** MaxEntScan score3ss (maximum-entropy 3' splice-acceptor model)
**Status:** ☑ Validated — independent Stage A/B re-validation complete (2026-06-25)
**Last Updated:** 2026-11-25

> Validated under the two-stage protocol. Stage A PASS, Stage B PASS, State CLEAN.
> See `docs/Validation/reports/SPLICE-MAXENT3-001.md`.

---

## 1. Evidence Summary

| # | Source |
|---|--------|
| 1 | Yeo G, Burge CB (2004) *J Comput Biol* 11(2–3):377–394, DOI 10.1089/1066527041410418 — original MaxEntScan score3ss model |
| 2 | maxentpy (kepbod/maxentpy, MIT) `maxent.py::score3` + `score3_matrix.txt` — reference factorisation & tables (embedded as `Data/maxent_score3.txt`) |

## 2. Canonical Method(s)

`ScoreAcceptorMaxEnt(string window)`; `CalculateMaxEntScore(motif, SpliceSiteType.Acceptor)` (MCP `maxent_score`) delegates to it.

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs` (~1003–1227)
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/SpliceSitePredictor_AcceptorSite_Tests.cs` (region ME1–ME10)

## 3. Contract / Invariants

- Window = exactly **23 nt** (20 intron + 3 exon), AG at 0-based positions 18–19.
- Returns `log2( P_maxent / P_background )` in **bits** (finite, deterministic).
- T≡U; case-insensitive. Non-AG dinucleotide is **scored** (penalised), not rejected.
- Length ≠ 23, non-A/C/G/T(/U) → `ArgumentException`; null → `ArgumentNullException`.

## 4. Cross-check / Differential Oracle

- **Reference:** MaxEntScan score3 (maxentpy port) — documented worked examples:
  `…tAGgga` = 2.8867730651152104 (→2.89); `…cAGtgg` = 8.190965 (→8.19); `…tAGcaa` = −0.080278 (→−0.08).
- **Comparison:** exact score (C# reproduces all bit-for-bit; see B11 review below).

### B11 review (2026-11-25)

- **Oracle:** 200 025 windows (100 000 random + 100 000 AG 23-mers, 16 dinucleotides, examples, homopolymers):
  C# bit-identical to maxentpy `score3` (after F22: log2 = ln/ln2; before, 62 349 differed by 1 ulp) and
  identical at 2 dp to the score3.pl CPAN port (original me2x3acc1-9 tables = maxentpy table values).
- `CalculateMaxEntScore(motif, Acceptor)` and MCP `maxent_score` (Acceptor) delegate to `ScoreAcceptorMaxEnt`
  (F21; formerly a PWM log-weight sum, −26.81 for the 2.89 window): 23-nt A/C/G/T(/U) required, else
  `ArgumentException`; empty → 0.
- Tests: `ScoreAcceptorMaxEnt_ReferenceSweep_MatchesMaxEntScan` (ME10, 16 cases, exact double + 2 dp),
  `CalculateMaxEntScore_AcceptorSite_ReturnsMaxEntScanScore3`, `_AcceptorWrongLengthOrAlphabet_Throws`,
  mutation `_AcceptorBranch_IsMaxEntScanScore3`, MCP `MaxentScore_Acceptor_IsMaxEntScanScore3`,
  `_AcceptorNot23NtOrInvalid_Throws`.

## 5. Validation Outcome

- [x] Stage A — sources retrieved; factorisation, constants, tables confirmed against maxentpy + Yeo & Burge.
- [x] Stage B — implementation line-for-line equivalent to maxentpy `score3`; differential oracle (8 windows) matches to <5e-7.
- [x] Full unfiltered `dotnet test Seqeron.sln -c Debug` — Failed: 0 (Genomics 18762 passed).
- [x] ☑ in `ALGORITHMS_CHECKLIST_V2.md`. (docs/checklists/*.md untouched — separate campaigns.)
