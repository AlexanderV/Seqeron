# Test Specification: CHROM-ALPHASAT-001

**Test Unit ID:** CHROM-ALPHASAT-001
**Area:** Chromosome
**Algorithm:** Alpha-Satellite Monomer Detection
**Status:** ☑ Complete — independently validated 2026-06-25 (Stage A ✅ / Stage B ✅ / CLEAN)
**Last Updated:** 2026-06-25

> Independently re-validated under the two-stage protocol against external first sources retrieved this
> session. See `docs/Validation/reports/CHROM-ALPHASAT-001.md`.

---

## 1. Evidence Summary

| # | Source | Confirms |
|---|--------|----------|
| 1 | Willard HF (1985); Waye JS & Willard HF (1987); review Hartley & O'Neill (2019), PMC6121732 | 171-bp AT-rich alphoid monomer; 50–70% intra-array monomer identity |
| 2 | Masumoto H et al. (1989), J Cell Biol 109(4):1963–1973; PMC4843215 | 17-bp CENP-B box consensus `YTTCGTTGGAARCGGGA` (Y=C/T, R=A/G), core `TTCG…CGGG` |

Note: PMC6121732 renders the box as `5'-T/CTCGTTGGAAA/GCGGGA-3'` (16 bp, a typo dropping one T);
the canonical Masumoto/PMC4843215 17-bp form is implemented.

## 2. Canonical Method(s)

`DetectAlphaSatellite(string)`, `FindCenpBBoxes(string)`; constants `AlphaSatelliteMonomerLength` (171),
`CenpBBoxConsensus` (`YTTCGTTGGAARCGGGA`).

- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs`
- **Test fixture:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_AlphaSatellite_Tests.cs`

## 3. Contract / Invariants

- R: monomer period = 171 bp (scanned window 166–176); perfect tandem → periodicity 1.0.
- R: `IsAlphaSatellite` ⇔ periodicity ≥ 0.50 AND AT content > 0.50 (AT over ACGT bases only).
- R: CENP-B box = 17-bp IUPAC match (Y=C/T, R=A/G, fixed positions exact); 0-based ascending positions.
- D: deterministic; case-insensitive; non-ACGT excluded from AT denominator.

## 4. Cross-check / Differential Oracle

Independent Python re-implementation (consensus + periodicity) reproduced every expected value:
period 171 / periodicity 1.0 / AT 100/171 / 10 boxes in a 10-copy array / four IUPAC corners at index 0 /
non-ACGT AT = 100/166. See report for the full table.

## 5. Validation Checklist

- [x] Stage A: external first sources retrieved; 171-bp monomer + 17-bp consensus + AT-richness confirmed.
- [x] Stage B: implementation realises the description; cross-checked vs independent oracle.
- [x] Full unfiltered `dotnet test Seqeron.sln` — Failed: 0 (Genomics 18779 passed).
- [x] Flipped `☐ → ☑` in `ALGORITHMS_CHECKLIST_V2.md` (registry + catalog) and the validation tracker.

## 6. Documented boundary

Suprachromosomal-family / chromosome-specific HOR family assignment is **data-blocked** — it requires
curated reference HOR libraries (T2T/CHM13) not in the repo. Out of scope for this unit (monomer
periodicity + AT-richness + CENP-B box detection); acceptable per the validation criteria.
