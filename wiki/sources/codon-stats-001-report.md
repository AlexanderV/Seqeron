---
type: source
title: "Validation report: CODON-STATS-001 (Codon Usage Statistics — CodonUsageAnalyzer.GetStatistics + CalculateCai, positional GC1/GC2/GC3/GC3s)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/CODON-STATS-001.md
sources:
  - docs/Validation/reports/CODON-STATS-001.md
source_commit: 518339cc81914ec10d51a388064f311afb0abd4f
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CODON-STATS-001

The two-stage **validation write-up** for test unit **CODON-STATS-001** (Codon Usage Statistics —
the codon-family **aggregation / reporting** method `CodonUsageAnalyzer.GetStatistics` together with
`CalculateCai`), validated 2026-06-15. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's independent **verdict** on both the algorithm
description (Stage A) and the shipped code (Stage B). The wider campaign is
[[validation-and-testing]]; the aggregated measures are synthesized in
[[relative-synonymous-codon-usage]], [[codon-adaptation-index]] and [[effective-number-of-codons]],
and [[test-unit-registry]] defines the unit. Distinct from the pre-implementation
[[codon-stats-001-evidence]] artifact (sourced from `docs/Evidence/`), which is the fuller
description of what `GetStatistics` bundles; this page records the *re-validation verdict* and the
independent numeric cross-checks.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: ✅ CLEAN.** **No algorithm defect** — the
code realises the validated CAI/GC/RSCU formulas exactly and every externally-sourced expected value
was reproduced. Two **test-quality** defects were found and **completely fixed in-session**: a
bounds-only `S6` CAI assertion (`>=0`, `<=1`) strengthened to the exact geometric mean
`0.47706538020472955` (`Within(1e-10)`), and a documented but **untested** "non-ACGT codon skipped"
edge case closed with `GetStatistics_NonAcgtCodon_IsSkipped` (`CTGNNNGTT` → 2 codons, no `NNN` key).
Full unfiltered suite **6528 passed / 0 failed** (1 benchmark skipped by design), `dotnet build`
0 errors, changed test file warning-free.

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs`:

- `CalculateCaiCore` (`:142–184`) — builds `w = ref / family-max`, skipping stops (`'*'`) and
  single-codon families (`Count == 1`), then `exp(logSum / count)` over codons with `w > 0`.
  Realises `CAI = exp[(1/L) Σ ln w_i]` and the seqinr/CodonW exclusion rule (non-synonymous +
  termination codons excluded).
- `GetStatisticsCore` (`:389–443`) — per-position GC counts; **GC3s** via
  `IsSynonymousAtThirdPosition` (degeneracy > 1, excludes `'*'`) — matches the Peden 1999 §1.8.2.1.3
  definition ("G or C at the third position of synonymous codons, excluding Met, Trp and stop").
- `EColiOptimalCodons` (`:195–213`) — 64 entries = 61 Biopython `SharpEcoliIndex` `w` + 3 stops at
  0.0. `HumanOptimalCodons` (`:224–242`) — Kazusa-derived RSCU. Both verified by independent recompute.
- Variant contract: `string` and `DnaSequence` overloads delegate to the same `*Core`; string
  overloads short-circuit null/empty to zeroed stats / CAI 0, `DnaSequence` overloads throw
  `ArgumentNullException` on null seq or null reference.

## Stage A — description (algorithm faithfulness)

Every formula was confirmed against independently-fetched primary sources this session:

- **CAI** `w_i = f_i/max(f_j)`, `CAI = exp[(1/L) Σ ln w_i]` — matches Wikipedia *CAI* + seqinr `cai`
  exactly. Non-synonymous (single-codon Met/Trp) + termination codons excluded — confirmed
  **verbatim** by both seqinr **and** CodonW.
- **GC3s** — "frequency of G or C at the third position of synonymous codons (excluding Met, Trp and
  termination codons)" confirmed **verbatim** against the fetched **Peden 1999 thesis §1.8.2.1.3**
  (the exact quoted string was located in the PDF via `pdftotext`), corroborated by CodonW
  `Indices.html`.
- **GC1/GC2/GC3** — per-codon-position G/C content over **all** codons (EMBOSS `cusp` "1st/2nd/3rd
  letter GC").
- **RSCU** `RSCU_j = n·x_j/Σx` (Sharp et al. 1986).

**Independent hand/Python cross-checks (from the sources, not code output):** M2 `√(1×0.122) =
0.3492849839314596`; M3 `∛(0.007×0.003×0.066) = 0.011149474795453503`; every requested E. coli `w`
value (TTT 0.296, CTG 1, GCC 0.122, GCA 0.586, CGT 1, GTC 0.066, …) matched the fetched Biopython
`SharpEcoliIndex`; human RSCU **recomputed from the fetched Kazusa per-thousand frequencies** —
CTG = 6·39.6/100.2 = 2.37126, GCC = 4·27.7/69.3 = 1.59885, GTG = 4·28.1/60.7 = 1.85173 — all match
the implementation to 4 dp.

**Findings / divergences (→ PASS-WITH-NOTES), all documented, none affecting a sourced value:**
1. **GC3s reported as a percentage (×100)** vs CodonW's fraction in [0,1] — display units only; the
   synonymous subset in numerator/denominator is exactly per Peden.
2. **Zero-`w` codons skipped**, not floored to 0.01 (Bulmer 1988, used by seqinr/EMBOSS) — affects
   only a gene using a codon entirely absent from the reference; with the bundled tables no
   synonymous `w` is 0, so all worked examples are unaffected. (Same guard the standalone
   [[codon-adaptation-index|CODON-CAI-001]] records from the "codon absent but family present" angle
   as a `1e-6` clamp.)
3. **GC3s 6-fold subtlety** — CodonW treats 6-fold families (Leu/Ser/Arg) per-block; the
   implementation and Peden's own *definition* use the simpler "exclude Met/Trp/stop" set. The test
   cases use the 4-fold Ala family, so the divergence is not exercised; the validated *definition* is
   realised faithfully. Noted, not a defect.

Edge semantics all sourced: only-Met/Trp/stop → CAI 0, GC3s 0 (empty scorable/synonymous set);
empty/null string → zeroed stats / CAI 0; null `DnaSequence`/reference → `ArgumentNullException`
(input contract, explicitly documented as not literature-specified); non-ACGT codon skipped, trailing
partial codon ignored (EMBOSS `cusp` valid-codon counting).

## Stage B — implementation

The reviewed code paths compute the Stage-A formulas exactly and every cross-checked value was
reproduced by the full suite:

| Case | Expected (sourced) | Result |
|------|--------------------|--------|
| M1 all-optimal `CTG…AAA` | CAI 1.0 | pass |
| M2 `GCTGCC` | √(1×0.122)=0.3492849839314596 | pass |
| M3 `CTAATAGTC` | ∛(·)=0.011149474795453503 | pass |
| M4 Met/Trp/stop only | CAI 0.0 | pass |
| M5 `ATGGCA` | GC3s 0 / GC3 50 | pass |
| M6 `GCCGCA` | GC3s 50 | pass |
| M7 `CTGGTTAAA` | GC1/2/3 = 66.667/0/33.333 | pass |
| M8 counts/total | 4, CTG=2 | pass |
| M9 E. coli `w` | exact (Biopython) | pass |
| M10 human RSCU | exact (Kazusa) | pass |
| M11 `RSCU(CTG)` | 6.0 | pass |
| S6 arbitrary CAI | 0.47706538020472955 (recomputed) | pass (strengthened) |

**Test-quality audit (HARD gate) — PASS after fix.** M1–M11 assert exact externally-sourced or
hand-derived values within `1e-10`; the E. coli (M9) and human (M10) tables are asserted against
values recomputed from Biopython and Kazusa this session, not echoed from code. Two defects found and
fixed: (1) `S6` was bounds-only (`>=0`,`<=1`) — strengthened to the exact geometric mean while keeping
the INV-1 bounds; (2) the documented "non-ACGT codon skipped" edge (INV-5 / doc §6.1) was untested —
added `GetStatistics_NonAcgtCodon_IsSkipped`. No assertion weakened, no tolerance widened, no test
skipped, no expected value bent. Coverage: every public method/overload and every Stage-A formula path
(CAI geo-mean, exclusions, GC1/2/3, GC3s synonymous subset, RSCU) plus edge cases (empty, null seq,
null ref, partial codon, non-ACGT, only-Met/Trp/stop) exercised.

## Findings

No algorithm defect. Two **test-quality** defects (weak `S6` bounds-only assertion; missing non-ACGT
coverage) found and **completely fixed** in-session. Stage-A PASS-WITH-NOTES is entirely the three
documented unit/edge-case choices above. **Contradictions:** none — Sharp & Li 1987 (+ Biopython),
Wikipedia, seqinr, CodonW/Peden, EMBOSS `cusp` and Kazusa agree on the formulae and the
synonymous-codon exclusion set. No open follow-ups.
