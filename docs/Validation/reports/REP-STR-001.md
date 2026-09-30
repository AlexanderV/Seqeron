# Validation Report: REP-STR-001 — Microsatellite / Short Tandem Repeat (STR) Detection

- **Validated:** 2026-09-30 (review campaign 2026-09, batch B04: F5–F7 + completeness audit WP3 F28–F33); first passes 2026-06-24/25 superseded
- **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindMicrosatellites(DnaSequence|string, minUnitLength=1, maxUnitLength=6, minRepeats=3)`
  (+ `CancellationToken` / `IProgress<double>` overloads) — `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs:40/67/86/104`,
  single core `FindMicrosatellitesCore` :240
- **MISA additions (audit WP3):** `FindMicrosatellites(DnaSequence|string, IReadOnlyDictionary<int,int> minRepeatsByUnitLength, …)` (:161/179),
  `MisaDefaultMinRepeats`, `MisaDefaultMaxInterruption`; `FindCompoundMicrosatellites` (:389–431), `AssembleCompoundMicrosatellites` (:458)
- **Related (own reports / specs):** approximate detector `FindApproximateTandemRepeats` / `ComputeBernoulliStatistics` → REP-APPROX-001;
  summary / motif classes → REP-TANDEM-001
- **Stage A verdict:** PASS-WITH-NOTES (reporting convention sourced: MISA, pytrf/Krait, maximal repetitions)
- **Stage B verdict:** FAIL → fixed (B04 F5–F7); MISA thresholds / compound SSRs / progress contract added and reference-verified (WP3, F28–F31)
- **State:** FIXED

## Stage A — Description

### Sources opened
- MISA `misa.pl` v1.0 (Thiel et al. 2003, Theor Appl Genet 106:411; raw GitHub mirror `cfljam/SSR_marker_design`; run with
  `perl`): per motif size the regex `(([acgt]{p})\2{k-1,})` scanned left to right (`/ig`), non-primitive ("false type")
  motifs rejected after matching; `misa.ini` `definition(unit_size,min_repeats): 1-10 2-6 3-5 4-5 5-5 6-5`,
  `interruptions(max_difference_for_2_SSRs): 100`; compound assembly loop (types `c`, `c*`, notation
  `(M)r<interruption>(M)r`, `*` after an overlapping component); `.misa` columns ID, SSR nr., type, SSR, size, start, end.
- pytrf 1.5.0 `src/str.c` (Krait engine; Du et al. 2018): run seeded at its start, `repeat = length / p`, `N` skipped.
- Kolpakov & Kucherov 1999 (maximal repetitions; cited). Wikipedia "Microsatellite" examples (TATATATATA, GTCGTCGTCGTCGTC, GGAT).

### Definition / conventions confirmed
- Per unit length p, each **maximal primitive run** of A/C/G/T is reported once at its left end with ⌊run/p⌋ complete copies;
  runs of different unit lengths may overlap. Case-insensitive; 0-based positions.
- MISA per-unit-size thresholds are a tool parameter (`minRepeatsByUnitLength`); the single `minRepeats` overload is the uniform case.
- Compound SSR (MISA): consecutive SSRs in start order chained when `start(i+1) − end(i) ≤ maxInterruption`
  (comparison with the previous SSR's end; compound end = last component's end); `c*` when any consecutive pair overlaps.

### Findings
- The 2026-06 description ("contained intervals suppressed, rotations AC×5 / CA×5 both reported") contradicted MISA / pytrf —
  corrected (F5). The 2026-06 TestSpec claimed "progress reporting not implemented"; the cancellable overloads do take
  `IProgress<double>` — contract now specified and tested (F31).

## Stage B — Implementation

### Defects found and fixed
1. (F5) the same locus reported once per rotation (phase); (F6) runs of `N` reported (string overload);
   (F7) `minUnitLength = 0` hang, cancellable overload unvalidated. Details: `docs/Validation/review-2026-09/B04.md`.

### Additions (audit WP3)
- Per-unit-size thresholds: one core over (unit length, minRepeats) pairs; eager validation (null / empty map, unit < 1, copies < 2).
- Compound SSRs: `AssembleCompoundsCore` reproduces misa.pl's loop (stable start order; interruption lower-cased).
- Progress: non-decreasing values in [0, 1) every 1000 visited run starts, final 1.0; token checked at the same points.

### Cross-verification (0 unexplained mismatches)
| Check | Reference | Cases |
|---|---|---|
| `FindMicrosatellites` (uniform) | brute-force maximal primitive runs; MISA regex; pytrf | 3 132 (F5): 3 132/3 132 brute force; MISA/pytrf differences all documented conventions |
| `FindMicrosatellites` (per-size map) | brute-force maximal primitive runs with per-size thresholds | 6 048 sequences, 6 `misa.ini` configurations, 72 974 SSRs: 0 |
| `AssembleCompoundMicrosatellites` fed misa.pl's SSR list (misa.pl order) | real `perl misa.pl` `.misa` rows | 6 048 sequences, 12 330 compounds: 0 |
| `FindCompoundMicrosatellites` end to end | real `perl misa.pl` | rows identical wherever the SSR lists agree; SSR-list differences: 1 220 SSRs same-size overlap < p (MISA truncates), 1 365 SSRs consumed by a rejected non-primitive match, 53 equal-start tie orders (Perl hash order) — 0 unexplained |

Locked values (misa.pl output): `ACGTATATATATATATccgtGAGAGAGAGAGAGAtttttAAAAAAAAAAAAT` → `c (TA)6tccgt(GA)7ttttt(A)12 48 4 51`;
`nnnACACACACACACACAGCAGCAGCAGCAGCAG` → `c* (AC)7(CAG)6* 31 4 34`; 100 interrupting bases joined, 101 not; interruptions 0/1;
`CCCCCCCCCGCCCCCCCCCC` → only `(C)10 11–20`; documented differences `(AAGAA)8@45` vs misa.pl `(AGAAA)8@46`, `(CTTAAA)7@129` vs `(TAAACT)6@131`.

### Tests
`RepeatFinder_Microsatellite_Tests` (F5–F7), `RepeatFinder_MisaCompound_Tests` (thresholds, compounds, progress / cancellation),
`RepeatFinderMutationTests`, `Properties/B04ComplexityRepeatsProperties` + `RepeatFinderProperties` (maximal runs),
heavy tier `Properties/RepStrMisaProperties` + `Metamorphic/RepStrMisaMetamorphicTests` (WP3); MCP `FindMicrosatellitesTests`
(`misaThresholds`, `maxCompoundInterruption`).

## Verdict & follow-ups
- Stage A: PASS-WITH-NOTES. Stage B: FAIL → fixed. **State: FIXED**; MISA thresholds, compound SSRs and the progress contract implemented and misa.pl-verified.
- MCP: `find_microsatellites` delegates (`misaThresholds`, `maxCompoundInterruption` optional parameters, additive).
- Cross-batch: `GenomicAnalyzer.FindTandemRepeats` (B09) and `GenomeAnnotator` (B11) private tandem engines should delegate to `FindMicrosatellites` (B04 report).
