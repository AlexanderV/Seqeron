# Validation Report: REP-TANDEM-001 — Tandem Repeat Detection / Tandem Repeat Summary

- **Validated:** 2026-09-30 (review campaign 2026-09, batch B04: F8–F9 + completeness audit WP3); first pass 2026-06-24 superseded
- **Area:** Repeats
- **Methods:** B04 part — `RepeatFinder.GetTandemRepeatSummary(DnaSequence, int minRepeats = 3)` and (WP3)
  `GetTandemRepeatSummary(DnaSequence, IReadOnlyDictionary<int,int> minRepeatsByUnitLength)`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs:3376/3397`), `GetCanonicalMotifClass` / `GetCanonicalMotifFrequencies`
  (:571/621), `GetStandardMotif` / `GetStandardMotifFrequencies`. The general detector `GenomicAnalyzer.FindTandemRepeats`
  is owned by batch B09 (GENOMIC-TANDEM-001) — see "Cross-batch" below.
- **Stage A verdict:** PASS-WITH-NOTES (summary fields sourced: MISA `.statistics`, Krait `statistics.py`; motif classes: MISA, Krait `motif.py`)
- **Stage B verdict:** FAIL → fixed (B04 F8 penta/hexa classes, F9 null `LongestRepeat`); MISA per-size thresholds and canonical motif tables added and reference-verified (WP3)
- **State:** FIXED

## Stage A — Description

### Sources opened
- MISA `misa.pl` v1.0 (Thiel et al. 2003; raw GitHub mirror, run with `perl`): `.statistics` — "Total number of identified SSRs",
  "Distribution to different repeat type classes" (one row per unit size), "Frequency of identified SSR motifs" (raw motifs),
  "Frequency of classified repeat types (considering sequence complementary)": class `X/Y`, X and Y = smallest rotation of the
  motif and of its reverse complement, smaller first (AC, CA, GT, TG → `AC/GT`); default `misa.ini` `1-10 2-6 3-5 4-5 5-5 6-5`.
- Krait (Du et al. 2018; lmdu/krait `src/statistics.py`, `src/motif.py`, `src/workers.py`, `src/widgets.py`): type table Mono…Hexa,
  "Length (bp)" = `SUM(length)`; `StandardMotif(level).standard()` = first of the motif's equivalence set sorted by
  `motif_to_number` (A < T < C < G); level 1 rotations, 2 + reverse complement, 3 + complement, 4 + reverse; GUI default level 3.
- Wikipedia "Tandem repeat" / Benson 1999 (definition, k ≥ 2; unchanged from the first pass).

### Definition / conventions confirmed
- Summary = aggregation of `FindMicrosatellites` (1–6 bp, maximal primitive runs): six class counts sum to `TotalRepeats`;
  `TotalRepeatBases` = Σ lengths; `PercentageOfSequence` = union coverage; `LongestRepeat` null when none; `MostFrequentUnit` = raw unit.
- Per-unit-size thresholds (MISA `definition`) are a tool parameter; canonical tables are separate results (the raw
  `MostFrequentUnit` is kept, as MISA keeps its raw motif table next to the classified one).

## Stage B — Implementation

### Defects found and fixed
1. (F8) penta/hexa class counts missing — classes did not sum to the total. (F9) `LongestRepeat` a default record instead of null.

### Additions (audit WP3)
- `GetTandemRepeatSummary(DnaSequence, IReadOnlyDictionary<int,int>)` (unit lengths 1–6, copies ≥ 2, eager validation), shared
  aggregation with the uniform overload; `MisaDefaultMinRepeats`.
- `GetCanonicalMotifClass` (MISA class name), `GetStandardMotif(motif, level 0–4)` (Krait), frequency tables over any SSR list
  (reverse complement via canonical `DnaSequence.GetReverseComplementString`, complement via `SequenceExtensions.GetComplementBase`).

### Cross-verification (0 mismatches)
| Check | Reference | Cases |
|---|---|---|
| summary fields (uniform) | Python brute-force maximal runs + aggregation | 9 000 random (F8) |
| per-class totals | real `perl misa.pl` | minRepeats 4, 5 identical (F8); per-size configurations: differences only via REP-STR-001 run conventions (default ini: misa.pl 6 639 SSRs, this library 6 641) |
| MISA class name | real `perl misa.pl` `.statistics` (one run per motif) | all 5 356 primitive motifs of 1–6 bp |
| classified table from misa.pl's own SSR list | misa.pl `.statistics` | 6 `misa.ini` configurations, 6 048 sequences |
| Krait standard motif | Krait `StandardMotif.standard()` in Python (fresh cache per level — Krait's cache is a class attribute shared by all levels) | all 5 460 motifs of 1–6 bp × levels 0–4 |

Locked values: test sequence with (AC)6 ×2, (GT)6, (TG)6, (ACAT)5, (A)13, (T)13 → misa.pl classes A/T 2, AC/GT 4, ACAT/ATGT 1;
Krait A 2, AC 4, ATAC 1; Krait levels e.g. CTG → CTG / TGC / AGC / ACG / ACG; MCP doc example 3 → A/T 4, C/G 4, AAAAG/CTTTT 1, AACCCT/AGGGTT 1 (misa.pl).

### Tests
`RepeatFinderTests` (summary), `RepeatFinderMutationTests`, `RepeatFinder_MisaCompound_Tests` (D10–D12), property tests
(`RepeatFinderProperties`, `B04ComplexityRepeatsProperties`, `RepStrMisaProperties`), MCP `TandemRepeatSummaryTests`
(`misaThresholds`, `canonicalMotifCounts`).

## Verdict & follow-ups
- Stage A: PASS-WITH-NOTES. Stage B: FAIL → fixed. **State: FIXED**; MISA thresholds and canonical motif tables implemented and reference-identical.
- MCP: `tandem_repeat_summary` delegates; optional `misaThresholds`, additive `canonicalMotifCounts` output.
- Cross-batch (B09, `GenomicAnalyzer.cs`): `FindTandemRepeats` validates its parameters but reports non-primitive units
  (`ATATATAT` → `AT×4` and `ATAT×2`) with an O(n²·m) substring scan; delegation to `RepeatFinder.FindMicrosatellites` is requested in the B04 report.
