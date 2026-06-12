# Validation Report: REP-TANDEM-001 — Tandem Repeat Detection

- **Validated:** 2026-06-12   **Area:** Repeats
- **Canonical method(s):** `GenomicAnalyzer.FindTandemRepeats(DnaSequence, int minUnitLength=2, int minRepetitions=2)`; summary delegate `RepeatFinder.GetTandemRepeatSummary(DnaSequence, int minRepeats=3)` (wraps `FindMicrosatellites`).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Tandem repeat** (fetched 2026-06-12): definition is exactly *"A pattern of one or more nucleotides is repeated and the repetitions are directly adjacent to each other"* (example ATTCG ATTCG ATTCG). Confirms: a tandem repeat is **adjacent (consecutive)** copies of a unit; constitutes ~8% of the human genome; >50 diseases (Huntington's etc.). No formal minimum copy count is mandated by the article, but "repeated… adjacent to each other" implies a minimum of **two** adjacent copies (one copy is not a repeat).
- **Wikipedia — Microsatellite** (fetched 2026-06-12): microsatellite / STR repeating element is "generally ten nucleotides or less" (commonly stated as 1–6 bp); minisatellites are the "longer cousins"; telomere motif TTAGGG (vertebrate, hexanucleotide); "typically repeated 5–50 times"; **forensic STRs are all tetra- or penta-nucleotide repeats**.
- **Richard et al. (2008), MMBR 72(4):686–727** — cited as the comparative review of eukaryotic DNA repeats; consistent with the micro/minisatellite unit-length classification above. (Cited in spec; class boundaries used here come from the two Wikipedia pages with their primary refs.)

### Definition / convention check
- **Tandem repeat = ≥2 directly adjacent copies of a unit.** The implementation's `minRepetitions` (default 2) realises "two or more adjacent copies"; the spec worked example uses ×3.
- **Period/unit & copy number reporting:** result reports `Unit` (the period string), `Position` (start), `Repetitions` (copy count), with `TotalLength = Unit.Length × Repetitions`. Matches the standard period+count reporting.
- **Coordinate base:** 0-based start position (test M9 confirms position 2 for `CCATGATGATGCC`).
- **Microsatellite vs minisatellite:** unit 1–6 bp = microsatellite; longer = minisatellite. The general `FindTandemRepeats` is unit-length agnostic (any `minUnitLength`); REP-STR (`FindMicrosatellites`) is the microsatellite-specific path (1–6 bp). The summary delegate restricts to 1–6 bp (microsatellite classes), which is correct for a "tandem repeat summary" focused on STRs.
- **Perfect vs approximate:** this algorithm detects **perfect** tandem repeats only (exact unit equality, `seq.Substring(pos, unitLen) == unit`). This matches the spec, which makes no approximate-repeat claim.

### Hand-computed worked example
`ATGATGATG`, minUnitLength=3, minRepetitions=3:
- unitLen=3, start=0: unit="ATG"; pos=3 "ATG"=match (reps 2), pos=6 "ATG"=match (reps 3), pos=9 out of bounds → reps=3 ≥ 3 → yield **(Unit="ATG", Position=0, Repetitions=3, TotalLength=9, FullSequence="ATGATGATG")**.
- No other unitLen yields a result for this input/params.
Matches spec expected values (unit "ATG" ×3) exactly.

### Edge-case semantics
- No repeat → empty (M5 `ACGT`); empty input → empty (M6).
- Two adjacent copies = minimal tandem (default minRepetitions=2). Confirmed by trace and M7.
- Longer unit lengths (tetra M4, penta M13, hexa S4/telomere C2) all sourced and supported.

**Stage A verdict: PASS** — definition, classification, coordinate base, and worked example match authoritative sources.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs:77-103` (`FindTandemRepeats`), struct `TandemRepeat` at `:354-368`.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs:439-470` (`GetTandemRepeatSummary` → `FindMicrosatellites`).

### Formula realised correctly?
The algorithm iterates `unitLen` from `minUnitLength` up to `seq.Length / minRepetitions`, and for each start greedily extends a maximal run of exact-equal units, emitting `(unit, start, repetitions)` when `repetitions >= minRepetitions`, then skips past the consumed run (`start = pos - unitLen`, then loop `start++`). This is the correct period-by-period perfect-tandem detector.

### Worked-example recompute vs code
| Input | params | Code result (traced) | Expected |
|-------|--------|----------------------|----------|
| ATGATGATG | (3,3) | (ATG,0,3), TotalLength 9 | ✔ unit ATG ×3 |
| CACACACA | (2,4) | (CA,0,4) | ✔ |
| AAAAA | (1,5) | (A,0,5) | ✔ |
| GATAGATAGATA | (4,3) | (GATA,0,3), TL 12 | ✔ |
| AAAGAAAAGAAAAGA | (5,3) | (AAAGA,0,3), TL 15 | ✔ |
| CCATGATGATGCC | (3,3) | (ATG,2,3) — 0-based | ✔ |
| ATATAT | (2,2) | (AT,0,3) only | ✔ |
| ACGT | (2,2) | empty | ✔ |

### Invariants verified
- `TotalLength = Unit.Length × Repetitions` (property `TotalLength`, M11). ✔
- `FullSequence` reconstructs `Unit × Repetitions` (M12). ✔
- `Position + TotalLength ≤ seqLength` (P3): proven — `while` only extends while a full unit fits in bounds, so `pos ≤ seq.Length`, and `start + TotalLength = pos`. ✔
- All results satisfy `Repetitions ≥ minRepetitions` (P1) and `Unit.Length ≥ minUnitLength` (P2). ✔

### Variant / delegate consistency
`GetTandemRepeatSummary` aggregates `FindMicrosatellites(seq, 1, 6, minRepeats)` by `RepeatType`, summing `TotalLength`, computing `% of sequence`, and identifying longest repeat / most-frequent unit. Test D3 (`AAACAGCAGCAGCAGCAGCAGAAA`) → longest CAG ×6 = 18 bp. ✔ D4 (`AAAAAATTTTTGGGGGCCCCC`) → 4 mononucleotide runs. ✔ Note: the summary uses the microsatellite (REP-STR) path, not `FindTandemRepeats`; this is the documented contract and is appropriate since a "tandem repeat summary" is STR-focused (1–6 bp classes).

### Notes (design choices, not defects)
- **Per-period enumeration:** a perfectly periodic run is reported once per applicable `unitLen` (e.g. `AAAA` with minUnitLength=1 yields both `A×4` and `AA×2`). This is mathematically correct — both are genuine tandem repeats — and tests isolate periods via `minUnitLength`. No deduplication is claimed or required; not a defect.

### Test quality audit
27-item spec realised as 42 executable tests (the fixture splits some MUST items and includes property tests). Assertions check exact sourced values (units, positions, counts, total lengths) — not "no-throw" tautologies. Edge cases (empty, no-repeat, thresholds, boundaries) covered. Tests are deterministic.

**Stage B verdict: PASS** — code faithfully realises the validated definition; all traced values match.

## Verdict & follow-ups
- **STATE: CLEAN.** No defect found. Build succeeds; Tandem filter = 42 passed; full `Seqeron.Genomics.Tests` suite = 4461 passed, 0 failed (baseline). No code changes required.
