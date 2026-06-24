# Validation Report: REP-TANDEM-001 — General Tandem Repeat Detection

- **Validated:** 2026-06-24   **Area:** Repeats
- **Canonical method(s):** `GenomicAnalyzer.FindTandemRepeats(DnaSequence, int minUnitLength = 2, int minRepetitions = 2)` → `IEnumerable<TandemRepeat>` (`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs:115`); summary delegate `RepeatFinder.GetTandemRepeatSummary(DnaSequence, int minRepeats = 3)` (`RepeatFinder.cs:460`).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** CLEAN
- **Relationship:** Same implementation previously validated under GENOMIC-TANDEM-001 (documented duplicate registry entry, 2026-06-16). This is an independent re-validation in a fresh context with sources re-retrieved this session.

## Distinction from REP-STR-001

- **REP-STR-001** = `RepeatFinder.FindMicrosatellites` — short-tandem-repeat (STR/microsatellite) detector. It restricts to short units, classifies by `RepeatType` (mono/di/tri/tetra…), and applies **redundant-unit filtering** (e.g. `ATAT` is reported as `AT`×2, not `ATAT`×1).
- **REP-TANDEM-001** = `GenomicAnalyzer.FindTandemRepeats` — a **general** brute-force tandem detector across an arbitrary period range (`minUnitLength .. seq.Length/minRepetitions`). It does **not** classify and does **not** deduplicate periods: every period interpretation that satisfies the definition is reported (e.g. `AAAA` → both `A`×4 and `AA`×2). The two methods are genuinely distinct.
- The `GetTandemRepeatSummary` "delegate" listed in the TestSpec actually wraps `FindMicrosatellites` (confirmed `RepeatFinder.cs:466`), exactly as the TestSpec §"Summary Tests" states ("wraps FindMicrosatellites"). It is a summary over the STR detector, not over `FindTandemRepeats`; not a defect, but noted.

## Stage A — Description

### Sources opened this session (URLs + extracted numbers)

1. **Wikipedia, "Tandem repeat"** — https://en.wikipedia.org/wiki/Tandem_repeat (accessed 2026-06-24).
   - Definition (verbatim): "a pattern of one or more nucleotides is repeated and the repetitions are directly adjacent to each other."
   - Worked example (verbatim): "ATTCG ATTCG ATTCG, in which the sequence ATTCG is repeated three times." → unit `ATTCG`, period 5, copy number 3, total length 15, 0-based start 0.
   - Classification by unit length: microsatellite/STR (< ~10 nt), minisatellite (10–60 nt), macrosatellite (~1000 nt).
   - Detection: "can be efficiently detected using suffix trees or suffix arrays."
2. **Benson, G. (1999), Tandem Repeats Finder, NAR 27(2):573–580** — https://academic.oup.com/nar/article/27/2/573/1061099 (accessed 2026-06-24). DOI 10.1093/nar/27.2.573.
   - Formal definition (verbatim): "two or more contiguous, **approximate** copies of a pattern of nucleotides." → minimum copy number **k ≥ 2**.
   - Period size = "the most common matching distance between corresponding characters in the alignment" (= unit length for an exact repeat); copy number = number of aligned pattern copies; consensus size may differ from period size.
   - Allows approximate (non-identical) copies: "typically, only approximate tandem copies are present"; uses percent identity / indel frequency, not exact equality.

### Formula / model check

The algorithm detects an **exact** tandem repeat: `S[p .. p + k·|U|) = U^k` with `k ≥ 2`, reporting period = `|U|` and copy number = `k`. This is the exact-only specialization of Benson's definition; over exact tandems Benson's "approximate copies" model and this model coincide, so the simplification is a sound subset (not a wrong rule). The minimum copy threshold `k ≥ 2` matches Benson exactly. `TotalLength = period × copies` and "block within bounds" are genuine mathematical invariants.

### Edge-case semantics

- Empty → no hits; no-tandem → empty; whole-sequence tandem → one spanning result.
- A region satisfying multiple periods is reported once per qualifying period (e.g. `AAAA` → `A`×4 and `AA`×2). Each is itself a valid tandem repeat under the formal definition, so reporting all is sound. (Contrast REP-STR-001, which deduplicates by reducing to the primitive period.)
- Guards: `minUnitLength ≥ 1` (a 0-length unit never terminates the scan) and `minRepetitions ≥ 2` (Benson's k ≥ 2; also avoids div-by-zero in the unit-length bound). Both throw `ArgumentOutOfRangeException`; `null` sequence throws `ArgumentNullException`. These are sourced/justified, not implementation-defined.

### Independent cross-check (hand computation, traced to sources)

`ATTCGATTCGATTCG` (Wikipedia worked example), minUnitLength 5, minRepetitions 2 — hand trace:
- unitLen 5, bound `5 ≤ 15/2 = 7` ✓; start 0 unit `ATTCG`; pos 5,10 match, pos 15 stops → 3 copies ≥ 2 → yield `(ATTCG, 0, 3)`; restart `start = 15 − 5 = 10`, loop end exceeds bound.
- unitLen 6/7 yield no period match.
- Result: exactly one `(ATTCG, pos 0, 3 copies, total 15)` — **matches Wikipedia exactly**.

`AAAA`, minUnit 1: periods 1 and 2 both satisfy the definition → `A`×4 and `AA`×2 both genuine.

**Stage A verdict: PASS.** Definition, period/copy-number conventions, 0-based coordinate, k ≥ 2 threshold, and edge cases are all correct and source-backed.

## Stage B — Implementation

### Code path reviewed

`FindTandemRepeats` (`GenomicAnalyzer.cs:115–126`) → `FindTandemRepeatsCore` (`:128–154`); `struct TandemRepeat` (`:568–582`). Brute force: for `unitLen` in `minUnitLength .. seq.Length/minRepetitions`, for each `start ≤ seq.Length − unitLen·minRepetitions`, take candidate unit, count consecutive exact copies, yield if `≥ minRepetitions`, then `start = pos − unitLen` to skip to the end of the matched block within that unit-length pass. `TotalLength = Unit.Length × Repetitions`; `FullSequence = Unit repeated Repetitions times`.

### Formula realised correctly?

Yes. Direct execution this session (throwaway harness referencing the real Analysis assembly, since removed) against externally-derived expectations:

| Input | minUnit/minReps | Source/derived expectation | Code output | Match |
|-------|-----------------|----------------------------|-------------|-------|
| `ATTCGATTCGATTCG` | 5 / 2 | ATTCG, pos 0, 3 copies, total 15 (Wikipedia) | ATTCG/0/3/15 | ✅ |
| `AAAA` | 1 / 2 | A×4 and AA×2 (both valid periods) | A/0/4/4 + AA/0/2/4 | ✅ |
| `ATGATGCCATGATGATG` | 3 / 2 | ATG×2 @0 and ATG×3 @8 (skip+restart) | ATG/0/2/6 + ATG/8/3/9 | ✅ |
| `AAATTT` | 1 / 3 | A×3 @0 and T×3 @3 (adjacent distinct runs) | A/0/3/3 + T/3/3/3 | ✅ |
| `` (empty) | 2 / 2 | empty | empty | ✅ |

### Variant / delegate consistency

No `*Fast` variant. `GetTandemRepeatSummary` is a summary over `FindMicrosatellites` (REP-STR-001 surface), as the TestSpec states; it is internally consistent (`TotalRepeatBases` = sum of repeat lengths; `PercentageOfSequence` uses the union of covered spans to stay ≤ 100, `CountCoveredBases` `RepeatFinder.cs:504`). `TotalLength`/`FullSequence` derived properties verified consistent with `Unit`/`Repetitions`.

### Numerical robustness

Integer indexing only; outer/inner loops bounded by `seq.Length`; restart index can only decrease the work; no overflow or precision concerns on stated ranges. Thresholds are guarded against the non-terminating/div-by-zero cases.

### Test quality audit (against sources, not code)

Canonical fixture `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_TandemRepeat_Tests.cs` — 24 tests, all green. Correctness assertions use **exact** sourced values (unit, 0-based position, copy count, TotalLength, FullSequence); both `minUnitLength` and `minRepetitions` filters exercised; empty/no-repeat boundaries; long repeat; whole-sequence; adjacent distinct runs; case normalization; three invariant property tests. The Wikipedia verbatim worked example `ATTCGATTCGATTCG` is locked (M14, `:296`). The lone range assertion is a performance bound (`< 30 s`), not a correctness claim — legitimate. No Greater/AtLeast/Contains weakening on correctness checks. Fuzz coverage exists (`Fuzzing/GenomicTandemFuzzTests.cs`, added in `ebc32e9f`).

### Findings / defects

None. No production-code defect. One documentation observation (not a defect): the TestSpec header says `GetTandemRepeatSummary` "wraps FindMicrosatellites" — accurate; the summary aggregates STR results, not `FindTandemRepeats` output. The "approximate" word elision in some prose is harmless given the deliberately documented exact-only simplification.

## Verdict & follow-ups

- **Stage A: PASS.** **Stage B: PASS.** **End-state: CLEAN.**
- No code changed. Canonical fixture (24 tests) passes; build clean (0 warnings).
- Independent cross-check coords confirmed: `ATTCGATTCGATTCG` → unit `ATTCG`, 0-based start 0, copy number 3, total length 15 (Wikipedia worked example); `ATGATGCCATGATGATG` → `ATG`×2 @0 + `ATG`×3 @8; minimum copy threshold k ≥ 2 enforced (Benson 1999).
