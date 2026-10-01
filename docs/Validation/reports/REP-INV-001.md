# Validation Report: REP-INV-001 — Inverted Repeat Detection

- **Re-validated:** 2026-09-29 / 2026-09-30 (review campaign 2026-09, batch B04: fixes F10, F21–F24; supersedes the 2026-06-24 PASS/PASS report, whose "report every (i, j, arm) tuple" rule was wrong)   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindInvertedRepeats(DnaSequence | string, minArmLength = 4, maxLoopLength = 50, minLoopLength = 3, maxMismatches = 0, maxArmLength = int.MaxValue, allowWobble = false)`; scored variant `RepeatFinder.FindInvertedRepeatsScored(string | DnaSequence, gapPenalty = 12, threshold = 50, matchScore = 3, mismatchScore = −4, maxRepeatLength = 2000)`; MCP `find_inverted_repeats` (delegates, exposes all options).
- **Stage A verdict:** PASS-WITH-NOTES (reporting rule re-sourced: EMBOSS `palindrome`; options sourced from `palindrome.c` / `einverted.c`)
- **Stage B verdict:** FAIL → fixed (F10), options implemented (F21–F24)

## Stage A — Description

### Sources opened
- **EMBOSS `palindrome.c`** (M. Faller; raw GitHub kimrutherford/EMBOSS, 6.6.0): stems start at a complementary outer pair and are extended inward while `mismatches <= maxmismatches && ic < ir`; trailing mismatches trimmed (`count -= mismatchAtEnd`); kept when `count >= minpallen && gap <= gaplimit`; with `-overlap Y` a stem inside an earlier one in both halves (`palindrome_AInB`) is dropped; outer pairs limited to `rev <= current + 2·maxpallen + gaplimit`, and `palindrome_Print` skips stems longer than `maxpallen` (they remain in the list and suppress sub-stems). ACD: minpallen 10, maxpallen 100, gaplimit 100, nummismatches 0, both lengths clamped to len/2.
- **EMBOSS `einverted.c`** (Durbin & Thierry-Mieg 1993) and `einverted.acd` (6.6.0): outward-growing local alignment of the sequence with its reverse complement (match 3, mismatch −4, linear gap 12, threshold 50, maxrepeat 2000 — Durbin's comment mentions a 4000 compile-time value); only a/c/g/t score as a match; ring buffer of `maxrepeat` rows, per-row best (`localMax`), left-start table (`back`), deferred reporting with window clearing and trace-back order same-row gap → previous-row gap → diagonal.
- **Wikipedia — Inverted repeat / Stem-loop / Palindromic sequence** (2026-06 session): IR = arm + spacer ≥ 0 + reverse complement; loops < 3 nt sterically impossible (`CanFormHairpin`).
- **G·U wobble:** Crick 1966 (J Mol Biol 19:548–555); Varani & McClain 2000, EMBO Rep 1:18–23 (WebSearch: G·U "fundamental unit of RNA secondary structure", thermodynamic stability comparable to Watson–Crick, nearly isomorphic). Canonical pair set: `RnaSecondaryStructure.CanPair` (ViennaRNA default pairs A·U, G·C, G·U; T read as U).

### Conventions
0-based coordinates; `LoopLength = RightArmStart − (LeftArmStart + ArmLength)`; `TotalLength = 2·ArmLength + LoopLength`; only A/C/G/T pair by default (einverted convention; `palindrome` lets n pair n); `minLoopLength` is an extension (0 = `palindrome`). Scored results: 0-based inclusive arm coordinates, einverted report order.

## Stage B — Implementation

- **F10 (2026-09-29):** every sub-stem was reported (`GAATTCAAAAGAATTC` → 6 hits; EMBOSS 1); N/IUPAC paired; lazy validation → maximal non-nested stems, ACGT via canonical `GetComplementBase`, eager validation.
- **F21 `maxMismatches`:** per-anti-diagonal implementation of the `palindrome` candidate walk + exact nesting test (a container on diagonal D ± m is ≥ m longer, m ≤ loop − minLoop). `Mismatches` added to `InvertedRepeatResult` (init property).
- **F22 `FindInvertedRepeatsScored`:** einverted scan reimplemented (own code) incl. ring reuse and stale-`localMax` behaviour; `ScoredInvertedRepeatResult` record.
- **F23 `maxArmLength`:** `-maxpallen` rule (start-span bound + print filter; longer stems not split).
- **F24 `allowWobble`:** pairing via `RnaSecondaryStructure.CanPair`; maximal-stem rule unchanged.

### Cross-verification (0 mismatches everywhere)
| Check | Cases | Result |
|---|---|---|
| Default exact stems vs `palindrome` binary (minLoop 0) / Python brute force | 1000 / 3000 | 0 mismatches (F10) |
| Options vs `palindrome` 6.6.0 (k 0–6, ≈ half with `-maxpallen` < len/2, n 8–400) | 5100 cases, 144 727 stems (113 676 with mismatches) | 0 |
| Options vs literal `palindrome.c` transcription generalised to minLoop 0–6 / wobble / N / U / lowercase | 6000 | 0 |
| Exact stems vs all-stems brute force (Python, ≈ half with the RNA wobble pair set; + 300 wobble cases in the C# test) | 1653 | 0 |
| Scored vs `einverted` 6.6.0 binary (planted mismatched/gapped repeats, n 20–4000, gap 0–16, threshold 0–70, match 1–5, mismatch 0…−6, maxrepeat 2–2000) | 4700 (89 aborted by einverted with SIGFPE) — 11 100 repeats | 0 (coordinates, score, matches, mismatches, gaps, all three alignment rows) |
| Scored vs einverted built from 6.6.0 source with only the `(100*nmatch)/(nmatch+nmis)` division guarded | 3400, incl. the aborting parameter sets; 8 560 repeats | 0 |

### Worked values (locked in tests)
- `palindrome -minpallen 4 -gaplimit 10`: `GAATTCAGGAAAACCTCAATTC` k=1 → (0,13,9, 1 mm); k=2 → (0,8,6)(0,10,5)(0,12,5)(0,13,9)(12,17,4). `GGGGGGGGGGAAACCCCCCCCCC` → (0,13,10); `-maxpallen` 4/5/6 → nothing. `TTGCATGCAAAAAATTTTTTTGCATGCAA` unbounded → (0,5,5)(0,15,14)(8,14,6)(19,24,5); maxpallen 6 → drops (0,15,14); 5 → (0,5,5)(19,24,5).
- `einverted` defaults on `CCCAACCCATGCGTACGTTAGCCTAGGATCCATTTTTTTTTGGATACTAGGCAACGTACGCATGGGAAGGG` → "Score 60: 28/31 (90%) matches, 1 gaps", 1..32 / 71..41 (= 28·3 − 3·4 − 12).

### Tests
Unit: `RepeatFinder_InvertedRepeat_Tests` (F10) + `RepeatFinder_InvertedRepeatOptions_Tests` (options, scored, literal oracle, brute force); heavy tier: `Properties/RepInvOptionsProperties`, `Metamorphic/RepInvOptionsMetamorphicTests` (revcomp mirror, wobble reverse mirror, scored case / N-prefix shift), `Fuzzing/RepInvOptionsFuzzTests`; MCP `FindInvertedRepeatsTests`, `FindInvertedRepeatsScoredTests` (tool `find_inverted_repeats_scored` → `FindInvertedRepeatsScored`, B04 F49).

## Verdict
**FIXED + extended.** Exact default output unchanged by the options; `palindrome` and `einverted` behaviour reproduced with 0 mismatches. Evidence: [REP-INV-001-Evidence.md](../../Evidence/REP-INV-001-Evidence.md).
