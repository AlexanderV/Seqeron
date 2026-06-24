# Validation Report: REP-DIRECT-001 — Direct Repeat Detection

- **Validated:** 2026-06-24   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindDirectRepeats(DnaSequence, minLength=5, maxLength=50, minSpacing=1)` + `string` overload — `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs:369-451` (core: `FindDirectRepeatsCore` at :416-451)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Direct repeat"** (accessed 2026-06-24): "a direct repeat occurs when a sequence is repeated with the same pattern downstream. There is no inversion and no reverse complement associated with a direct repeat." Direct repeat = identical nucleotide sequences appearing multiple times in the **same orientation**, possibly with intervening nucleotides. Tandem = repeated copies that lie directly adjacent (zero spacer) and can be direct or inverted; thus *tandem* describes positioning, *direct* describes orientation.
- **Wikipedia — "Repeated sequence (DNA)"** (accessed 2026-06-24): "Direct repeats occur when a nucleotide sequence is repeated with the same directionality" (e.g. CATCAT→CATCAT); "Inverted repeats occur when a nucleotide sequence is repeated in the inverse direction" (reverse complement); tandem = "directly adjacent"; interspersed = same/similar sequence at non-adjacent locations. Confirms the orientation distinction: direct ≠ inverted, and direct with spacer > 0 = interspersed.
- **Ussery et al. (2009)** and **Richard (2021) PMC8145212** (cited in spec) — consistent with same-strand same-orientation definition and supply the disease-relevant CAG context (S5).

### Definition / conventions confirmed
- Direct repeat = the **same** subsequence occurring ≥ 2× on the **same strand**, **same orientation** (literal match, NOT reverse complement), separated by a spacer; spacer 0 = tandem direct repeat.
- Parameters: min/max repeat-unit length, min spacer. Exact (perfect) matching — no mismatches allowed (spec lists no mismatch tolerance; the code matches literally).
- Per pair: `FirstPosition`, `SecondPosition`, `RepeatSequence`, `Length`, `Spacing = SecondPosition − FirstPosition − Length`.
- Coordinate base: **0-based**, matching the rest of the library and test expectations.

### Worked example (hand computed, independent)
Designed direct repeat with a non-palindromic motif so orientation is unambiguous:
`ATCGGG` + `NNNN`(4bp spacer) + `ATCGGG` = `ATCGGGNNNNATCGGG`. revcomp(`ATCGGG`) = `CCCGAT` ≠ `ATCGGG`, so a same-orientation literal match is a true *direct* repeat (not inverted). Copy1 at i=0, copy2 at j=10, len=6 → Spacing = 10 − 0 − 6 = **4**. Matches the `Spacing = j − i − len` formula.

### Edge-case semantics
No repeat → empty (M3); empty input → empty (M4); spacer 0 = tandem (M2); minSpacing filter excludes sub-threshold pairs (M13); short seq (< 2·minLength) → empty (M14); ≥3 copies → all pairwise pairs (S1). All defined and sourced.

### Findings / divergences
None. Description is biologically correct and consistent with authoritative sources; direct vs tandem vs inverted distinctions are exactly as the spec states.

## Stage B — Implementation

### Code path reviewed
`FindDirectRepeatsCore` (`RepeatFinder.cs:416-451`): builds a suffix tree over the sequence; for each `len = minLength..maxLength` and each start `i ≤ seq.Length − 2·len − minSpacing`, extracts `repeat = seq[i..i+len]`, finds all literal occurrences `j` via `suffixTree.FindAllOccurrences(repeat)`, keeps `j` with `j > i + len − 1 + minSpacing` (⇔ `Spacing = j − i − len ≥ minSpacing`), dedups on `(i, j, len)`, and emits `DirectRepeatResult`.

### Formula realised correctly?
- **Same orientation, literal:** matches the literal substring — NOT a reverse complement. Distinct from `FindInvertedRepeats` (`:320` uses `DnaSequence.GetReverseComplementString`) and `FindPalindromes` (`:597`). Confirmed this is a true direct repeat.
- **min/max unit length:** outer loop `len = minLength..maxLength` enforces both bounds; validation throws for `minLength < 2` and `maxLength < minLength` (both DnaSequence and string overloads).
- **min spacer:** filter `p > i + len − 1 + minSpacing` ⇔ `Spacing ≥ minSpacing`; `minSpacing = 0` admits tandem pairs; `minSpacing ≥ 1` excludes self/overlapping occurrences.
- **Coordinates:** 0-based; `Spacing = j − i − len`.
- **Outer bound** `i ≤ seq.Length − 2·len − minSpacing` is exactly the largest `i` admitting a valid second copy; drops no reachable repeat.

### Cross-verification table recomputed vs code (21 tests run, all pass)
| Case | Input | Params | Expected (hand) | Code |
|------|-------|--------|-----------------|------|
| M1 | ACGTATTTTACGTA | 5,10,1 | (0,9) len5 spacing4 | match |
| M2 | ACGTAACGTA | 5,10,0 | (0,5) spacing0 | match |
| M13 | ACGTAACGTA | 5,10,1 | empty (tandem filtered) | match |
| S1 | ACGTATTACGTATTACGTA | 5,5,1 | (0,7),(0,14),(7,14) | match |
| S4 | CCCGGGCCC+20bp+CCCGGGCCC | 9,9,1 | 1 pair, spacing20 | match |
| C1 | AAAAAATTTTAAAAAA | 4,6,1 | 9+4+1 = 14 pairs | match |
| (own) | ATCGGGNNNNATCGGG | 6,6,1 | (0,10) len6 spacing4 | match (formula) |

C1 recomputed: len4 {0,1,2}×{10,11,12}=9 (min spacing 10−2−4=4 ≥ 1); len5 {0,1}×{10,11}=4; len6 {0}×{10}=1; total 14. Agrees with code/test.

### Variant/delegate consistency
String overload normalises via `ToUpperInvariant()` then calls the same core; validation hoisted into an eager wrapper so exceptions surface at call time. S2 asserts identical results vs the DnaSequence overload; S3 confirms case-insensitivity. Consistent.

### Test quality audit
21 canonical tests (14 MUST, 5 SHOULD, 2 COULD). Assertions check exact sourced positions, lengths, spacings, and pair counts (not mere no-throw). Edge cases — empty, no-repeat, tandem, min/max length, min spacing, multi-copy, overlap, large — all covered. M5–M7 lock parameter validation.

### Findings / defects
None. Code faithfully realises the validated definition; direct repeats are distinguished from inverted (revcomp) and tandem (spacer) correctly.

## Verdict & follow-ups
- Stage A: PASS. Stage B: PASS. **State: CLEAN** — no defects.
- 21 DirectRepeat tests pass. No code or test changes required.
