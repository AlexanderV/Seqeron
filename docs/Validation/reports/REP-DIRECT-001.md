# Validation Report: REP-DIRECT-001 — Direct Repeat Detection

- **Validated:** 2026-06-12   **Area:** Repeats
- **Canonical method(s):** `RepeatFinder.FindDirectRepeats(DnaSequence, minLength=5, maxLength=50, minSpacing=1)` and the `string` overload (`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs:366-430`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — "Direct repeat"** (https://en.wikipedia.org/wiki/Direct_repeat): "a direct repeat occurs when a sequence is repeated with the same pattern downstream. There is no inversion and no reverse complement associated with a direct repeat." Two or more copies of a specific sequence in the **same orientation** on the **same strand**, with possible intervening nucleotides. A **tandem** direct repeat = copies lying directly adjacent (zero spacer).
- **Wikipedia — "Repeated sequence (DNA)"** (https://en.wikipedia.org/wiki/Repeated_sequence_(DNA)): direct repeat = "a nucleotide sequence is repeated with the same directionality"; inverted repeat = "repeated in the inverse direction" (reverse complement); tandem = "directly adjacent". This confirms the orientation distinction: direct ≠ inverted/palindromic.
- **Ussery et al. (2009), Computing for Comparative Microbial Genomics (Springer, ch. 8)** — cited in spec as the technical reference for direct/repeat detection in microbial genomes; consistent with the same-strand same-orientation definition with a spacer between copies.
- **Richard (2021) PMC8145212** — trinucleotide (e.g. CAG) repeat context; supports the disease-relevant S5 test (CAGCAG copies).

### Definition / conventions confirmed
- Direct repeat = the **same** subsequence occurring twice (or more) on the **same strand** in the **same orientation** (NOT reverse complement — that would be an inverted repeat), separated by a spacer; spacer length 0 = tandem direct repeat.
- Parameters: minimum repeat-unit length, maximum repeat-unit length, minimum spacer/gap (perfect/exact matching here).
- Reported per pair: first copy position, second copy position, repeat sequence, length, spacer = `SecondPosition − FirstPosition − Length`.
- Coordinate base: **0-based** (matches the rest of the library and the test expectations).

### Worked example (hand computed)
"ATGCG…spacer…ATGCG": copy1 at i, copy2 at j → Spacing = j − i − 5. For `"ACGTATTTTACGTA"` (M1): "ACGTA" at i=0 and j=9 ⇒ Spacing = 9 − 0 − 5 = **4**. For tandem `"ACGTAACGTA"` (M2): j=5 ⇒ Spacing = 5 − 0 − 5 = **0**. Both match the spec expected values.

### Edge-case semantics
- No direct repeat → empty (M3); empty input → empty (M4); spacer 0 = adjacent/tandem (M2); minSpacing filter excludes pairs below threshold (M13); short sequence (< 2·minLength) → empty (M14); ≥3 copies → all pairwise pairs (S1). All defined and sourced.

### Findings / divergences
None. Description is biologically correct and consistent with authoritative sources.

## Stage B — Implementation

### Code path reviewed
`FindDirectRepeatsCore` (`RepeatFinder.cs:395-430`): builds a suffix tree, iterates lengths `minLength..maxLength`, for each start `i` extracts `repeat = seq[i..i+len]`, finds all occurrences `j` of that substring, keeps `j` with `j > i + len − 1 + minSpacing` (i.e. `Spacing = j − i − len ≥ minSpacing`), dedups on `(i,j,len)`.

### Validated against definition
- **Identity, same orientation:** matches the literal substring (`suffixTree.FindAllOccurrences(repeat)`) — NOT a reverse complement. Distinct from `FindInvertedRepeats` (uses `GetReverseComplementString`) and `FindPalindromes`. Correct: this is a true direct repeat, same strand, same orientation.
- **min/max unit length:** outer length loop `len = minLength..maxLength` enforces both bounds; validation throws for `minLength < 2` and `maxLength < minLength`.
- **min spacer:** filter `p > i + len − 1 + minSpacing` ⇔ `Spacing ≥ minSpacing`; `minSpacing=0` admits tandem pairs. Correct.
- **Coordinate base:** 0-based; `Spacing = j − i − len` set in `DirectRepeatResult`.
- **Outer-loop bound:** `i <= seq.Length − len*2 − minSpacing` is exactly the largest `i` that can admit a valid second copy (`i + 2·len + minSpacing ≤ seq.Length`); it does not drop any reachable repeat. Verified.

### Cross-verification table recomputed vs code
| Case | Input | Params | Expected (hand) | Code result |
|------|-------|--------|-----------------|-------------|
| M1 | ACGTATTTTACGTA | 5,10,1 | (0,9) len5 spacing4 | matches |
| M2 | ACGTAACGTA | 5,10,0 | (0,5) spacing0 | matches |
| M13 | ACGTAACGTA | 5,10,1 | empty (only tandem) | matches |
| S1 | ACGTATTACGTATTACGTA | 5,5,1 | (0,7),(0,14),(7,14) | matches |
| S4 | CCCGGGCCC+20bp+CCCGGGCCC | 9,9,1 | 1 pair, spacing 20 | matches |
| C1 | AAAAAATTTTAAAAAA | 4,6,1 | 9+4+1 = 14 pairs | matches |

C1 recomputed by hand: len4 → i∈{0,1,2} × second-block {10,11,12} filtered by `Spacing≥1` = 9; len5 → 4; len6 → 1; total 14. Agrees with code and test.

### Variant/delegate consistency
String overload normalizes via `ToUpperInvariant()` then calls the same core; S2 asserts identical results to the `DnaSequence` overload; S3 confirms case-insensitivity. Consistent.

### Test quality audit
21 canonical tests (14 MUST, 5 SHOULD, 2 COULD) plus property/snapshot/metamorphic coverage. Assertions check exact sourced positions, lengths, spacings and pair counts (not mere no-throw). Edge cases (empty, no-repeat, tandem, threshold filters, multi-copy, overlap, large) all covered.

### Findings / defects
None. Code faithfully realises the validated definition.

## Verdict & follow-ups
- Stage A: PASS. Stage B: PASS. State: CLEAN — no defects.
- Build + tests green: 32 DirectRepeat-filtered tests pass; full `Seqeron.Genomics.Tests` suite = 4461 passed, 0 failed.
- No code or test changes were required.
