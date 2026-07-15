---
type: concept
title: "Edit distance (Levenshtein) and indel-tolerant approximate search"
tags: [alignment, algorithm]
mcp_tools:
  - edit_distance
  - find_with_edits
sources:
  - docs/algorithms/Pattern_Matching/Edit_Distance.md
source_commit: 1ab8b7c126cdaca6ce7f9cf835a65ffbf4997441
created: 2026-07-15
updated: 2026-07-15
graph:
  relationships:
    - predicate: alternative_to
      object: concept:approximate-pattern-matching-mismatches
      source: edit-distance
      evidence: "Edit distance models insertion, deletion AND substitution events and 'provides a more appropriate notion of distance than Hamming distance when indels are possible' (Edit_Distance.md §2.1); the Hamming PAT-APPROX family is substitution-only. Same problem — approximate matching — solved with different edit models."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: edit-distance
      evidence: "Test Unit ID: PAT-APPROX-002 (Edit_Distance.md header table), the area-prefixed unit under which Levenshtein distance is validated."
      confidence: high
      status: current
---

# Edit distance (Levenshtein) and indel-tolerant approximate search

**Edit (Levenshtein) distance** is the minimum number of single-character
**insertions, deletions, and substitutions** — each cost 1 — needed to transform one
string into another. Because it allows indels, it is the *indel-tolerant* member of the
approximate-matching family: it captures the difference between two sequences of
**different lengths**, which the fixed-length
[[approximate-pattern-matching-mismatches|Hamming-distance model]] cannot. In
bioinformatics it models the insertion/deletion/substitution events that separate real
sequences, so it is "a more appropriate notion of distance than Hamming distance when
indels are possible" (Levenshtein 1966; Wagner & Fischer 1974; Navarro 2001). Validated
under test unit **PAT-APPROX-002**; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern.

## The recurrence (Wagner–Fischer DP)

For strings `a`, `b`, `lev(a,b)` is defined recursively: an empty operand forces `|a|` or
`|b|` edits; matching heads recurse on the tails at no cost; otherwise take `1 +` the
minimum of the three edit branches:

```
lev(a,b) =
  |a|                                        if |b| = 0
  |b|                                        if |a| = 0
  lev(tail(a), tail(b))                      if head(a) = head(b)
  1 + min( lev(tail(a), b),                  (deletion)
           lev(a, tail(b)),                  (insertion)
           lev(tail(a), tail(b)) )           (substitution)   otherwise
```

Seqeron computes this with a **space-optimized two-row** dynamic program: only the
previous and current DP rows are kept, so distance for lengths `m`, `n` costs
**O(m·n) time, O(n) space** rather than the full O(m·n) matrix. The returned value is the
last cell of the final row. This is the same DP fabric as
[[global-alignment-needleman-wunsch]] global alignment run with unit costs and a *minimize
cost* objective instead of *maximize similarity score* — edit distance is the metric dual
of a unit-penalty NW alignment — but the core routine returns only the scalar distance, no
traceback (see *No alignment traceback* below).

## Invariants

| ID | Invariant | Why |
|----|-----------|-----|
| INV-01 | `d(a,b) = 0` iff the strings are identical (under the comparison `EditDistance` uses) | Zero cost only when every aligned char matches and lengths agree |
| INV-02 | `d(a,b) ≥ ∣len(a) − len(b)∣` | At least the length gap must be closed by insertions/deletions |
| INV-03 | `d(a,b) ≤ max(len(a), len(b))` | Deletions plus substitutions always suffice to convert one into the other |
| — | equal-length strings: `d ≤ Hamming(a,b)` | Substitutions alone realize the Hamming path; indels can only help |

The last row is the bridge to the Hamming model: on **equal-length** inputs Levenshtein
distance never exceeds the [[approximate-pattern-matching-mismatches|Hamming distance]],
and equals it when no indel shortens the edit script.

## Two surfaces: `EditDistance` and `FindWithEdits`

`ApproximateMatcher` (in `Seqeron.Genomics.Alignment`) exposes two entry points:

- **`EditDistance(string, string)`** — the two-row Levenshtein distance. **Case-sensitive**:
  it compares characters as-is, so callers needing normalized comparison must uppercase
  first. Throws `ArgumentNullException` on a null operand.
- **`FindWithEdits(sequence, pattern, maxEdits)`** — brute-force **approximate search**. It
  **uppercases** both sequence and pattern (case-insensitive, unlike the core distance
  method), then compares the pattern against every window whose length ranges from
  `pattern.Length − maxEdits` to `pattern.Length + maxEdits` — the variable-length window
  band is what lets a match absorb indels — and yields each window whose `EditDistance` to
  the pattern is `≤ maxEdits`. A `DnaSequence` overload is a thin typed wrapper.

Each hit reports `Position` (0-based window start), `MatchedSequence`, `Distance`
(observed edit distance), and a `MismatchType`: **`Substitution`** when the edit distance
equals the Hamming distance on an equal-length window, otherwise **`Edit`** (an indel was
involved).

**Complexity.** `EditDistance` is O(m·n) / O(n). `FindWithEdits` is
`O(s·(2e+1)·p·(p+e))` time, `O(p+e)` space, where `s` = sequence length, `p` = pattern
length, `e` = `maxEdits` — every one of the `2e+1` window lengths at every start invokes
`EditDistance` on a window of length `p−e … p+e`.

## Contract and corner cases

- `EditDistance` throws `ArgumentNullException` on a null string.
- `FindWithEdits(string, …)` returns **no matches** when the sequence or pattern is null or
  empty; it throws `ArgumentOutOfRangeException` when `maxEdits < 0`.
- The **`DnaSequence` overload has no null guard** — it dereferences `sequence.Sequence`,
  so a null `DnaSequence` throws (a `NullReferenceException`) rather than returning empty.
  This is the same sharp edge the Hamming
  [[approximate-pattern-matching-mismatches|`FindWithMismatches`]] typed wrappers carry.
- `"" vs "abc"` → distance `3`. `kitten → sitting` → distance `3` (substitute k→s,
  e→i, insert g).

## Deliberate simplifications

- **No alignment traceback.** `FindWithEdits` reports the edit distance and matched window
  but does not reconstruct insertion/deletion coordinates; for indel edits
  `MismatchPositions` is returned **empty**. Callers needing a full edit script / aligned
  columns must use the [[global-alignment-needleman-wunsch|global alignment]] family, which
  performs traceback.
- **Brute-force scan.** The variable-window search is easy to reason about but not tuned
  for very large texts or very permissive `maxEdits` — no banded/bit-parallel
  (Myers/Ukkonen) acceleration.
- **No transpositions.** Damerau-style adjacent-swap edits and other extended operations
  are **not** implemented; the model is the classic three-operation Levenshtein metric.

## Relation to the rest of the pattern-matching family

Edit distance is the **indel-tolerant `alternative_to`** the substitution-only
[[approximate-pattern-matching-mismatches|Hamming approximate matcher]]: both answer "does
this pattern approximately occur here?", but the Hamming model requires equal-length
windows (mismatches only) while edit distance admits length-changing indels. It shares the
`ApproximateMatcher.cs` home and the `ApproximateMatchResult` shape (`Position` /
`MatchedSequence` / `Distance` / `MismatchType`) with that Hamming surface. For a full
alignment (aligned columns, gap placement, affine gaps) rather than a scalar distance or a
match list, escalate to [[global-alignment-needleman-wunsch]] or its ends-free
[[semi-global-alignment-fitting]] variant.

## Reference sources

`docs/algorithms/Pattern_Matching/Edit_Distance.md` (PAT-APPROX-002, "Simplified" status).
Primary literature: **Levenshtein (1966)** (original metric), **Wagner & Fischer (1974)**
(the string-to-string correction DP), **Navarro (2001)** (approximate string matching
survey), **Berger, Waterman & Yu (2021)** (Levenshtein distance in biological database
search); plus Wikipedia (Levenshtein / Edit distance) and Rosetta Code. Oracle:
`kitten → sitting = 3`. No deviations from the cited theory; the search layer is an
intentional brute-force simplification without traceback.
