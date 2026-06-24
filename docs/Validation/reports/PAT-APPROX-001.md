# Validation Report: PAT-APPROX-001 — Approximate Matching (Hamming Distance / k-mismatches)

- **Validated:** 2026-06-24   **Area:** Pattern Matching
- **Canonical method(s):**
  - `ApproximateMatcher.HammingDistance(string, string)` — `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs:172`
  - `ApproximateMatcher.FindWithMismatches(string, string, int[, CancellationToken])` — `ApproximateMatcher.cs:26,40`
  - `ApproximateMatcher.FindWithMismatches(DnaSequence, …)` overloads — `ApproximateMatcher.cs:92,101`
  - `SequenceExtensions.HammingDistance(ReadOnlySpan<char>, ReadOnlySpan<char>)` — `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:392`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Distance model
Hamming distance (substitutions only, **no indels**). Approximate occurrence at 0-based text
offset `i` ⇔ the equal-length window `text[i..i+m]` has Hamming distance ≤ `k` from the pattern.
Returns `ApproximateMatchResult(Position [0-based, in [0, n−m]], MatchedSequence, Distance,
MismatchPositions [pattern-relative 0-based indices], MismatchType.Substitution)`.

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Hamming distance** (fetched 2026-06-24). Verbatim: *"The Hamming distance between
  two equal-length strings of symbols is the number of positions at which the corresponding symbols
  are different."* Confirms the **equal-length precondition** and that it is a **metric**
  (non-negativity; identity d=0 iff identical; symmetry; triangle inequality). Worked examples
  reproduced: karolin/kathrin = 3, karolin/kerstin = 3, kathrin/kerstin = 4, 0000/1111 = 4,
  2173896/2233796 = 3.
- **Rosalind HAMM** (fetched 2026-06-24). Defines d_H(s,t) = number of corresponding symbols that
  differ between two equal-length strings. Sample `GAGCCTACTAACGGGAT` / `CATCGTAATGACGGCCT` → **7**.
- **Navarro (2001), Gusfield (1997)** (reference per spec). Pattern matching with k mismatches:
  at each text offset the Hamming distance of the length-m window to the pattern must be ≤ k,
  no insertions/deletions. Brute force O(n·m).

### Formula check
`d_H(s,t) = Σ_{i=0}^{n−1} 1[s[i] ≠ t[i]]` for equal-length s,t — matches TestSpec §1.2 and the
Wikipedia equation exactly. k-mismatch occurrence definition matches Navarro (no indels, 0-based
offsets in [0, n−m]).

### Edge-case semantics check
- Unequal lengths → undefined per source; impl throws `ArgumentException`. Correct.
- Identical → 0; completely different → length; empty/empty → 0. Match the definition.
- k=0 ⇒ exact matching; k ≥ m ⇒ every window matches. Both follow directly from the definition.
- Case-insensitivity (A1) is a documented implementation choice (DNA convention), not from the
  pure metric; standard for this domain.

### Independent cross-check (hand computation)
Rosalind HAMM, position-by-position (0-based):
`G≠C(1) A=A G≠T(2) C=C C≠G(3) T=T A=A C≠A(4) T=T A≠G(5) A=A C=C G=G G=G G≠C(6) A≠C(7) T=T` → **7** ✓.
k=1 check: pattern `ACGG` in `ACGTACGT` → offsets 0 and 4 are `ACGT`, d=1 (mismatch at pattern idx 3);
offsets 1,2,3 have d ≥ 2 → match set {0,4} ✓.

### Findings / divergences
None. Description matches Wikipedia, Rosalind, and Navarro/Gusfield exactly.

## Stage B — Implementation

### Code path reviewed
- `HammingDistance` (`ApproximateMatcher.cs:172-187`): null guards → `ArgumentNullException`;
  length mismatch → `ArgumentException("…equal length…")`; case-insensitive per-position count via
  `char.ToUpperInvariant`. Realises the formula exactly.
- `SequenceExtensions.HammingDistance` span (`SequenceExtensions.cs:392-405`): length guard →
  `ArgumentException`; same case-insensitive count. Consistent with the string overload (spans can't
  be null so no null guard).
- `FindWithMismatches` (`ApproximateMatcher.cs:40-87`): null/empty seq or pattern → empty; negative
  `maxMismatches` → `ArgumentOutOfRangeException`; uppercases both; pattern longer than seq → empty;
  0-based loop `i ∈ [0, n−m]`; per-window substitution count with early break once `> maxMismatches`;
  yields `(Position=i, MatchedSequence, Distance=mismatches, MismatchPositions [pattern-relative],
  Substitution)` when `mismatches ≤ maxMismatches`. Cancellation checked every 1000 offsets. Correct.

### Formula realised correctly?
Yes. Direct per-position comparison; each length-m window gated by Hamming ≤ k with no indels;
0-based offsets. Matches the validated Navarro/Wikipedia definition.

### Cross-verification table recomputed vs code (tests executed)
| Case | Input | Expected (source) | Code |
|------|-------|-------------------|------|
| Rosalind HAMM | GAGCCTACTAACGGGAT / CATCGTAATGACGGCCT | 7 | 7 ✓ |
| Wikipedia | karolin / kathrin | 3 | 3 ✓ |
| Wikipedia | karolin / kerstin | 3 | 3 ✓ |
| Wikipedia | kathrin / kerstin | 4 | 4 ✓ |
| Identical | ACGT / ACGT | 0 | 0 ✓ |
| All diff | AAAA / TTTT | 4 | 4 ✓ |
| Unequal len | ACGT / ACG | throw | ArgumentException ✓ |
| k=0 exact | ACGTACGT, ACGT, 0 | [0,4] | [0,4] ✓ |
| k=1 | ACGTACGT, ACGG, 1 | offsets 0,4 (mismatch idx 3) | matches ✓ |
| too many | ACGT, TGCA, 2 | [] (d=4>2) | [] ✓ |
| k≥m | XXXX, AB, 2 | [0,1,2] | [0,1,2] ✓ |
| MismatchPositions | AXGX vs ACGT, k=2 | {1,3} | {1,3} ✓ |

### Variant/delegate consistency
DnaSequence overloads delegate to the string implementation (`ApproximateMatcher.cs:95,107`); span and
string Hamming overloads agree (`HammingDistance_SpanApi_MatchesStringApi`). Verified.

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/ApproximateMatcher_HammingDistance_Tests.cs` (38 tests) asserts
exact sourced values (Rosalind 7; Wikipedia 3/3/4), the equal-length exception, null/negative
contracts, 0-based positions, INV-7/8/9, MismatchPositions correctness, symmetry, triangle inequality,
span-API consistency, DnaSequence overload equivalence, and cancellation. Assertions are concrete and
deterministic; all Stage-A edge cases covered.

### Findings / defects
None.

## Verdict & follow-ups
- Stage A: **PASS** — description matches Wikipedia, Rosalind, and Navarro/Gusfield exactly; worked
  examples reproduce (Rosalind = 7; Wikipedia 3/3/4).
- Stage B: **PASS** — implementation faithfully realises the Hamming/k-mismatch model; all edge cases
  enforced; cross-check values reproduce against the code.
- **State: CLEAN.** No defects. Hamming tests: 38 passed, 0 failed. No code changes required (full
  suite not re-run since no source was touched).
