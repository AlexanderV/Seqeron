# Validation Report: PAT-APPROX-001 — Approximate Matching (Hamming Distance)

- **Validated:** 2026-06-12   **Area:** Pattern Matching
- **Canonical method(s):**
  - `ApproximateMatcher.HammingDistance(string, string)` — `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/ApproximateMatcher.cs:167`
  - `ApproximateMatcher.FindWithMismatches(string, string, int[, CancellationToken])` — `ApproximateMatcher.cs:21,35`
  - `ApproximateMatcher.FindWithMismatches(DnaSequence, ...)` overloads — `ApproximateMatcher.cs:87,96`
  - `SequenceExtensions.HammingDistance(ReadOnlySpan<char>, ReadOnlySpan<char>)` — `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs:264`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: Hamming distance** (fetched 2026-06-12). Exact wording: *"The Hamming distance between two equal-length strings of symbols is the number of positions at which the corresponding symbols are different."* Confirms the **equal-length precondition** and that the distance is a **metric** (non-negativity; identity d=0 iff strings identical; symmetry; triangle inequality). Examples confirmed: karolin/kathrin = 3, karolin/kerstin = 3, kathrin/kerstin = 4, 0000/1111 = 4.
- **Rosalind HAMM** (fetched 2026-06-12). Defines d_H(s,t) as the number of corresponding symbols that differ between two equal-length strings. Sample dataset `GAGCCTACTAACGGGAT` / `CATCGTAATGACGGCCT` → **7**.
- **Navarro (2001), Gusfield (1997)** (reference, per spec). Pattern matching with k mismatches = at each text offset the Hamming distance of the length-m window to the pattern must be ≤ k (no indels). Brute force O(n·m).

### Formula check
d_H(s,t) = Σ_{i=0}^{n-1} 1[s[i] ≠ t[i]] for equal-length s,t. Matches the cited Wikipedia equation in the spec (TestSpec §1.2). Approximate match at offset i ⇔ d_H(text[i..i+m], pattern) ≤ d, with 0-based offsets in [0, n−m]. Confirmed against Navarro definition.

### Edge-case semantics check
- Unequal lengths for Hamming distance → undefined per source; implementation throws `ArgumentException`. Correct per Wikipedia constraint.
- Identical → 0; completely different → length; empty/empty → 0. All match the mathematical definition.
- d=0 → exact matching only; d ≥ m → every window matches. Both follow directly from the definition.

### Independent cross-check (hand computation)
Rosalind HAMM, position-by-position (0-based):
`G≠C(1) A=A G≠T(2) C=C C≠G(3) T=T A=A C≠A(4) T=T A≠G(5) A=A C=C G=G G=G G≠C(6) A≠C(7) T=T` → **7**. Matches Rosalind sample output.

Small "≤1 mismatch" check: pattern `ACGG` in `ACGTACGT`. Windows: offset 0 `ACGT` (d=1, mismatch at pattern idx 3), offset 4 `ACGT` (d=1) → both match at d≤1; offsets 1,2,3 have d≥2. Matches code output and spec M14.

### Findings / divergences
None. Description is mathematically and biologically correct and matches all authoritative sources.

## Stage B — Implementation

### Code path reviewed
- `ApproximateMatcher.HammingDistance` (`ApproximateMatcher.cs:167-182`): null guards → `ArgumentNullException`; length mismatch → `ArgumentException("…equal length…")`; counts case-insensitive differing positions via `char.ToUpperInvariant`. Realises the formula exactly.
- `SequenceExtensions.HammingDistance` span overload (`SequenceExtensions.cs:264-277`): length guard → `ArgumentException`; same case-insensitive count. Consistent with string overload (spans cannot be null, so no null guard needed).
- `FindWithMismatches` (`ApproximateMatcher.cs:35-82`): empty/null seq or pattern → empty; negative `maxMismatches` → `ArgumentOutOfRangeException`; uppercases both; pattern longer than seq → empty; 0-based loop `i ∈ [0, n−m]`; per-window substitution count with early break once exceeded; yields `(Position=i, MatchedSequence, Distance=mismatches, MismatchPositions [pattern-relative], Substitution)` when `mismatches ≤ maxMismatches`. Correct.

### Formula realised correctly?
Yes. Hamming count is a direct per-position comparison; FindWithMismatches gates each length-m window by Hamming distance ≤ d with no indels, returning 0-based offsets. Matches the validated Navarro/Wikipedia definition.

### Cross-verification table recomputed vs code (tests executed)
| Case | Input | Expected (source) | Code result |
|------|-------|-------------------|-------------|
| Rosalind HAMM | GAGCCTACTAACGGGAT / CATCGTAATGACGGCCT | 7 | 7 ✓ |
| Wikipedia | karolin / kathrin | 3 | 3 ✓ |
| Wikipedia | karolin / kerstin | 3 | 3 ✓ |
| Wikipedia | kathrin / kerstin | 4 | 4 ✓ |
| Identical | ACGT / ACGT | 0 | 0 ✓ |
| All diff | AAAA / TTTT | 4 | 4 ✓ |
| Unequal len | ACGT / ACG | throw | ArgumentException ✓ |
| d=0 exact | ACGTACGT, ACGT, 0 | [0,4] | [0,4] ✓ |
| d=1 | ACGTACGT, ACGG, 1 | offsets 0,4 (mismatch idx 3) | matches ✓ |
| too many | ACGT, TGCA, 2 | [] (d=4>2) | [] ✓ |
| d≥len | XXXX, AB, 2 | [0,1,2] | [0,1,2] ✓ |

### Variant/delegate consistency
DnaSequence overloads delegate to the string implementation (`ApproximateMatcher.cs:90,102`); span and string Hamming overloads agree (test `HammingDistance_SpanApi_MatchesStringApi`). Verified.

### Test quality audit
`tests/Seqeron/Seqeron.Genomics.Tests/ApproximateMatcher_HammingDistance_Tests.cs` (55 tests in the `~Hamming` filter) asserts exact sourced values (Rosalind 7, Wikipedia 3/3/4), the equal-length exception, null/negative contracts, 0-based positions, INV-7/8/9 invariants, MismatchPositions correctness, and symmetry/triangle-inequality. Assertions are concrete and deterministic; edge cases from Stage A are covered.

### Findings / defects
None.

## Verdict & follow-ups
- Stage A: **PASS** — description matches Wikipedia, Rosalind, and Navarro/Gusfield exactly; worked examples reproduce (Rosalind = 7; Wikipedia 3/3/4).
- Stage B: **PASS** — implementation faithfully realises the formula; all edge cases enforced; cross-check values reproduce against the code.
- **State: CLEAN.** No defects. Hamming tests: 55 passed. Full suite: 4461 passed, 0 failed (baseline preserved). No code changes required.
