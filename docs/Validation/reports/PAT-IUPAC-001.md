# Validation Report: PAT-IUPAC-001 — IUPAC Degenerate Motif Matching

- **Validated:** 2026-06-24   **Area:** Pattern Matching
- **Canonical method(s):** `MotifFinder.FindDegenerateMotif(DnaSequence, string)` (+ cancellable / string variants, `FindDegenerateMotifCore`); `IupacHelper.MatchesIupac(char, char)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Nucleic acid notation** (https://en.wikipedia.org/wiki/Nucleic_acid_notation), fetched live this session. Table 1 returns the full ambiguity code with the exact bases and mnemonics:
  R=A,G (purine); Y=C,T (pyrimidine); S=C,G (strong); W=A,T (weak); K=G,T (keto); M=A,C (amino); B=C,G,T (not A); D=A,G,T (not C); H=A,C,T (not G); V=A,C,G (not T); N=A,C,G,T (any).
- This matches the NC-IUB (1984) / IUPAC-IUB (1970) standard cited in the TestSpec, with zero divergence.

### Matching model / direction (asymmetry)
The model is **asymmetric and pattern-degenerate**: the *pattern* carries IUPAC ambiguity codes; the *sequence* is concrete A/C/G/T (`DnaSequence` validates to ACGT). At each position the concrete sequence base must be a **member** of the set the pattern code represents. The sequence side is literal — a text base does NOT match a pattern code unless it is in that code's set, and a literal pattern base (A/C/G/T) matches only itself. This is the standard direction used by EMBOSS fuzznuc / Biopython ambiguity search. INV-9 ("symmetric for standard bases") is the only symmetry claimed and is trivially true; no over-claim.

### Positions / conventions
0-based start positions; result range [0, len−pattern]; case-insensitive (pattern uppercased, `DnaSequence` already normalized). Edge cases all defined & sourced: empty pattern/sequence and pattern-longer-than-sequence → no matches; null sequence → ArgumentNullException; non-IUPAC pattern char → ArgumentException; U/gap not supported (DNA-only).

### Independent cross-checks (hand-computed)
- **Table membership for the negation-derived codes** (the focus of this validation): B={C,G,T} excludes A; D={A,G,T} excludes C; H={A,C,T} excludes G; V={A,C,G} excludes T; N={A,C,G,T}. All confirmed against Wikipedia.
- **Pattern `RYSW`** (R={A,G},Y={C,T},S={C,G},W={A,T}) vs text:
  - `AGCT`: R/A✓ Y/G✗ → no match (asymmetry: G is not in Y's set).
  - `ACGA`: R/A✓ Y/C✓ S/G✓ W/A✓ → match at 0.
  - `GTCT`: R/G✓ Y/T✓ S/C✓ W/T✓ → match at 0.
- M50 `R`/`ATGC` → [0,2]; M51 `Y`/`ATGC` → [1,3]; M53 `RTG`/`ATGCGTGC` → [0,4]; M54 `CANNTG`/`CAGCTG` → [0]. All confirmed.

### Findings
None. Stage A PASS.

## Stage B — Implementation

### Code paths reviewed
- `IupacHelper.cs:16-37` — `MatchesIupac`: 15-arm switch, positive-set membership for every code incl. B/D/H/V/N (lines 22, 29-32). Unknown code → ArgumentOutOfRangeException.
- `MotifFinder.cs:46-63` — `IupacCodes` dictionary: B="CGT", D="AGT", H="ACT", V="ACG", N="ACGT" — verbatim-correct.
- `MotifFinder.cs:90-119` — canonical `FindDegenerateMotif(DnaSequence, string)`: null-check, empty short-circuit, `ValidateIupacPattern`, O(n×m) scan, membership via `IupacCodes[motifChar].Contains(seqChar)`.
- `MotifFinder.cs:148-201` — `FindDegenerateMotifCore` (string + cancellable variants): same scan via switch. **The B/D/H/V/N arms are now positive-set membership** (lines 181-185: `seqChar is 'C' or 'G' or 'T'`, etc.), NOT the old negation form. This is the key item flagged for this validation, and it is correct.

### Verification of all 15 codes incl. negation-derived ones
The three encodings of the table — `IupacHelper.MatchesIupac` switch, `IupacCodes` dictionary, and `FindDegenerateMotifCore` switch — agree exactly for all 15 codes. Critically, B/D/H/V/N now use explicit positive membership in all three, so a non-ACGT sequence character correctly fails to match (the prior negation bug, where `B => seqChar != 'A'` would have spuriously matched `X`/`-`/`N`, is gone).

### Cross-verification table recomputed vs code
| Test | Input | Expected | Code result |
|------|-------|----------|-------------|
| M50 | R / ATGC | [0,2] | ✅ |
| M51 | Y / ATGC | [1,3] | ✅ |
| M53 | RTG / ATGCGTGC | [0,4] | ✅ |
| M54 | CANNTG / CAGCTG | [0] | ✅ |
| S8  | BBB / CGT | [0] | ✅ |
| S9  | DDD / AGT | [0] | ✅ |
| S10 | HHH / ACT | [0] | ✅ |
| S11 | VVV / ACG | [0] | ✅ |
| H10 | B/D/H/V/N vs "X" | empty | ✅ |
| H12 | B vs "C-GT" | [0,2,3] (gap excluded) | ✅ |

### Variant/delegate consistency
Dictionary path (canonical DnaSequence overload) and the Core switch (string/cancellable) produce identical results for all ACGT sequences; for non-ACGT input the Core path now also correctly rejects, matching the spirit of the validated table.

### Test quality audit
`IupacMotifMatchingTests.cs` asserts exact positions, matched substrings, and exception types — deterministic, no tautologies. Hardening tests H1–H14 lock the positive-set behaviour of B/D/H/V/N through the string overload, including the bug-regression cases (H10: `X` not matched; H11: sequence-`N` not matched; H12: gap excluded) plus positive regressions (H13 GAAY, H14 N). 245 tests under the `~Iupac` filter pass.

### Findings
None. Stage B PASS.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS**, **State: CLEAN.**
- Degeneracy table verified live against Wikipedia; all 15 codes including the negation-derived B/D/H/V/N are implemented as positive-set membership across all three table encodings and agree. Asymmetric pattern-degenerate match direction confirmed by hand against `RYSW`.
- No code changes this session (prior fix already in place and correct). Iupac-filtered tests: 245 passed, 0 failed. Build clean.
