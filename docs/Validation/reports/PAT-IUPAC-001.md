# Validation Report: PAT-IUPAC-001 — IUPAC Degenerate Motif Matching

- **Validated:** 2026-06-12   **Area:** Pattern Matching
- **Canonical method(s):** `MotifFinder.FindDegenerateMotif(DnaSequence, string)` (+ cancellable / string variants); `IupacHelper.MatchesIupac(char, char)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Bioinformatics.org — IUPAC codes** (https://www.bioinformatics.org/sms/iupac.html): full ambiguity table fetched and read.
- **Wikipedia — Nucleic acid notation** (https://en.wikipedia.org/wiki/Nucleic_acid_notation): degenerate symbol table fetched; explicitly attributes the codes to **NC-IUB (1984)** with the foundational IUPAC notation dating to **1970** (Biochemistry 9(20):4022–4027). This matches the spec's cited references.

### Degeneracy table check (verified EXACTLY against both sources)
| Code | Represents | Bioinformatics.org | Wikipedia |
|------|------------|--------------------|-----------|
| R | A, G | ✅ "A or G" | ✅ |
| Y | C, T | ✅ "C or T" | ✅ |
| S | G, C | ✅ "G or C" | ✅ |
| W | A, T | ✅ "A or T" | ✅ |
| K | G, T | ✅ "G or T" | ✅ |
| M | A, C | ✅ "A or C" | ✅ |
| B | C, G, T | ✅ "not A" | ✅ |
| D | A, G, T | ✅ "not C" | ✅ |
| H | A, C, T | ✅ "not G" | ✅ |
| V | A, C, G | ✅ "not T" | ✅ |
| N | A, C, G, T | ✅ "any base" | ✅ |

All 11 ambiguity codes + 4 standard bases confirmed; zero divergence.

### Match-direction (CRITICAL semantic question)
The spec defines matching **asymmetrically**: the **pattern** carries IUPAC ambiguity codes, and at each position the concrete sequence base (always one of A/C/G/T, since `DnaSequence` validates input) must be a member of the set the pattern code represents. The sequence itself is *not* degenerate. This is the standard direction used by reference tools for motif/restriction-site search (EMBOSS `fuzznuc`, Biopython ambiguity search, REBASE site lookup), and is correct and defensible. INV-9 ("symmetric for standard bases") only claims symmetry for the literal A/C/G/T cases, which is trivially true; it does not claim full symmetry, so there is no over-claim.

### Edge-case semantics
Defined and sourced: standard base = exact match; N = any; empty pattern / empty sequence / pattern-longer-than-sequence = no matches; null sequence = `ArgumentNullException`; non-IUPAC pattern char = `ArgumentException` (justified by IUPAC-IUB: only 15 codes defined). U/T equivalence is intentionally **not** supported (DNA-only system) — documented in spec §6.3, not a defect.

### Independent cross-check (hand-computed)
- `GAAY` vs `GAAC` → Y={C,T} ∋ C → **match**; vs `GAAT` → ∋ T → **match**; vs `GAAA` → A ∉ {C,T} → **no match**. ✅
- `N` matches A, C, G, T. ✅
- M50: `R` on `ATGC` → positions 0 (A) and 2 (G) = [0,2]. ✅
- M53: `RTG` on `ATGCGTGC` → [0]=ATG, [4]=GTG = [0,4]. ✅
- M54: `CANNTG` on `CAGCTG` → [0]. ✅

### Findings
None. Stage A PASS.

## Stage B — Implementation

### Code paths reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/IupacHelper.cs:16-37` — `MatchesIupac`: 15-arm switch, table verbatim-correct, unknown/lowercase code → `ArgumentOutOfRangeException`.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs:46-63` — `IupacCodes` dictionary, table verbatim-correct.
- `MotifFinder.cs:90-119` — canonical `FindDegenerateMotif(DnaSequence, string)`: null-check, empty-motif short-circuit, pattern validation, O(n×m) scan over `i ∈ [0, len-m]`, membership via `IupacCodes[motifChar].Contains(seqChar)`.
- `MotifFinder.cs:148-201` — `FindDegenerateMotifCore` (string + cancellable variants): same logic via a switch; B/D/H/V expressed as `seqChar != 'A'` etc. and `N => true`.

### Formula realised correctly?
Yes. The canonical `DnaSequence` overload looks up the exact base set and tests membership — equivalent to the validated set-membership semantics. `DnaSequence` (Core, `ValidateSequence`) guarantees the sequence is uppercase A/C/G/T only, so the lookup is sound. Case handling: pattern is upper-cased (`ToUpperInvariant`); sequence is already normalised by `DnaSequence`, and the string-Core variant upper-cases the sequence too.

### Minor robustness note (not a defect for the canonical path)
In `FindDegenerateMotifCore`, the negation-based arms (`'B' => seqChar != 'A'`, `'N' => true`, etc.) would accept a non-ACGT sequence character if the *string* overload were ever called with an unvalidated sequence. This cannot occur on the canonical `(DnaSequence, string)` contract (sequence is validated to ACGT) and the string overload is not contracted to validate sequence chars. No behavioural divergence on any in-contract input; recorded as a note only.

### Cross-verification table recomputed vs code (representative)
| Test | Input | Expected | Code result |
|------|-------|----------|-------------|
| M50 | R / ATGC | [0,2] | ✅ |
| M51 | Y / ATGC | [1,3] | ✅ |
| M52 | N / ACGT | [0,1,2,3] | ✅ |
| M53 | RTG / ATGCGTGC | [0,4] | ✅ |
| M54 | CANNTG / CAGCTG | [0] | ✅ |
| M59 | null / ATG | ArgumentNullException | ✅ |
| M60 | atgc / ATGC | [0] | ✅ |
| S6 | RRG / AAGAG | [0,2] | ✅ |
| S7 | AXG | ArgumentException | ✅ |

### Variant/delegate consistency
Canonical `DnaSequence` overload (dictionary membership) and Core switch (negation form) produce identical results for all ACGT sequences — the two encodings of the same 15-code table agree.

### Test quality audit
`IupacMotifMatchingTests.cs` (74 canonical + others, 117 under the `~Iupac` filter) asserts exact sourced values (positions, matched substrings, exception types), is deterministic, and covers every Stage-A edge case (empty seq/pattern, longer pattern, null, invalid char, case-insensitivity, all 11 codes both positive and negative). Property/snapshot tests strengthen coverage.

### Findings
None. Stage B PASS.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS**, **State: CLEAN.**
- Degeneracy table verified EXACTLY against two authoritative sources; asymmetric (pattern-degenerate) match direction is correct and sourced; all worked examples hand-confirmed and reproduced by the code.
- No code changes. Iupac-filtered tests: 117 passed. Full suite: 4461 passed, 0 failed (baseline preserved).
- Optional future hardening (non-blocking): validate sequence characters in the string-based `FindDegenerateMotifCore` so its negation-based arms cannot accept non-ACGT input. No in-contract impact today.

## Fix applied (2026-06-12)

The optional hardening noted in Stage B was applied. In `MotifFinder.FindDegenerateMotifCore`
(`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs`), the negation-based ambiguity
arms were converted to **positive set membership**, matching the verified IUPAC table exactly:

| Code | Before (negation) | After (positive set) |
|------|-------------------|----------------------|
| B    | `seqChar != 'A'`  | `seqChar is 'C' or 'G' or 'T'` |
| D    | `seqChar != 'C'`  | `seqChar is 'A' or 'G' or 'T'` |
| H    | `seqChar != 'G'`  | `seqChar is 'A' or 'C' or 'T'` |
| V    | `seqChar != 'T'`  | `seqChar is 'A' or 'C' or 'G'` |
| N    | `true`            | `seqChar is 'A' or 'C' or 'G' or 'T'` |

This is now consistent with `IupacHelper.MatchesIupac` (already positive-membership) and with the
canonical `FindDegenerateMotif(DnaSequence, string)` dictionary path. No second divergent table was
introduced; match direction (pattern-degenerate vs sequence base) is unchanged.

**Behavioural impact:** For valid ACGT sequences the result is byte-for-byte identical — the change
only affects genuinely non-ACGT sequence characters, which previously could spuriously match
B/D/H/V/N via the string overload (callable with unvalidated input). Such characters now correctly
fail to match unless an exact literal equal applies.

**Tests:** Added 14 hardening/regression tests (H1–H14) to `IupacMotifMatchingTests.cs` exercising
the string overload: each of B/D/H/V/N matches only its positive set (incl. lowercase); non-ACGT
sequence chars (`X`, sequence `N`, gap `-`) are not matched by B/D/H/V/N; regression cases
(`GAAY`→`GAAC`/`GAAT`, `N`→any ACGT) still pass.

**Verification:** `~Iupac` filter: 135 passed (was 117). Full Genomics.Tests suite: 4484 passed,
0 failed, 0 regressions.
