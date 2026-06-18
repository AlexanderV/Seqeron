# Validation Report: ALIGN-STATS-001 — Pairwise Alignment Statistics (Identity / Similarity / Gaps) and Formatting

- **Validated:** 2026-06-15   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.CalculateStatistics(AlignmentResult, ScoringMatrix?)`, `SequenceAligner.FormatAlignment(AlignmentResult, int, ScoringMatrix?)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened this session
- **EMBOSS needle docs, rel 6.6** (WebFetch, 2026-06-15) — `https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html`. Verbatim worked-example block retrieved:
  ```
  # Length: 149
  # Identity:      65/149 (43.6%)
  # Similarity:    90/149 (60.4%)
  # Gaps:           9/149 ( 6.0%)
  # Score: 292.5
  ```
  Definitions retrieved verbatim: Identity = "the percentage of identical matches between the two sequences over the reported aligned region (including any gaps in the length)"; Similarity = "the percentage of matches … (including any gaps in the length)". The denominator (149) includes gap columns for all three metrics.
- **EMBOSS Alignment Formats (srspair/pair markup legend)** (WebFetch + WebSearch, 2026-06-15) — `https://emboss.sourceforge.net/docs/themes/AlignFormats.html`. Retrieved legend: `|` = identity (same residue regardless of score); `:` = similarity **scoring more than 1.0**; `.` = **any small positive score**; space = mismatch or gap.
- **pseqsid reference implementation** (WebFetch, 2026-06-15) — `https://github.com/amaurypm/pseqsid`. Confirms denominator options incl. `alignment` (= "includes gap-only columns"); similarity = identical + similar residues; similarity grouping is via user-defined groups / matrices.

### Formula check
The algorithm doc and TestSpec state, matching EMBOSS source verbatim:
- Identity% = M / L × 100 (L includes gap columns) ✓
- Similarity% = (M + Sim⁺) / L × 100, where a non-identical column is similar iff its substitution score is **positive** ✓ (EMBOSS/BLAST "positives" rule, corroborated by pseqsid)
- Gaps% = G / L × 100 ✓
- Partition M + X + G = L ✓

### Independent cross-check (numbers re-derived this session, Python)
| Quantity | Source counts | Computed % | Published % |
|----------|---------------|-----------|-------------|
| Identity | 65/149 | 43.624161… | 43.6 ✓ |
| Similarity | 90/149 | 60.402685… | 60.4 ✓ |
| Gaps | 9/149 | 6.040268… | 6.0 ✓ |
| Hand M5 Identity | 6/9 | 66.6667 | (formula) ✓ |
| Hand M5 Gap | 2/9 | 22.2222 | (formula) ✓ |

Hand-alignment `ACGT-ACGT` / `ACCTAAC-T` re-classified in Python: 6 match, 1 mismatch, 2 gap, L=9 — matches Evidence dataset.

### Edge-case semantics
Empty → `AlignmentStatistics.Empty` / `""`; null → `ArgumentNullException`; `lineWidth ≤ 0` → `ArgumentOutOfRangeException`; all-gap → Gaps==L, Identity 0%; perfect identity → 100/100/0. All defined and consistent with EMBOSS conventions (denominator undefined for L=0).

### Findings / divergences (Stage A)
- **NOTE-1 (markup threshold, documented):** EMBOSS marks `:` for substitution score **> 1.0** and `.` for small positive scores (0 < score ≤ 1.0). The description (algorithm doc §4.2/§5.3, TestSpec Assumption #1) declares the `.` tier "unreachable" for the scalar integer `Match`/`Mismatch` model and renders any positive non-identical column as `:`. This is a rendering-only simplification, explicitly disclosed, with no effect on the counted statistics. Honest scoping → PASS-WITH-NOTES, not FAIL. The description's claim "positive substitution score ⇒ similar" for the **counted** statistics is fully correct (EMBOSS/BLAST positives definition is "scores ≥/> 0" for the count; the >1.0 threshold is a *display* nuance specific to the srspair markup glyph).

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs:566-618` (`CalculateStatistics`) and `:630-673` (`FormatAlignment`).
- `src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/AlignmentTypes.cs:8-53` (`ScoringMatrix`, `AlignmentResult.Empty`, `AlignmentStatistics.Empty`).

### Formula realised correctly?
Yes. Single O(L) pass classifies each column: gap if either char is `-`; else identical if equal; else mismatch. A mismatch increments the similar counter iff `score.Mismatch > 0`. The three percentages use denominator `alignmentLength` (includes gaps). For the scalar DNA model every non-identical pair scores exactly `Mismatch`, so `score.Mismatch > 0` is **mathematically exact** (not an approximation) for "is this mismatch column's substitution score positive". The scalar `ScoringMatrix` has only `Match`/`Mismatch` (no per-pair matrix), so the same test applies uniformly. SimpleDna/BlastDna/HighIdentityDna all have `Mismatch < 0` ⇒ Similarity == Identity, as documented.

### Cross-verification table recomputed vs code (via the test suite)
| Case | Expected (sourced) | Result |
|------|--------------------|--------|
| EMBOSS 149-col (PositiveMismatch, 65 id + 25 sim-mm + 59 gap) | Identity 43.6%, Similarity 60.4%, L=149, Matches=65 | PASS |
| SimpleDna ACGT/ACCT | Identity 75%, Similarity==Identity | PASS |
| PositiveMismatch ACGT/ACCT | Identity 75%, Similarity 100% | PASS |
| Hand 9-col | M6 X1 G2 L9, Id 66.67%, Sim 66.67%, Gap 22.22% | PASS |
| All-gap | Gaps==L, Id 0%, Gap 100% | PASS |
| FormatAlignment SimpleDna ACGT/ACCT | markup `\|\| \|` | PASS |
| FormatAlignment PositiveMismatch | markup `\|\|:\|` | PASS |

### Variant/delegate consistency
`CalculateStatistics` with default scoring (M2b) matches explicit SimpleDna; `FormatAlignment` default `lineWidth`/`scoring` exercised. Both public methods and both optional-parameter paths covered.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** M1 reproduces the published 43.6/60.4 (counts 65/149, 90/149) — a deliberately-wrong denominator or similarity rule fails it. M3 locks 75%→100% (positive-score similarity). M5 locks exact 6/9 and 2/9. FormatAlignment markup strings are the literal srspair glyphs. The `Similarity Is.EqualTo(Identity)` assertions in M2/M2b/M5 are **anchored** by independent exact-value Identity asserts in the same test, so they are not tautologies.
- **No green-washing:** exact `.Within(1e-10)` / `1e-3` tolerances throughout; the only ranged asserts (M1 `.Within(0.05)` "rounds to 43.6") sit alongside the exact `65.0/149.0*100 .Within(1e-3)` assert and mirror the published 1-dp rounding — not a weakening. No skip/ignore.
- **Coverage:** all 7 MUST (M1–M7), all 6 SHOULD (S1–S6), both COULD (C1–C2) present; both public methods, default + explicit scoring, null/empty/boundary (lineWidth 0 and −1), all-gap, perfect-identity, line-wrap blocks, and INV-1/2/3 all exercised.
- **Result: PASS.** No green-washing, no weakened assertion, no skip. Note: M7 follows the documented `:`-for-any-positive simplification rather than strict EMBOSS `.`-for-score-1.0 — this is the disclosed rendering assumption, not a test defect.

### Findings / defects (Stage B)
- No algorithm defect. No code or test change required.
- **NOTE-1 (carry-through):** the srspair `:`/`.` threshold simplification (Stage A NOTE-1) is realised in `FormatAlignment` (`score.Mismatch > 0 ⇒ ':'`). Rendering-only; counted statistics unaffected. Documented in algorithm doc §5.3 and TestSpec Assumption #1.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (markup `:`/`.` threshold simplification, documented).
- **Stage B:** PASS-WITH-NOTES (same display nuance; statistics exact and sourced).
- **End-state:** ✅ CLEAN — no defect; full unfiltered suite **6536 passed / 0 failed**; build 0 errors.
- **Test-quality gate:** PASS.
- No logged defects (BY-DESIGN note recorded in FINDINGS_REGISTER §D).
