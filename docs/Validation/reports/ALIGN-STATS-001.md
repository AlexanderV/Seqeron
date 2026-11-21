# Validation Report: ALIGN-STATS-001 — Pairwise Alignment Statistics (Identity / Similarity / Gaps) and Formatting

- **Validated:** 2026-06-15; re-validated 2026-11-18 (B09)   **Area:** Alignment
- **Canonical method(s):** `SequenceAligner.CalculateStatistics(AlignmentResult, ScoringMatrix? | SubstitutionMatrix)`, `SequenceAligner.FormatAlignment(AlignmentResult, int, ScoringMatrix? | SubstitutionMatrix)`, `SubstitutionMatrix.Blosum62 / Pam250 / EdnaFull`, MCP `alignment_statistics`, `format_alignment`
- **Stage A verdict (B09 re-validation):** PASS-WITH-NOTES (2026-06 pass accepted a wrong markup claim — see below)
- **Stage B verdict (B09 re-validation):** PASS-WITH-NOTES → FIXED (F11–F14)

## B09 current state (2026-11-19) - supersedes the original verdict lines above and below

- **Stage A:** PASS-WITH-NOTES (the 2026-06 markup claim was refuted by needle / water 6.6.0 and `ajalign.c`). **Stage B:** PASS-WITH-NOTES. **State:** FIXED (not CLEAN; the 2026-06 'CLEAN' below is superseded).
- **Fixes (B09.md, F-numbers):** F11 case-insensitive residues, F12 unequal rows rejected, F13 gap-opening count + positives, F14 residue-pair similarity and the exact EMBOSS markup, F42 BLAST integer percentages / positives, F43 statistics of > 2 sequences (`MultipleAlignmentStatistics`).
- **Oracle numbers:** needle / water 301 alignments (12 439 columns): Length / Identity / Similarity / Gaps, mismatches, open gaps and the markup line equal; BLAST+ 2.12.0 320 / 320; infoalign 3 268 / 3 268; distmat 600 / 600 (+ 1 100 / 1 100 for every correction and option, A2-15); ClustalW 151 / 151; Biopython identity 77 / 77; matrices 576 / 576 / 256 cells equal the EMBOSS data files; 70 / 70 matrices of the EMBOSS data directory (A2-15).
- **Open items:** none. NOTE-1 of the original pass ('`:` for any positive score, `.` unreachable') and 'M7 follows the simplification' are REFUTED and no longer true: the code and M7 follow the EMBOSS rule exactly.
- The unit's TestSpec (`tests/TestSpecs/ALIGN-STATS-001.md`), Evidence (`docs/Evidence/ALIGN-STATS-001-Evidence.md`) and algorithm document carry the same state; where the original text below says CLEAN / no defects / no code changes, that was the state of the first pass only.

## B09 re-validation (2026-11-18) — summary

Stage A. Sources actually opened: EMBOSS 6.6 `ajax/core/ajalign.c` (raw.githubusercontent.com; `alignConsStats`, `alignSim`), the local needle/water 6.6.0 binaries and EMBOSS data files (`EBLOSUM62`, `EPAM250`, `EDNAFULL`), Biopython 1.88 (`Alignment.counts`, `substitution_matrices`). The 2026-06 pass relied on the AlignFormats web page (`:` > 1.0, `.` small positive, space mismatch) and on the claim that `.` is unreachable; the binary and the source refute both: needle writes `.` for every aligned non-similar column and space only for gap columns, for identical pairs whose own score is not positive, and (an EMBOSS quirk) for first-row ambiguity residues. The Similarity header counts `M(x,y) > 0 ∧ (M(x,x) > 0 ∨ M(y,y) > 0)`. Denominator = alignment length including gap columns (confirmed on 301 alignments).

Stage B. Findings (details in B09.md): **F11** case-sensitive residue comparison; **F12** unequal rows accepted; **F13** no gap-opening / positives output (BLAST `gapopen`, Biopython `open_gaps` parity); **F14** no residue-pair similarity (protein / IUPAC) — added `SubstitutionMatrix` (BLOSUM62, PAM250, EDNAFULL, every cell equal to the EMBOSS data files and Biopython), exact EMBOSS similarity/markup rule, MCP `matrix` / `match` / `mismatch` parameters, `positives`, `gapOpens`; MCP `format_alignment` description claimed `.` for mismatches while the code printed a space.

Oracle: `TestData/Alignment/needle_water_stats.tsv` (301 alignments, 12 439 columns; generator `scripts/alignment/generate_alignment_stats_oracle.py`, deterministic): Length, Identity, Similarity, Gaps (EMBOSS header), mismatches and open gaps (Biopython), and the markup line of every alignment equal. HBA/HBB published example: Length 149, 65 / 90 / 9 reproduced with the real BLOSUM62 (4 gap runs). Matrices: 576 / 576 / 256 cells equal the EMBOSS data files; 576 / 576 / 225 shared cells equal Biopython 1.88.

Scalar presets verified: `BlastDna` 2 / −3 / −5 / −2 = blastn defaults (reward 2, penalty −3, existence 5, extension 2 under `GapOpen + k·GapExtend`); `HighIdentityDna` match/mismatch = EDNAFULL on ACGT.

---

## Original 2026-06-15 pass (superseded where it conflicts with the above)


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
- **NOTE-1 (REFUTED 2026-11-18 - see the B09 re-validation above; kept as history):** EMBOSS marks `:` for substitution score **> 1.0** and `.` for small positive scores (0 < score ≤ 1.0). The description (algorithm doc §4.2/§5.3, TestSpec Assumption #1) declares the `.` tier "unreachable" for the scalar integer `Match`/`Mismatch` model and renders any positive non-identical column as `:`. This is a rendering-only simplification, explicitly disclosed, with no effect on the counted statistics. Honest scoping → PASS-WITH-NOTES, not FAIL. The description's claim "positive substitution score ⇒ similar" for the **counted** statistics is fully correct (EMBOSS/BLAST positives definition is "scores ≥/> 0" for the count; the >1.0 threshold is a *display* nuance specific to the srspair markup glyph).

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
- **Result: PASS.** No green-washing, no weakened assertion, no skip. Note (2026-06, REFUTED 2026-11-18): M7 was written to the `:`-for-any-positive simplification; the 2026-11 re-validation replaced that rule by the exact EMBOSS markup (needle / water 6.6.0 + `ajalign.c`) and M7 now follows it.

### Findings / defects (Stage B)
- No algorithm defect. No code or test change required.
- **NOTE-1 (carry-through, REFUTED 2026-11-18; `FormatAlignment` now implements the exact EMBOSS rule, F14):** the srspair `:`/`.` threshold simplification (Stage A NOTE-1) is realised in `FormatAlignment` (`score.Mismatch > 0 ⇒ ':'`). Rendering-only; counted statistics unaffected. Documented in algorithm doc §5.3 and TestSpec Assumption #1.

## Verdict & follow-ups
- **Stage A (2026-06):** PASS-WITH-NOTES (markup `:`/`.` threshold simplification - refuted 2026-11-18).
- **Stage B (2026-06):** PASS-WITH-NOTES (same display nuance; statistics exact and sourced) - the markup was fixed in F14.
- **End-state (2026-06):** CLEAN at the time; superseded by the B09 re-validation above (F11–F14 fixed; fast tier see B09.md).
- **Test-quality gate:** PASS.
- No logged defects (BY-DESIGN note recorded in FINDINGS_REGISTER §D).

## 2026-11 audit round 1 (A1-9, B09 F42/F43)

BLAST+ integer percentages / `BlastPositives` / `pident` text / report line on `AlignmentStatistics` and `MultipleAlignmentStatistics` (infoalign, distmat, ClustalW conservation, pairwise identity) with MCP `multiple_alignment_statistics`, `alignment_distance_matrix` and the BLAST fields of `alignment_statistics`; oracle BLAST+ 2.12.0 (320 HSPs), infoalign/distmat 6.6.0, ClustalW 2.1, Biopython 1.88. Details: B09.md F42/F43, algorithm doc §5.5-5.6.

- **Audit r3 (A3-13c, F148):** EMBOSS `showalign` ported (`MultipleAlignmentDisplay.ShowAlign`, MCP `showalign`): 600 / 600 texts equal; the consensus core of `cons` shared (gap option, case kept).

- **Audit A4-3 (F153):** the EMBOSS reader drops gap characters and every other non-letter (except `*`) from plain input sequences; one shared `EmbossSequenceReader` (Infrastructure) serves showalign: 240 oracle rows each equal (`emboss_gap_oracle.tsv`).

- **Audit A3-13c2 (F191):** `showalign -html` / `-highlight` ported (`EmbossRangeSpec`, `ShowAlignOptions.Html` / `Highlight`): 300 / 300 HTML texts equal.
- **Audit A4-1 (F157):** `MultipleAlignmentDisplay.ShowAlign` poll the ambient cancellation token (mid-run cancel tests); MCP tools take the request token.

- **Audit r5 / r6 (A5-3, A5-7, A6-2, A6-4; F169, F172, F177, F179):** the null / out-of-range contract of the `MultipleAlignmentStatistics` entry points (`SequenceAligner_NullContract_Tests`, MCP `NullContractTests`); `distmat` reads the sequence qualifiers (`-sbegin` / `-send` have no effect on it, `-sreverse` reverse-complements every row, the set is nucleotide when its first row is: `emboss_seqopts_oracle.tsv`, 180 rows), the qualifier TEXT of the `distmat` MCP tool is resolved as the program resolves it (`EmbossSequenceQualifiers.ParseInputs("distmat", text)`: 130 + 40 oracle rows).
