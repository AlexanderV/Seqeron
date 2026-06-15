# Validation Report: QUALITY-PHRED-001 — Phred Score Handling

- **Validated:** 2026-06-15   **Area:** Quality
- **Canonical method(s):** `QualityScoreAnalyzer.ParseQualityString(qualStr, encoding)`,
  `QualityScoreAnalyzer.ToQualityString(scores, encoding)`,
  `QualityScoreAnalyzer.ConvertEncoding(qualStr, from, to)`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.IO/QualityScoreAnalyzer.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened (this session) & what they confirm

1. **Cock, P.J.A. et al. (2010), *Nucleic Acids Research* 38(6):1767–1771** — the de-facto FASTQ
   format specification. Retrieved this session via
   `https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/`. Verbatim extracts:
   - Sanger / Phred+33: *"Sanger FASTQ files use ASCII 33–126 to encode PHRED qualities from 0 to 93
     (i.e. PHRED scores with an ASCII offset of 33)."* → offset 33, Q ∈ [0, 93], ASCII ∈ [33, 126].
   - Illumina 1.3+ / Phred+64: *"The Illumina 1.3+ FASTQ variant encodes PHRED scores with an ASCII
     offset of 64, and so can hold PHRED scores from 0 to 62 (ASCII 64–126)."* → offset 64,
     Q ∈ [0, 62], ASCII ∈ [64, 126].
   - Phred formula: *"Q = −10 log₁₀(P)"* where P is the probability of error.
   - Solexa (out of scope): *"Q_solexa = −10 log₁₀(P/(1−P))"*, offset 64, ASCII 59–126, scores −5 to 62.
   - Invariance: *"Conversion from 'fastq-illumina' to 'fastq-sanger' … is very straightforward since
     both variants use PHRED scores but with different offsets."*

2. **Biopython `Bio.SeqIO.QualityIO` docs** (`https://biopython.org/docs/1.75/api/Bio.SeqIO.QualityIO.html`)
   — independent reference implementation. Confirms: `fastq`/`fastq-sanger` = PHRED, offset 33;
   `fastq-illumina` = PHRED, offset 64; `fastq-solexa` = Solexa scores, offset 64. Phred score is
   preserved between sanger and illumina variants.

3. **Wikipedia, *FASTQ format*** (`https://en.wikipedia.org/wiki/FASTQ_format`) — independent
   cross-check of individual characters: under Phred+33, ASCII 95 (`_`) = Q62 and ASCII 126 (`~`) = Q93;
   under Phred+64, ASCII 126 (`~`) = Q62, ASCII 64 (`@`) = Q0. Offsets 33 and 64 confirmed.

### Formula check

- Decode `Q = ord(char) − offset`, encode `char = chr(Q + offset)` — matches Cock et al. exactly.
- Phred quality definition `Q = −10·log₁₀(P)` — matches.
- Cross-variant conversion = pure byte re-offset (shift ±31), Phred score invariant — matches.

### Edge-case semantics check

- Char below offset → negative Q → malformed for the variant (Phred Q ≥ 0). Sourced (Cock et al.).
- Phred+64 (Q 0–62) → Phred+33 (Q 0–93): always representable. Sourced.
- Phred+33 Q ∈ (62, 93] → Phred+64: not representable (target max is 62). Sourced.
- The .NET exception *types* (`ArgumentNullException`, `ArgumentOutOfRangeException`) are an API-shape
  choice declared as an assumption in the spec/Evidence; the *range bounds* themselves are source-backed.

### Independent cross-check (hand computation + Wikipedia + Biopython)

| Op | Input | Source-derived expected | Confirmed by |
|----|-------|-------------------------|--------------|
| Parse Phred+33 | `!~` | [0, 93] | Cock (33→0, 126→93); Wikipedia |
| Parse Phred+33 | `5?I` | [20, 30, 40] | ASCII 53/63/73 − 33 |
| Parse Phred+64 | `@h~` | [0, 40, 62] | Cock (64→0, 104→40, 126→62); Wikipedia |
| Encode Phred+33 | [0,20,30,40,93] | `!5?I~` | chr(Q+33) |
| Encode Phred+64 | [0,40,62] | `@h~` | chr(Q+64) |
| Convert P64→P33 | `@h~` (Q0,40,62) | `!I_` | chr(Q+33)=33,73,95; Wikipedia `_`=Q62@P33 |
| Convert P33→P64 | `!I` (Q0,40) | `@h` | chr(Q+64)=64,104 |
| Parse below offset | `" "`=32 @P33 | Q=−1 → reject | Cock corner case 1 |
| Convert P33→P64 | `~`=Q93 | overflow (max 62) → reject | Cock corner case 3 |
| Encode | [94] @P33 | >93 → reject | Cock max Q 93 |

Every non-trivial expected value traces to Cock et al. and/or Wikipedia/Biopython retrieved this session,
not to the implementation's output.

### Findings / divergences

None. The description (Phred_Score_Handling.md), TestSpec, and Evidence all match the primary source
and two independent references.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.IO/QualityScoreAnalyzer.cs`:
- `ParseQualityString` (L100–125): null → `ArgumentNullException`; empty → `Array.Empty<int>()`;
  per char `score = char − offset`; validates `score ∈ [0, maxScore]` else `ArgumentOutOfRangeException`.
- `ToQualityString` (L137–159): null → `ArgumentNullException`; empty → `string.Empty`;
  validates `score ∈ [0, maxScore]` else `ArgumentOutOfRangeException`; emits `(char)(score + offset)`.
- `ConvertEncoding` (L173–180): null check; `ParseQualityString(from)` then `ToQualityString(to)` —
  overflow on the target is surfaced naturally by `ToQualityString`'s range check.
- `GetEncodingParameters` (L81–86): Phred64 → (64, 62), else → (33, 93). Constants (L67–71) cite Cock et al.

### Formula realised correctly?

Yes. The code computes the exact validated formula (`ord−offset` / `chr(Q+offset)`), with the exact
sourced ranges, for all three canonical methods. No approximation.

### Cross-verification table recomputed vs code (via the test run)

All 10 MUST + 6 SHOULD + 1 COULD rows above pass against the actual code (full suite green, below).
M6 `@h~`→`!I_`, M7 `!I`→`@h`, the negative-Q reject, the Q93→Phred+64 overflow reject, and the
encode-94 reject all behave exactly as the source prescribes.

### Variant/delegate consistency

The legacy helpers `CharToPhred`/`PhredToChar`/`QualityStringToPhred`/`PhredToQualityString` use the same
offsets but perform **no** range validation (by design — used by statistics/trim helpers, out of scope
for this unit). The canonical trio adds the source-backed validation. No inconsistency in the offset math.

### Test quality audit (against the HARD gate)

- **Sourced expectations, not code echoes** — every assertion uses `Is.EqualTo` with the
  externally-confirmed exact values/strings (e.g. `[0,93]`, `"!I_"`, `"@h"`). A deliberately off-by-one
  offset would fail these. PASS.
- **No green-washing** — no `Greater`/`AtLeast`/`Contains`/range assertions; exception tests assert the
  exact type via `Assert.Throws<T>`; nothing skipped/ignored/widened. PASS.
- **Cover all logic** — all 3 public canonical methods exercised; both Phred+33 and Phred+64 paths;
  boundary + interior decode; encode; both conversion directions + identity; all three documented
  error cases (below-offset decode M8, encode-overflow M10, convert-overflow M9); null (S3/S4/S5) and
  empty (S1/S2) boundaries; round-trip property (C1). 17 tests = full spec. The only un-exercised branch
  is the `Auto` encoding path inside `ParseQualityString` — explicitly scoped out by spec Open Question #1
  (callers pass an explicit encoding; `Auto`/`DetectEncoding` is a separate out-of-scope method). Acceptable.
- **Honest green** — FULL unfiltered suite: **Failed: 0, Passed: 6510, Skipped: 0**; `dotnet build`
  0 warnings / 0 errors. PASS.

**Test-quality gate result: PASS.**

### Findings / defects

None. No code or test changes were required this session.

## Verdict & follow-ups

- **Stage A: PASS** — description/formula/ranges/edge-cases independently confirmed against Cock et al.
  (2010), Biopython, and Wikipedia.
- **Stage B: PASS** — implementation faithfully realises the validated formula; tests assert exact
  sourced values and cover every public method, both encodings, both conversion directions, and all
  documented error/boundary cases.
- **End-state: ✅ CLEAN** — no defect found; algorithm fully functional.
- No follow-ups. (Note: the `Auto`/`DetectEncoding` heuristic is a separate, out-of-scope method for
  this unit and is not validated here.)
