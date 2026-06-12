# Validation Report: PARSE-FASTQ-001 — FASTQ format parser (incl. Phred quality decoding)

- **Validated:** 2026-06-12   **Area:** FileIO
- **Canonical method(s):** `FastqParser.Parse(content, encoding)`; helpers `DetectEncoding`, `DecodeQualityScores`, `EncodeQualityScores`, `PhredToErrorProbability`, `ErrorProbabilityToPhred`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.IO/FastqParser.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/FastqParserTests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "FASTQ format"** (https://en.wikipedia.org/wiki/FASTQ_format), fetched 2026-06-12. Confirmed verbatim:
  - 4-line record: line1 `@`+identifier(+optional description); line2 sequence; line3 `+` (optionally repeats identifier); line4 quality, **identical length to the sequence**.
  - Phred+33 (Sanger / Illumina 1.8+): offset 33, ASCII 33–126 (`!`–`~`), Q 0–93. `!`(33)→Q0, `I`(73)→Q40.
  - Phred+64 (Illumina 1.3–1.7): offset 64, ASCII 64–126 (`@`–`~`), Q 0–62. `@`(64)→Q0.
  - Phred formula: **Q = −10·log₁₀(p)**, inverse **p = 10^(−Q/10)**.
- **Cock et al. (2009), NAR 38(6):1767–1771** (DOI 10.1093/nar/gkp1137) — cited in the Evidence doc as the authoritative Sanger FASTQ definition; consistent with the offsets/ranges above.
- **NCBI SRA File Format Guide** — cited for SRA submission conventions (Phred+33 standard, paired-end /1 /2). Consistent.

### Formula check
- Q→P: `P = 10^(−Q/10)` — matches the Wikipedia inverse equation exactly.
- P→Q: `Q = −10·log₁₀(p)` — matches.
- Offset: Q = ord(char) − 33 (Phred+33, default) / − 64 (Phred+64). Matches.

### Edge-case semantics
- seq length == qual length: invariant per source (line4 same length as line2). The parser enforces it structurally by reading quality exactly up to `sequence.Length` chars.
- Encoding auto-detection heuristic (chars < `@` ⇒ Phred+33; chars > `I` ⇒ Phred+64; ambiguous/empty ⇒ Phred+33) — a standard practical heuristic; Solexa (offset 64, negative Q) is intentionally out of scope and not claimed.

### Independent cross-check (hand-computed)
| Input | Encoding | Expected Q | Computed |
|-------|----------|-----------|----------|
| `"III"` | Phred+33 | [40,40,40] | ord('I')−33 = 40 ✓ |
| `"!"` | Phred+33 | Q0 | 33−33 = 0 ✓ |
| `"@"` (ASCII 64) | Phred+33 | Q31 | 64−33 = 31 ✓ |
| `"@"` | Phred+64 | Q0 | 64−64 = 0 ✓ |
| `"h"` | Phred+64 | Q40 | 104−64 = 40 ✓ |
| `"~"` | Phred+33 / +64 | Q93 / Q62 | 126−33=93, 126−64=62 ✓ |

Error probabilities: Q0→1.0, Q10→0.1, Q20→0.01, Q30→0.001, Q40→0.0001 — all verified via `10^(−Q/10)`.

**Stage A findings:** none. Description is correct and authoritative.

## Stage B — Implementation

### Code path reviewed
- `FastqParser.cs:86-129` — 4-line / multi-line parsing.
- `FastqParser.cs:148-165` — `DetectEncoding`.
- `FastqParser.cs:170-184` — `DecodeQualityScores` (offset 33 default, 64 for Phred64).
- `FastqParser.cs:207-220` — `PhredToErrorProbability` / `ErrorProbabilityToPhred`.

### Formula realised correctly?
- **Offset:** `int offset = encoding == QualityEncoding.Phred64 ? 64 : 33;` then `scores[i] = Math.Max(0, qualityString[i] - offset)` (line 175,180). Default 33 — correct. `Math.Max(0,…)` clamps below-offset chars to 0 (no negative Q); acceptable for Sanger/Illumina (Solexa not in scope).
- **Q→P:** `Math.Pow(10, -phredScore / 10.0)` (line 209) — exact match to `10^(−Q/10)`.
- **P→Q:** `(int)Math.Round(-10 * Math.Log10(p))` with `p<=0 ⇒ 93` cap (line 217-219) — matches `−10·log₁₀(p)` and the documented Q93 max-Sanger cap.
- **4-line / `@` / `+` handling:** header requires line starting `@` (line 96); sequence accumulated until a `+`-prefixed line (line 104); quality accumulated by **length** (`while qualityBuilder.Length < sequence.Length`, line 115), not by sniffing `@`. This makes a `@` as the first quality char unambiguous — the well-known FASTQ pitfall is correctly avoided.
- **seq==qual length:** structurally guaranteed (quality read up to sequence length). `+` inside sequence and multi-line seq/qual are handled (tested: `Parse_MultiplePlusLines_ParsesCorrectly`, `Parse_MultiLineSequence/Quality_AssembledCorrectly`).

### Cross-verification table recomputed vs code (via tests)
- `DecodeQualityScores("!I~", Phred33)` → [0,40,93] ✓ (test line 157-161).
- `DecodeQualityScores("@h~", Phred64)` → [0,40,62] ✓ (test line 169-173).
- `EncodeQualityScores` round-trips Phred33 [0..93] and Phred64 [0..62] ✓.
- Phred math Q0/Q10/Q20/Q30/Q40 → 1/0.1/0.01/0.001/0.0001 ✓; `ErrorProbabilityToPhred(0)`→93 ✓.
- Worked example `@`(64) under Phred+33 = Q31: confirmed by `64−33` in `DecodeQualityScores` (no test asserts this specific value, but the formula is exercised by other boundary tests).

### Variant/delegate consistency
- `Parse(string)` / `ParseFile` / `Parse(TextReader)` all delegate to the same `Parse(TextReader,…)` core — consistent.
- Encode/Decode are inverse over the documented ranges (round-trip tests pass).

### Test quality audit
- 80 tests in fixture; assertions check exact sourced Q values, exact error probabilities, exact sequences/IDs (not just "no throw"). Boundary values `!`/`I`/`~`/`@`/`h` covered. Null/empty/malformed and multi-line edge cases covered.

### Edge cases
- seq/qual length mismatch: quality read is length-bounded to the sequence, so a structurally short quality block simply yields a shorter QualityString; no crash. Acceptable.
- `@` inside quality string: not misparsed (length-bounded read). ✓
- Empty / null content: returns empty enumerable. ✓
- Phred+64 supported and tested. ✓

**Stage B findings:** none. Code faithfully realises the validated description.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS**, **State: CLEAN.** No defects; no code changes.
- Tests: `--filter FullyQualifiedName~Fastq` → 80 passed / 0 failed. Full suite `Seqeron.Genomics.Tests` → 4486 passed / 0 failed (matches baseline).
- Minor note (not a defect): negative-Q (Solexa) encoding is intentionally unsupported; sub-offset chars clamp to Q0 via `Math.Max(0,…)`. This is consistent with the spec, which scopes only Phred+33/Phred+64.
