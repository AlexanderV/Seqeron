# Validation Report: PARSE-FASTQ-001 — FASTQ format parser (incl. Phred quality decoding)

- **Validated:** 2026-06-24   **Area:** FileIO
- **Canonical method(s):** `FastqParser.Parse(content, encoding)`; helpers `DetectEncoding`, `DecodeQualityScores`, `EncodeQualityScores`, `PhredToErrorProbability`, `ErrorProbabilityToPhred`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.IO/FastqParser.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/FastqParserTests.cs` (+ `FastqParser_MutationKillers_Tests.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "FASTQ format"** (https://en.wikipedia.org/wiki/FASTQ_format), fetched fresh 2026-06-24. Confirmed verbatim against the page:
  - 4-line record: Field 1 starts `@` + identifier (+ optional description); Field 2 = raw sequence; Field 3 starts `+` (optionally repeats the identifier); Field 4 = quality, **"must contain the same number of symbols as letters in the sequence."**
  - Phred+33 (Sanger / Illumina 1.8+): offset 33, ASCII 33–126, **Q 0–93**.
  - Phred+64 (Illumina 1.3–1.7): offset 64, ASCII 64–126, **Q 0–62**.
  - Phred formula: **Q = −10·log₁₀(p)**, inverse **p = 10^(−Q/10)**.
- **Cock et al. (2009), NAR 38(6):1767–1771** (DOI 10.1093/nar/gkp1137) — cited in Evidence as the authoritative Sanger FASTQ definition; consistent with the offsets/ranges above.
- **NCBI SRA File Format Guide** — cited for SRA conventions (Phred+33 standard, paired-end /1 /2). Consistent.

### Formula check
- Q→P: `P = 10^(−Q/10)` — matches the Wikipedia inverse equation.
- P→Q: `Q = −10·log₁₀(p)` — matches.
- Offset: Q = ord(char) − 33 (Phred+33 default) / − 64 (Phred+64). Matches.

### Edge-case semantics
- seq length == qual length: an invariant per the spec ("same number of symbols"). Parser enforces it structurally (reads quality up to `sequence.Length` chars).
- Auto-detection heuristic (chars < `@` ⇒ Phred+33; chars > `I` ⇒ Phred+64; ambiguous/empty ⇒ Phred+33): a standard practical heuristic, documented as such in TestSpec/Evidence. Solexa (offset 64, negative Q) is intentionally out of scope.
- Empty / null content ⇒ empty enumerable. Sourced as defined behaviour.

### Independent cross-check (hand + fresh Wikipedia ASCII table)
| Input | Encoding | Expected Q | Wikipedia table | Code |
|-------|----------|-----------|-----------------|------|
| `!` (33) | Phred+33 | 0 | 0 ✓ | 33−33=0 ✓ |
| `I` (73) | Phred+33 | 40 | 40 ✓ | 73−33=40 ✓ |
| `~` (126) | Phred+33 | 93 | 93 ✓ | 126−33=93 ✓ |
| `@` (64) | Phred+33 | 31 | 31 ✓ | 64−33=31 ✓ |
| `@` (64) | Phred+64 | 0 | 0 ✓ | 64−64=0 ✓ |
| `h` (104) | Phred+64 | 40 | 40 ✓ | 104−64=40 ✓ |
| `~` (126) | Phred+64 | 62 | 62 ✓ | 126−64=62 ✓ |

Error probabilities via `10^(−Q/10)`: Q0→1.0, Q10→0.1, Q20→0.01, Q30→0.001, Q40→0.0001 — all verified.

**Stage A findings:** none. Description is correct and authoritative.

## Stage B — Implementation

### Code path reviewed
- `FastqParser.cs:86-129` — 4-line / multi-line parsing.
- `FastqParser.cs:131-139` — `ParseHeader` (first space splits Id/Description).
- `FastqParser.cs:148-165` — `DetectEncoding`.
- `FastqParser.cs:170-184` — `DecodeQualityScores`.
- `FastqParser.cs:189-220` — `EncodeQualityScores`, `PhredToErrorProbability`, `ErrorProbabilityToPhred`.

### Formula realised correctly?
- **Offset:** `int offset = encoding == QualityEncoding.Phred64 ? 64 : 33;` then `scores[i] = Math.Max(0, qualityString[i] - offset)` (lines 175,180). Default 33 — correct. `Math.Max(0,…)` clamps sub-offset chars to Q0 (no negative Q); consistent with Sanger/Illumina scope (Solexa explicitly out of scope).
- **Q→P:** `Math.Pow(10, -phredScore / 10.0)` (line 209) — exact `10^(−Q/10)`.
- **P→Q:** `(int)Math.Round(-10 * Math.Log10(p))`, with `p<=0 ⇒ 93` cap (lines 217-219) — matches `−10·log₁₀(p)` and the documented Q93 max-Sanger cap.
- **4-line / `@` / `+` handling:** header requires line starting `@` (line 96); sequence accumulated until a `+`-prefixed line (line 104); quality accumulated **by length** (`while qualityBuilder.Length < sequence.Length`, line 115), not by sniffing `@`. The well-known FASTQ pitfall (a `@` as first quality char) is correctly avoided.
- **seq==qual length:** structurally guaranteed (quality read up to sequence length). `+` inside sequence and multi-line seq/qual handled (tests `Parse_MultiplePlusLines`, `Parse_MultiLineSequence/Quality`).

### Cross-verification table recomputed vs code (via tests, re-run 2026-06-24)
- `DecodeQualityScores("!I~", Phred33)` → [0,40,93] ✓
- `DecodeQualityScores("@h~", Phred64)` → [0,40,62] ✓
- `EncodeQualityScores` round-trips Phred33 [0..93] and Phred64 [0..62] ✓
- Phred math Q0/Q10/Q20/Q30/Q40 → 1/0.1/0.01/0.001/0.0001 ✓; `ErrorProbabilityToPhred(0)`→93 ✓

### Variant/delegate consistency
- `Parse(string)` / `ParseFile` / `Parse(TextReader)` all delegate to the same core `Parse(TextReader,…)` — consistent.
- Encode/Decode are mutual inverses over the documented ranges (round-trip tests pass).

### Test quality audit
- 132 Fastq-related tests (canonical fixture + mutation-killer fixture) pass; assertions check exact sourced Q values, exact error probabilities, exact sequences/IDs (not just "no throw"). Boundary chars `!`/`I`/`~`/`@`/`h` covered; null/empty/malformed and multi-line edge cases covered.

### Findings / notes (not defects)
1. **`DetectEncoding` heuristic boundary.** Phred+64 is triggered by any char `> 'I'` (ASCII >73). A legitimately high-quality Phred+33 base (e.g. `J`=Q41 .. `~`=Q93) on a read whose first non-`@` char exceeds `I` would be misdetected as Phred+64. This is the inherent limitation of *all* range-based auto-detectors (the `@`–`I` overlap is genuinely ambiguous), and the heuristic is documented as the "standard auto-detection approach" in TestSpec/Evidence. It is a property of the heuristic, not a divergence from the FASTQ spec. Callers wanting determinism pass an explicit `encoding`. → drives the Stage B PASS-WITH-NOTES.
2. **Negative-Q (Solexa) intentionally unsupported.** Sub-offset chars clamp to Q0 via `Math.Max(0,…)`; consistent with the documented Phred+33/Phred+64-only scope.

## Verdict & follow-ups
- **Stage A: PASS**, **Stage B: PASS-WITH-NOTES**, **State: CLEAN.** No defects found; no code changes made.
- Tests: `--filter FullyQualifiedName~Fastq` → 132 passed / 0 failed. Build: succeeded, 0 warnings.
- The two notes are documented, sourced, deliberate design choices (heuristic detection + Phred+33/64-only scope), not implementation defects.
