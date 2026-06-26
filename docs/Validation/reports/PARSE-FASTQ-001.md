# Validation Report: PARSE-FASTQ-001 — FASTQ format parser (incl. Phred quality decoding)

- **Validated:** 2026-06-24; **re-validated:** 2026-06-26 (file-level encoding detector, fix `7977bdde`)   **Area:** FileIO
- **Canonical method(s):** `FastqParser.Parse(content, encoding)`; helpers `DetectEncoding`, `DecodeQualityScores`, `EncodeQualityScores`, `PhredToErrorProbability`, `ErrorProbabilityToPhred`. **Re-validation focus:** `QualityScoreAnalyzer.DetectEncoding(IEnumerable<string>)` (new file-level / multi-read Phred-offset detector + `EncodingConfidence`).
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.IO/FastqParser.cs`; `src/Seqeron/Algorithms/Seqeron.Genomics.IO/QualityScoreAnalyzer.cs`
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/FastqParserTests.cs` (+ `FastqParser_MutationKillers_Tests.cs`); `QualityScoreAnalyzer_DetectEncodingFile_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

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

---

## Re-validation 2026-06-26 — file-level Phred-offset detector (`QualityScoreAnalyzer.DetectEncoding(IEnumerable<string>)`, fix `7977bdde`)

A new file-level overload was added that aggregates the global min/max ASCII over all reads and applies the Cock et al. (2010) ranges, returning an `EncodingDetectionResult(Encoding, Confidence, MinAscii, MaxAscii, CharactersExamined)` with `EncodingConfidence { Definitive, Inferred, Ambiguous }`. This was independently re-validated this session (sources retrieved fresh — NOT lifted from the repo Evidence).

### Stage A — Description vs external sources (verdict: PASS)

Sources retrieved this session:
- **Cock et al. (2010), "The Sanger FASTQ file format…", NAR 38(6):1767–1771**, DOI 10.1093/nar/gkp1137 (fetched via academic.oup.com). Confirms the canonical table: **Sanger/Phred+33 = ASCII 33–126, offset 33, Q 0–93; Illumina 1.3+/Phred+64 = ASCII 64–126, offset 64, Q 0–62; Solexa = ASCII 59–126, offset 64, Q −5–62.** The paper states the variants "cannot be reliably distinguished" from content alone, and Solexa/Illumina-1.3+ share ASCII 64–126 — the explicit basis for "Phred+64 can only be inferred."
- **Wikipedia "FASTQ format"** (fetched) — same ASCII/offset/Q table for all five variants (Sanger, Solexa, Illumina 1.3+/1.5+/1.8+); Illumina 1.8+ reverts to offset 33 with typical raw range Q 0–41.
- **qcfail.com, "Incorrect encoding of Phred scores"** (fetched) — confirms the information-theoretic residual verbatim: "because the phred+33 and phred+64 ranges overlap it is possible to generate datasets which are compatible with both encoding schemes" (e.g. aggressively trimmed sets). This is exactly the LIMITATIONS §1 irreducible case.
- **FastQC heuristic** (web search) — real tools key off the lowest observed quality character: scores in the low (33–63) range ⇒ Sanger/Phred+33, otherwise Illumina 1.3–1.7; presence of 'J' ⇒ Illumina 1.8+ Phred+33. Matches the implementation's "lowest char decides" logic.
- **Biopython 1.85** (`Bio.SeqIO.QualityIO`): `SANGER_SCORE_OFFSET = 33`, `SOLEXA_SCORE_OFFSET = 64`. Independent confirmation of the offsets.

Boundary arithmetic independently recomputed (Biopython / hand):
- ASCII 63 = '?' (below 64 ⇒ outside Phred+64 ⇒ proves Phred+33); ASCII 64 = '@' (Phred+64 Q0, overlap start).
- **Illumina 1.8+ ceiling Q41 → chr(41+33) = 'J' = ASCII 74.** (The earlier Wikipedia auto-extract erroneously reported ')'/ASCII 41; corrected here: 33+41 = 74 = 'J'.) ASCII 75 = 'K' is therefore above the Phred+33 ceiling ⇒ infers Phred+64.

Judgement on the sharpened contract:
- **"A char < 64 PROVES Phred+33; Phred+64 can only be INFERRED" — information-theoretically correct.** Phred+64's valid ASCII range is 64–126, so any char < 64 is impossible under Phred+64 ⇒ Phred+33 is the only consistent encoding (a proof). Conversely every Phred+64 string (chars 64–126) is also a syntactically valid Phred+33 string, so Phred+64 can never be proven by range alone — only inferred from an implausibly high char. Asymmetry is correct.
- **The [64,74] overlap boundary is correct.** 64 = Phred+64 floor; 74 = 'J' = Phred+33's Illumina-1.8 ceiling (Q41). A char strictly above 74 is implausible under Phred+33 (modern Illumina caps at Q41) yet valid under Phred+64 ⇒ inferred Phred+64. The closed overlap [64,74] is the genuinely ambiguous band, matching qcfail.com and Cock et al.
- **LIMITATIONS §1 residual is faithful and irreducibly correct**: only a FASTQ where *every* read is confined to [64,74] remains undeterminable; no method can resolve it (both encodings decode it without error). The fix narrows the residual to exactly this case.

### Stage B — Implementation realises the description (verdict: PASS)

Code reviewed: `QualityScoreAnalyzer.cs:308–354` (the overload) + enum `EncodingConfidence` (30–51), constants `Phred64Offset=64`, `Phred33PlausibleMaxChar=74`.

Decision logic, hand-traced against the ASCII ranges:
- `min < 64` ⇒ Phred33 / **Definitive** (proof). ✓
- else `max > 74` ⇒ Phred64 / **Inferred**. ✓
- else (all chars in [64,74]) ⇒ Phred33 / **Ambiguous**. ✓
- empty / all-empty / all-null input ⇒ count 0 ⇒ Phred33 / Ambiguous, span (0,0). ✓ null enumerable ⇒ `ArgumentNullException`. ✓

Independent cross-checks performed this session:
- **Branch table** recomputed by hand for every test input; all match (Definitive at 63, Ambiguous at 64, Ambiguous at 74, Inferred at 75; low-char-beats-high-char [35,80,104] ⇒ Definitive Phred+33; across-reads resolution in both directions).
- **Across-reads disambiguation is real**: a lone overlap-only read (`@F`/`Ascii(64,70,73)`) returns per-read Phred33 by *default*, not proof; adding a read with a sub-64 char flips the file result to Definitive Phred+33 (test `OverlapOnlyRead_ResolvedByAnotherReadWithLowChar`), and adding a >74-char read flips it to Inferred Phred+64.
- **Single-read parity — proven, not sampled.** I exhaustively scanned all (min,max) ASCII extreme combinations (33..126) for a single read: the file-level overload's chosen `Encoding` equals `DetectEncoding(string)` in **all** cases (0 mismatches), despite the per-string detector using a `< 59` low threshold and the file-level one using `< 64` — they coincide on the final encoding because the 59–63 band defaults to Phred33 in the per-string detector and is Definitive Phred33 in the file-level one. The 6 parity TestCases are representative of this proven invariant.

Test-quality audit (`QualityScoreAnalyzer_DetectEncodingFile_Tests.cs`, 19 tests):
- Expectations are **sourced** (Cock et al. ASCII boundaries, 'J'=74 ceiling), built from raw ASCII codes via an `Ascii(...)` helper so the byte ranges under test are explicit — not code echoes.
- **No green-washing**: every assertion is exact `Is.EqualTo` on Encoding/Confidence/MinAscii/MaxAscii/CharactersExamined; no ranges, `Contains`, or weakened matchers.
- **All branches covered**: Definitive (63, across-reads low char, low-beats-high), Inferred (75, across-reads high char), Ambiguous (overlap-only, 64, 74), empty input, null/empty element skipping, all-empty, null enumerable throw, single-read parity, global span/count.

### Re-validation findings

- **No defect found.** Description, code, and tests are mutually consistent and faithful to the external sources. The asymmetry (proof of Phred+33 vs inference of Phred+64) and the [64,74] overlap boundary are information-theoretically correct. The single-read parity claim is provably true for all inputs. Full suite green.
- **State: CLEAN.** Stage A PASS, Stage B PASS. No code/test changes were needed; this re-validation edited only validation docs.
- Full unfiltered suite: `Failed: 0, Passed: 18880` (Seqeron.Genomics.Tests) plus all other projects green (`dotnet test Seqeron.sln -c Debug`, 2026-06-26).


## Runtime enforcement (LimitationPolicy)

This unit's guarded branch — encoding auto-detection on overlap-confined input (every quality char in ASCII 64–74) — has **minimum access mode `Permissive`** (`Seqeron.Genomics.Core.LimitationCatalog`). Under the default `LimitationPolicy.DefaultMode = Moderate` it throws `SeqeronLimitationException` (this guarded branch is allowed only under `Permissive`); see [LIMITATIONS.md](../LIMITATIONS.md) › Runtime enforcement. Additive policy layer; the validated contract and `✅ CLEAN` verdict are unchanged.
