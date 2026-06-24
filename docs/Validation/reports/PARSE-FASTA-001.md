# Validation Report: PARSE-FASTA-001 — FASTA format parser

- **Validated:** 2026-06-24   **Area:** FileIO
- **Canonical method(s):** `FastaParser.Parse`, `FastaParser.ParseFile`, `FastaParser.ParseFileAsync`, `FastaParser.ToFasta`, `FastaParser.WriteFile` (in `src/Seqeron/Algorithms/Seqeron.Genomics.IO/FastaParser.cs`); record type `FastaEntry`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "FASTA format"** (https://en.wikipedia.org/wiki/FASTA_format, fetched 2026-06-24):
  - "A sequence begins with a greater-than character (`>`) followed by a description of the sequence (all in a single line)."
  - Sequence is "interleaved, or on multiple lines ... but may also be sequential, or on a single line." Multi-FASTA "would be obtained by concatenating several single-sequence FASTA files in one file"; each new sequence must start with `>`.
  - "Anything other than a valid character would be ignored (including spaces, tabulators, asterisks, etc...)."
  - "lower-case letters are accepted and are mapped into upper-case."
  - Line length "typically no more than 80 characters."
- **Biopython SeqIO FASTA** (https://biopython.org/wiki/SeqIO and FastaIO docs, fetched 2026-06-24): with the default title parser, "the entire title line will be used as the **description**, and the **first word as the id** and name." i.e. the title (after `>`) is split at the **first whitespace** → id = first word, description = remainder. This is the canonical reference-implementation convention.
- **NCBI FASTA / BLAST** (per Evidence doc): defline begins with `>`; "Blank lines are not allowed in the middle of FASTA input"; lowercase mapped to uppercase.

### Format rules (validated)
1. Record starts with `>` + a single-line header (defline).
2. **id = first whitespace-delimited token after `>`; description = remainder** (space or tab). Matches Biopython default.
3. Sequence on following lines until next `>` or EOF; multi-line sequence concatenated.
4. Multiple records per file (multi-FASTA = concatenation of single records).
5. Whitespace within sequence lines ignored.
6. Lowercase letters mapped to uppercase.
7. Output line length conventionally ≤ 80.
8. Blank lines disallowed by spec; tolerant parsers skip them.

### Edge-case semantics
Empty input → no records. Header with no following sequence → not a valid record (entry = header + sequence). CRLF/LF both tolerated. Trailing newline ignored.

### Independent cross-check (hand trace)
Input `>seq1 First desc\nAAAA\nCCCC\n>seq2 Second\nGGGGTTTT`:
- Record 1: id=`seq1`, description=`First desc`, sequence=`AAAACCCC` (two lines concatenated).
- Record 2: id=`seq2`, description=`Second`, sequence=`GGGGTTTT`.
- Record count = 2. This matches Biopython's `SeqIO.parse(..., "fasta")` semantics (id = first word, description = remainder, concatenated sequence) and the code (verified below).

### Findings / divergences
None at the description level. TestSpec/Evidence accurately capture the sourced rules.

## Stage B — Implementation

### Code path reviewed
`FastaParser.cs`: `ParseReader` (107–139, shared by `Parse`/`ParseFile`); `ParseFileAsync` (44–77, structurally identical, async); `CreateEntry` (141–149); `ToFasta` (82–97); `WriteFile` (102–105); `FastaEntry` (155–185). Uppercase normalisation occurs in `DnaSequence` ctor (`DnaSequence.cs:30`, `ToUpperInvariant()`); character validation `DnaSequence.cs:112` (`ValidateSequence`).

### Rules realised correctly? (evidence)
- **Header detection** — `line.StartsWith('>')` begins a record (115). Rule 1. ✓
- **id/description split** — header is `Substring(1).Trim()` then `Split({' ','\t'}, 2)`: `parts[0]` = id, `parts[1]` = description, else `null` (122, 144–146). Matches "first whitespace separates id from description" (space **or** tab). ✓
- **Multi-line concatenation** — non-header lines append non-whitespace chars; record emitted on next `>` or EOF (117–137). Rules 3, 5. ✓
- **Multi-record** — iterator yields each record. Rule 4. ✓
- **Whitespace stripped** — `if (!char.IsWhiteSpace(c))` drops internal/leading/trailing whitespace incl. residual `\r` (127–131). Rule 5. ✓
- **Lowercase → uppercase** — by `DnaSequence` ctor. Rule 6. ✓
- **Blank lines** — a blank non-header line contributes 0 chars → effectively skipped. Rule 8. ✓
- **Header w/o sequence** — guarded by `sequenceBuilder.Length > 0` (117, 135); empty-sequence headers not yielded. ✓
- **CRLF** — `ReadLine()` strips `\r\n`; per-char filter drops any stray `\r`. ✓
- **ToFasta** — emits `>`+Header, wraps at `lineWidth` (default 80) (87–94). Rules 1, 7. ✓

### Cross-verification table (recomputed vs code, via test run)
| Case | Input | Expected | Code result |
|------|-------|----------|-------------|
| id/desc split | `>NM_001 Homo sapiens gene X (GeneX), mRNA` | id=`NM_001`, desc=`Homo sapiens gene X (GeneX), mRNA` | ✓ (M5) |
| NCBI pipes | `>gi\|12345\|gb\|AAA00000.1\| hypothetical protein` | id=`gi\|12345\|gb\|AAA00000.1\|`, desc=`hypothetical protein` | ✓ (M12) |
| multi-line concat | `>g\nAAAA\nCCCC\nGGGG\nTTTT` | `AAAACCCCGGGGTTTT` | ✓ (M3) |
| multi-record | 3 records | count=3, ids seq1/seq2/seq3 | ✓ (M2) |
| whitespace in seq | `AT GC\nGG\tCC` | `ATGCGGCC` | ✓ (M7) |
| lowercase | `acgtacgt` | `ACGTACGT` | ✓ |
| CRLF, 2 recs | CRLF input | no `\r` in seq, count=2 | ✓ |
| header no seq | `>empty\n>has\nATGC` | count=1, id=`has` | ✓ |

### Variant/delegate consistency
`Parse` (StringReader) and `ParseFile` (StreamReader) both delegate to `ParseReader`; `ParseFileAsync` duplicates the same logic over `ReadLineAsync` (verified line-by-line equivalent). `WriteFile` = `File.WriteAllText(ToFasta(...))`. Consistent.

### Test quality audit
`FastaParserTests.cs` — 29 tests with **exact** `Is.EqualTo` assertions on id, description, concatenated sequence, record count, wrap widths, and round-trips; covers every Stage-A edge case (empty, whitespace-only, blank lines, CRLF, header-without-sequence, tab delimiter, special chars, lowercase, line-width boundaries). `Properties/FastaRoundTripProperties.cs` — 5 property tests. Total 34, all passing (`--filter ~FastaParser`: 34 passed / 0 failed). Tests are real and deterministic.

### Findings / notes
1. **Scope: DNA-only.** `FastaParser` constructs `DnaSequence`, whose `ValidateSequence` (`DnaSequence.cs:112`) accepts **only A/C/G/T** and throws `ArgumentException` on anything else. Therefore any otherwise-valid FASTA containing IUPAC ambiguity codes (N, R, Y, …), RNA (U), gaps (`-`), `*`, or protein residues will throw — even though Wikipedia says invalid characters should be *ignored* and IUPAC codes are listed in the Evidence doc's character table. This is a **scope boundary of the unit** (it is a DNA FASTA reader, and the TestSpec/Evidence frame it as such), not a defect in the FASTA-parsing mechanics. The id/description/concatenation/multi-record/whitespace logic — the actual subject of this unit — is correct. Drives the PASS-WITH-NOTES.
2. **Description leading whitespace (cosmetic).** Because the header is `Trim()`-ed and then split with limit 2 on a *single* delimiter char, runs of whitespace between id and description leave the extra whitespace at the start of the description (e.g. `>seq1  desc` → description=` desc`). Biopython lstrips this. Cosmetic only; no test exercises double-spacing, and round-trip is unaffected (single space re-emitted). Not a defect.

## Verdict & follow-ups
- **Stage A: PASS.** Format model (id=first word, description=remainder, multi-line concat, multi-record, whitespace stripped, lowercase→upper, 80-col output) is correct and matches Wikipedia + the Biopython reference convention.
- **Stage B: PASS-WITH-NOTES.** Code faithfully realises the validated FASTA model; all 34 tests pass with exact-value assertions. Notes: (1) DNA-only scope — non-ACGT FASTA throws via `DnaSequence` validation (scope boundary, not a parse bug); (2) cosmetic leading-whitespace retention in description on multi-space headers.
- **End state: CLEAN.** No code change required; no half-fix. The two notes are scope/cosmetic, not defects in the parsing unit.
