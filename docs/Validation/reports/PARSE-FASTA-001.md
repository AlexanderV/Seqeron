# Validation Report: PARSE-FASTA-001 — FASTA format parser

- **Validated:** 2026-06-25 (re-validated; original 2026-06-24)   **Area:** FileIO
- **Canonical method(s):** `FastaParser.Parse`, `FastaParser.ParseFile`, `FastaParser.ParseFileAsync`, `FastaParser.ToFasta`, `FastaParser.WriteFile`, plus the `SequenceAlphabet` overloads (in `src/Seqeron/Algorithms/Seqeron.Genomics.IO/FastaParser.cs`); record types `FastaEntry`, `FastaRecord`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS  _(see "Re-validation 2026-06-25" below; supersedes the earlier PASS-WITH-NOTES once the DNA-only scope note was addressed by the opt-in alphabets)_

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

## Re-validation 2026-06-25 — fresh Stage A/B (post opt-in multi-alphabet mode)

Independent re-validation of the extended unit (canonical parse surface + `SequenceAlphabet`
overloads). External sources retrieved THIS session; repo TestSpec/Evidence/tests not trusted.

- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

### Sources retrieved this session
- **Wikipedia "FASTA format"** (fetched 2026-06-25): record begins with `>` + single-line header;
  multi-FASTA = concatenated records each starting `>`; "lower-case letters are accepted and are
  mapped into upper-case"; "anything other than a valid character would be ignored (including
  spaces, tabulators, asterisks…)"; lines "typically no more than 80 characters". First-word-as-id
  is the Biopython default title convention.
- **Wikipedia "Nucleic acid notation"** (NC-IUB 1985; fetched 2026-06-25): IUPAC nucleotide set =
  A C G T U + W S M K R Y B D H V N, plus gap `-`. Confirmed letter-for-letter.
- **bioinformatics.org SMS2 IUPAC** (fetched 2026-06-25): amino-acid one-letter codes = 20 standard
  (A R N D C Q E G H I L K M F P S T W Y V) + B (Asx) Z (Glx) J (Xle) X (Xaa) + U (Sec) O (Pyl) +
  `*` (termination); `-`/`.` listed as gap.

### Stage A — alphabet membership confirmed vs source
| Alphabet | Code set in source | `FastaParser` set | Match |
|----------|--------------------|--------------------|-------|
| StrictDna | A C G T | A C G T | ✓ |
| Rna | A C G U | A C G U | ✓ |
| IupacNucleotide | A C G T U R Y S W K M B D H V N + `-` | identical | ✓ |
| Protein | 20 + B Z J X + U O + `*` | identical | ✓ |
Note: Protein set excludes the gap `-`/`.`; defensible for a strict-residue alphabet (no defect).

### Stage B — independent cross-check (run against built library, 2026-06-25)
Parsed small inputs directly through the compiled overloads:
| Mode | Input | Parsed result | Verdict |
|------|-------|---------------|---------|
| default | `>seq1 First desc\r\nAA aa\r\n\r\nCCCC\r\n>seq2\tTab Desc\nGGGGtttt` | rec1 id=`seq1` desc=`First desc` seq=`AAAACCCC`; rec2 id=`seq2` desc=`Tab Desc` seq=`GGGGTTTT`; count=2 | CRLF stripped, blank line skipped, intra-line ws ignored, lc→UC, space+tab header split — byte-identical to pre-change ✓ |
| default strict | `AUGC` / `ACGTNRY` / `MWYK` | all throw `ArgumentException` (U/N/M rejected) | ✓ |
| Rna | `>rna1 messenger\naugcAUGC\nGGUU` | id=`rna1` desc=`messenger` seq=`AUGCAUGCGGUU` (U kept) | ✓; T rejected ✓ |
| Protein | `>prot1 sample\nARNDCQEGHILKMFPSTWYV\nBZXUOJ*` | seq=`ARNDCQEGHILKMFPSTWYVBZXUOJ*` (27 chars, all preserved) | ✓; `@` rejected ✓ |
| IupacNucleotide | `>amb degenerate\nACGTNRYSWKMBDHV-U` | seq=`ACGTNRYSWKMBDHV-U` (ambiguity+gap+U) | ✓; `E` rejected ✓ |

### Stage B — test-quality audit
Full unfiltered `dotnet test Seqeron.sln -c Debug` = **0 failed** (Seqeron.Genomics.Tests: 18819
passed). FASTA coverage: `FastaParserTests.cs` (default path: multi-record, multi-line concat,
empty/whitespace-only, header±desc, NCBI pipes, tab delimiter, blank lines, CRLF, header-without-
sequence, lowercase, line-width boundaries, round-trip, file/async I/O); `FastaParser_Alphabet_Tests.cs`
(14 tests: RNA-U preserved, protein W/Y/`*`/B/Z/J/X/U/O, IUPAC ambiguity+gap, per-mode reject, and
the StrictDna/default regressions); `FastaParser_MutationKillers_Tests.cs` (async flush guards +
exact wrap output). Assertions are exact `Is.EqualTo`/`Should().Be` on source-traced values, not
code echoes; every public method/overload and every Stage-A edge case is covered.

### Re-validation verdict
**Stage A PASS / Stage B PASS / End-state CLEAN.** Alphabets trace letter-for-letter to NC-IUB 1985
and the IUPAC amino-acid codes; record structure traces to the FASTA spec; the default DNA-only path
is confirmed byte-identical to the pre-change behaviour. No defect; no code change required.

---

## Update 2026-06-24 — opt-in multi-alphabet mode (limitation fix)

Note (1) above (DNA-only scope) is now addressed by an **opt-in** API, while the default path is kept strict DNA-only and byte-for-byte unchanged.

### Alphabets retrieved (this session, verbatim from citable sources)
- **IUPAC nucleotide** — A C G T U R Y S W K M B D H V N + gap `-`. Source: NC-IUB (1985), "Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences", *Nucleic Acids Research* 13(9):3021–3030, retrieved via Wikipedia "Nucleic acid notation" (https://en.wikipedia.org/wiki/Nucleic_acid_notation, fetched 2026-06-24); cross-confirmed at https://www.bioinformatics.org/sms/iupac.html.
- **RNA** — A C G U (bioinformatics.org IUPAC nucleotide table, fetched 2026-06-24).
- **Protein (IUPAC amino-acid one-letter codes)** — 20 standard A R N D C Q E G H I L K M F P S T W Y V + ambiguity B (Asx) Z (Glx) J (Xle) X (Xaa) + rare U (Sec) O (Pyl) + stop `*`. Source: https://www.bioinformatics.org/sms2/iupac.html (fetched 2026-06-24); cross-confirmed at NCBI BLAST topics (https://blast.ncbi.nlm.nih.gov/doc/blast-topics/, fetched 2026-06-24).

### API added
`SequenceAlphabet` enum {`StrictDna` (default), `IupacNucleotide`, `Rna`, `Protein`} and overloads `FastaParser.Parse(string, SequenceAlphabet)`, `ParseFile(string, SequenceAlphabet)`, `ParseFileAsync(string, SequenceAlphabet)` returning a new `FastaRecord` (id, description, raw uppercased sequence string, alphabet). Out-of-alphabet characters throw `ArgumentException`. The default no-alphabet overloads (returning `FastaEntry`/`DnaSequence`) are unchanged.

### Tests
`FastaParser_Alphabet_Tests.cs` (14 tests): RNA-with-U preserved (`AUGCAUGCGGUU`), protein with W/Y/`*` and ambiguity (`MWYXBZJUO*` preserved), IUPAC-ambiguous DNA preserved (`ACGTNRYSWKMBDHV-U`), per-mode out-of-alphabet rejection, and the default strict-DNA still-rejects regression (U / protein residues / N,R,Y all throw on the default `Parse`). All 48 FastaParser tests (34 existing + 14 new) pass.

- **Residual note:** multi-space-header `Description` keeps a single leading space (by-design header-split contract; unchanged).
- **Checklist:** Status reset `☑`→`☐` in the root registry pending re-validation of the extended unit.
