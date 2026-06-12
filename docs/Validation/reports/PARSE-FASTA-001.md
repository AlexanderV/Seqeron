# Validation Report: PARSE-FASTA-001 — FASTA format parser

- **Validated:** 2026-06-12   **Area:** FileIO
- **Canonical method(s):** `FastaParser.Parse`, `FastaParser.ParseFile`, `FastaParser.ParseFileAsync`, `FastaParser.ToFasta`, `FastaParser.WriteFile` (in `src/Seqeron/Algorithms/Seqeron.Genomics.IO/FastaParser.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "FASTA format"** (https://en.wikipedia.org/wiki/FASTA_format):
  - "A sequence begins with a greater-than character (`>`) followed by a description of the sequence (all in a single line)."
  - "The lines immediately following the description line are the sequence representation" — multi-line; multi-FASTA is "obtained by concatenating several single-sequence FASTA files."
  - Line length "typically no more than 80 characters in length."
  - "lower-case letters are accepted and are mapped into upper-case."
  - "Anything other than a valid character would be ignored (including spaces, tabulators, asterisks, etc...)."
  - Legacy `;` comment convention: a leading `;` line "was taken as a comment ... ignored by software" — now **deprecated**.
- **NCBI BLAST Help — Query Input** (https://blast.ncbi.nlm.nih.gov/doc/blast-topics/):
  - "The description line (defline) is distinguished from the sequence data by a greater-than (`>`) symbol at the beginning."
  - "It is recommended that all lines of text be shorter than 80 characters in length."
  - "Blank lines are not allowed in the middle of FASTA input."
  - "lower-case letters are accepted and are mapped into upper-case."
- **NCBI GenBank FASTA** (https://www.ncbi.nlm.nih.gov/genbank/fastaformat/): defline must begin with `>` followed by a unique SeqID "that should not contain any spaces."
- **scikit-bio FASTA format docs** (reference implementation, confirms ID/description split): "all characters after the `>` and before the first whitespace character (if any) are taken as the sequence ID. If a description is present, it is taken as the remaining characters that follow the sequence ID and initial whitespace(s)."
- **Lipman & Pearson (1985)**, Science 227:1435–41 — origin of the FASTA tool/format (historical anchor; format detail comes from the specs above).

### Format rules (validated)
1. Record starts with `>` + a single-line header (defline).
2. **ID = first whitespace-delimited token after `>`; description = remainder after the first whitespace** (space or tab). Confirmed by scikit-bio reference wording.
3. Sequence on following lines until next `>` or EOF; multi-line sequence is concatenated.
4. Multiple records per file (multi-FASTA = concatenation of single records).
5. Whitespace and invalid characters within sequence lines are ignored.
6. Lowercase letters mapped to uppercase.
7. Line length conventionally ≤ 80 chars on output.
8. Blank lines "not allowed" by spec; tolerant parsers skip them.
9. `;` comment convention is legacy/deprecated.

### Edge-case semantics
- Empty input → no records. Record with no sequence → invalid (entry = header + sequence). Sequence with no header → not a valid record. CRLF/LF both must be tolerated. Trailing newline ignored. `>` appearing inside a header (not at line start) stays in the description.

### Findings / divergences
None material. The Evidence/TestSpec docs accurately capture the sourced rules.

## Stage B — Implementation

### Code path reviewed
`FastaParser.cs`: `ParseReader` (lines 107–139, shared by `Parse`/`ParseFile`), `ParseFileAsync` (44–77, structurally identical), `CreateEntry` (141–149), `ToFasta` (82–97). `FastaEntry` (155–185). Normalization to uppercase happens in `DnaSequence` ctor (`DnaSequence.cs:30`, `ToUpperInvariant()`).

### Rules realised correctly? (evidence)
- **Header detection** — `line.StartsWith('>')` begins a record (115). Matches rule 1.
- **ID/description split** — `header.Split({' ','\t'}, 2)`: `parts[0]` = ID, `parts[1]` = description, else `null` (144–146). Exactly matches the validated "first whitespace separates ID from description" rule (space **or** tab).
- **Multi-line concatenation** — non-header lines append non-whitespace chars to `sequenceBuilder` (127–131); record emitted on next `>` or EOF (117–119, 135–137). Matches rules 3, 5.
- **Multi-record** — loop yields each record. Matches rule 4.
- **Whitespace/invalid in sequence** — `if (!char.IsWhiteSpace(c))` strips all internal/leading/trailing whitespace (127–131). Matches rule 5.
- **Lowercase** — normalized to uppercase by `DnaSequence`. Matches rule 6.
- **CRLF tolerance** — `TextReader.ReadLine()` strips both `\n` and `\r\n`; additionally the per-char whitespace filter would drop any residual `\r`. Verified behaviorally: no `\r` leaks into the sequence.
- **Blank lines** — a blank line is a non-header line contributing zero chars → effectively skipped. Matches rule 8.
- **ToFasta** — emits `>` + Header, wraps sequence at `lineWidth` (default 80). Matches rules 1, 7.

### Worked example (recomputed against the actual code)
Input (CRLF): `">seq1 desc here\r\nACGT\r\nTTGG\r\n>seq2\r\nGGCC\r\n"`
- Record 1: ID=`seq1`, Description=`desc here`, Sequence=`ACGTTTGG` (no `\r`).
- Record 2: ID=`seq2`, Description=`null`, Sequence=`GGCC`.
- Count = 2. Confirmed by executing the parser.

Additional traced edge cases (executed):
- `>id1 a > b description\nACGT` → ID=`id1`, Description=`a > b description` (mid-header `>` preserved in description).
- Leading sequence with no header (`ACGT\nGGCC\n>real\nTTTT`) → only `real` yielded (orphan sequence silently dropped).
- Header without sequence (`>empty\n>has\nATGC`) → only `has` yielded.
- Trailing newline → no extra/empty record.

### Variant/delegate consistency
`Parse`, `ParseFile`, and `ParseFileAsync` share identical parsing logic (sync vs async copy). `ToFasta`/`WriteFile` agree (WriteFile = `File.WriteAllText(ToFasta(...))`). Round-trip parse→format→parse preserves data (tested).

### Test quality audit
`FastaParserTests.cs` — now 35 tests (was 34 in this project; suite was 4485). Assertions check exact sourced values (`Is.EqualTo`), cover all Stage-A edge cases: single/multi/multiline, ID-only & ID+description, tab separator, special chars (`gi|...|`), whitespace stripping, blank-line skip, header-without-sequence, lowercase normalization, line-width wrapping, round-trip, file I/O, async. Added **`Parse_CrlfLineEndings_NoCarriageReturnInSequence`** to lock CRLF tolerance (previously untested, a flagged risk).

### Findings / defects
No defects. Minor non-defect notes (documented, sourced behavior — not bugs):
- Orphan sequence lines (no preceding header) and headers with no sequence are silently dropped rather than raising an error. This matches the Evidence doc's documented intent ("entry = header + sequence; header without sequence is not yielded"). A strict parser could surface a diagnostic, but the lenient behavior is a defensible, sourced choice.
- Legacy `;` comment lines are treated as sequence lines (their non-whitespace chars would be fed to `DnaSequence` and could throw on invalid chars). This convention is deprecated and out of scope per the TestSpec; modern FASTA does not use it.

## Verdict & follow-ups
- **Stage A: PASS** — format rules and ID/description split confirmed against Wikipedia, NCBI BLAST, NCBI GenBank, and the scikit-bio reference implementation.
- **Stage B: PASS** — implementation faithfully realises every validated rule; CRLF is tolerated with no `\r` leakage; worked examples recomputed against the code match.
- **State: CLEAN** — no defect; added one CRLF regression test. Build green, full suite 4486 passed / 0 failed.
- No logged defects.
