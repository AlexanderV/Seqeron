# Validation Report: PARSE-GFF-001 — GFF3/GTF Format Parser (reading)

- **Validated:** 2026-06-12   **Area:** FileIO
- **Canonical method(s):** `GffParser.Parse(content/reader)`, `GffParser.ParseFile(path)`; supporting `ParseLine`, `ParseAttributes`, `UnescapeGff`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.IO/GffParser.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/GffParserTests.cs` (106 tests in fixture)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

> Scope note: this unit is GFF3 **reading/parsing**. The sibling writing path (ANNOT-GFF-001)
> was validated separately; `ToGff3`/`WriteToStream`/`EscapeGff` are exercised here only as
> needed for round-trip context.

## Stage A — Description

### Sources opened & what they confirm
- **Sequence Ontology GFF3 Specification v1.26** (Lincoln Stein, 18 Aug 2020),
  https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md — fetched and
  read. Confirms every claim below.
- Wikipedia "General feature format" and UCSC Genome Browser FAQ (FAQformat) — corroborate the
  9-column layout, 1-based inclusive coordinates, strand/frame conventions, and GTF
  `key "value";` attribute syntax (cited in the Evidence doc; consistent with SO spec).

### Spec facts verified EXACTLY against GFF3 v1.26
1. **Nine TAB-separated columns**, in order: `seqid, source, type, start, end, score, strand,
   phase, attributes`. ✓
2. **Coordinates are 1-based and INCLUSIVE on both ends.** Spec: coordinates are "1-based
   integer coordinates" and "Start is always less than or equal to end." Therefore
   **feature length = end − start + 1** (contrast BED, which is 0-based half-open). ✓
3. **strand** ∈ {`+`, `-`, `.`, `?`}; `.` = non-stranded, `?` = unknown. ✓
4. **phase** ∈ {`0`, `1`, `2`} and is REQUIRED for CDS; `.` = undefined (non-CDS). ✓
5. **score** `.` = undefined. ✓
6. **Column 9**: `tag=value` pairs separated by `;`. Reserved tags are **capitalized**
   (`ID`, `Parent`, `Name`, …) and **attribute names are case-sensitive**. Reserved characters
   (`,`, `=`, `;`, `&`, plus tab/newline/CR/`%`/controls) are **percent-encoded** and must be
   **percent-DECODED on read**. Multiple values for one tag are comma-separated; **`Parent`
   may be multi-valued** (e.g. `Parent=mRNA1,mRNA2`). ✓
7. **Directives** begin `##` (notably `##gff-version 3`), **comments** begin `#`, and a
   `##FASTA` directive ends the annotation section. ✓

### Hand-constructed worked example (1-based inclusive)
Line (canonical SO gene, tab-separated):
`ctg123 . gene 1000 9000 . + . ID=gene00001;Name=EDEN`
- seqid=`ctg123`, source=`.`, type=`gene`, start=**1000**, end=**9000**, score=`.`→null,
  strand=`+`, phase=`.`→null, attributes `{ID=gene00001, Name=EDEN}`.
- Length = 9000 − 1000 + 1 = **8001** bp (NOT 8000). ✓ confirms 1-based inclusive.

CDS example `... CDS 1201 1500 . + 0 ID=cds00001;Parent=mRNA00001`:
length = 1500 − 1201 + 1 = **300** bp = 100 codons. ✓ (divisible by 3, consistent with phase 0).

### Findings / divergences
None. Description matches GFF3 v1.26 exactly.

## Stage B — Implementation

### Code path reviewed
`GffParser.cs`: `Parse` (lines 84–117) → `ParseLine` (119–152) → `ParseAttributes`
(154–197) → `UnescapeGff` (199–202).

### Conformance checks
- **9-column TAB split**: `line.Split('\t')` (line 121). Columns 0–7 mapped to
  seqid/source/type/start/end/score/strand/phase; column 8 (when present) → attributes. ✓
- **1-based INCLUSIVE coordinates preserved verbatim**: `start`/`end` parsed as raw integers
  (lines 129–132); no off-by-one applied at parse time. Length recovered as `End − Start + 1`.
  `ExtractSequence` correctly converts to 0-based half-open (`Start − 1` .. `End`, line 495–496)
  only when slicing, which is the right place for the conversion. ✓
- **score**: `.` → null; otherwise invariant-culture float parse (lines 134–139). ✓
- **strand**: first char of column 7; `+`/`-`/`.`/`?` all retained (line 141). ✓
- **phase**: `.` → null; else int parse 0/1/2 (lines 143–145). ✓
- **Attributes (GFF3)**: split on `;`, then split each on the **first** `=` via
  `IndexOf('=')` and slice (lines 183–193). Splitting on first `=` means an `=` inside an
  (encoded) value does not break parsing. Key and value are each individually
  percent-DECODED **after** structural splitting — correct ordering, so an encoded `%3B`/`%3D`
  in a value is never mistaken for a separator. ✓
- **Percent-DECODING**: `Uri.UnescapeDataString` (line 201) implements RFC 3986 single-pass
  decode (no double-decode): `%3B`→`;`, `%253B`→`%3B`. ✓
- **Case sensitivity**: GFF3 attribute dictionary uses `StringComparer.Ordinal`
  (case-sensitive per spec); GTF/GFF2 uses `OrdinalIgnoreCase` (line 156–158). ✓
- **Multi-valued Parent**: kept as the raw comma-joined string on parse
  (`Parent=mRNA1,mRNA2`); `BuildGeneModels` splits on `,` to assign children to all parents
  (lines 285–294). ✓
- **Directives / comments / blank lines**: `##gff-version` sets format (3 → GFF3, else GFF2,
  lines 98–105); other `##` directives skipped; `#` comments skipped; blank lines skipped
  (lines 92–111). ✓
- **Version mis-detection guard**: checks the digit after `##gff-version` is `'3'`, so
  `##gff-version 2.3` → GFF2 (not GFF3). ✓

### Cross-verification recomputed vs code (tests trace exact sourced values)
| Case | Input | Expected (spec) | Test |
|------|-------|-----------------|------|
| 1-based inclusive length | `gene 1 100` | len = 100 − 1 + 1 = 100; `exon 50 75` → 26 | `Parse_1BasedCoordinates_Validated` asserts `End−Start+1` |
| Phase | CDS 0/1/2, gene `.` | 0,1,2,null | `Parse_Phase_ParsedCorrectly` |
| Strand | `+ - . ?` | all retained | `Parse_Strand_AllValidValues` |
| Score `.` | undefined | null | `Parse_NoScore_ScoreIsNull` |
| Percent decode | `Name=Test%3BGene` | `Test;Gene` | `Parse_SpecialCharacters_Unescaped` |
| No double-decode | `Name=Test%253BGene` | `Test%3BGene` | `Parse_DoubleEncodedPercent_DecodedCorrectly` |
| Multi-Parent | `Parent=mRNA1,mRNA2` | raw `mRNA1,mRNA2`, split in gene models | `Parse_MultipleParentValues`, `BuildGeneModels_MultiParentExons...` |
| Case-sensitive attrs | `ID` ≠ `id` | only `ID` present | `Parse_AttributeCaseSensitive` |
| Comments/blank/directive | `#`/blank/`##` | skipped | `Parse_SkipsComments`, `Parse_SkipsEmptyLines`, `Parse_DetectsGFF3Version` |
| GFF2 not GFF3 | `##gff-version 2.3` | GFF2 (case-insensitive attrs) | `Parse_VersionDetection_DoesNotMisdetectGFF2AsGFF3` |
| Malformed | `< 8` fields | line skipped | `Parse_MalformedLine_Skips` |

### Test quality audit
Assertions check exact sourced values (start/end integers, `End−Start+1` lengths, decoded
strings, null score/phase), are deterministic, and cover all Stage-A edge cases. The earlier
weak `GreaterThanOrEqualTo`/palindrome assertions noted in the TestSpec are already replaced
with exact `EqualTo` checks.

### Findings / minor notes (not defects)
- **Lenient column count**: `ParseLine` rejects only `< 8` fields (line 122); a valid GFF3
  feature line has 9 columns, but a line with exactly 8 (no attributes column) is parsed with
  an empty attribute set. This is deliberately tolerant of real-world files that omit the
  trailing column and does not corrupt any field — it accepts a slightly malformed line rather
  than producing a wrong value. Not a correctness defect.
- **Strand stored as the first char**: a multi-character strand field would keep only `[0]`.
  In practice strand is always a single character; the valid set is fully covered.

## Verdict & follow-ups
- Stage A: **PASS** — description matches GFF3 v1.26 exactly (9 columns; 1-based inclusive,
  length = end − start + 1; strand/phase/score conventions; attribute percent-decoding;
  multi-valued Parent; directives/comments).
- Stage B: **PASS** — parser realises the validated spec; 1-based inclusive coordinates kept
  verbatim, attributes split structurally before single-pass percent-DECODE, case sensitivity
  per format. Worked examples recomputed against the code.
- **End state: CLEAN** — no defect found; no code change required.
- Build + tests: `Seqeron.Genomics.Tests` builds clean; `~Gff` filter = 106 passed;
  full suite = **4486 passed, 0 failed**.
